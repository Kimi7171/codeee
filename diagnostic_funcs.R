
## ======================
##    diag_all_1stepEM
## ======================

if (!requireNamespace("MASS", quietly = TRUE))
  stop("請先安裝套件 MASS 以使用 ginv()。")

# ───── 小工具（與主檔一致；如已存在則沿用） ───── #
.logistic <- if (exists(".logistic")) get(".logistic") else function(z, clamp = 30) { z <- pmin(pmax(z, -clamp), clamp); 1/(1 + exp(-z)) }
.safe_inv <- if (exists(".safe_inv")) get(".safe_inv") else function(A) { tryCatch(solve(A), error = function(e) MASS::ginv(A)) }

# 若主檔未先載入，提供 fallback 的預算與打分函式
if (!exists(".chol_psd")) {
  .chol_psd <- function(S, jitter = 1e-8) {
    S <- (S + t(S)) / 2
    tryCatch(chol(S), error = function(e) {
      eg  <- eigen(S, symmetric = TRUE)
      lam <- pmax(eg$values, jitter)
      chol(eg$vectors %*% diag(lam) %*% t(eg$vectors))
    })
  }
}
if (!exists(".precompute_post_gauss")) {
  .precompute_post_gauss <- function(SigmaX, W, sigma2_mat) {
    W <- as.matrix(W); sigma2_mat <- as.matrix(sigma2_mat)
    n <- nrow(W); p <- ncol(W)
    mi <- matrix(NA_real_, n, p)
    Li <- vector("list", n)
    for (i in seq_len(n)) {
      Di <- diag(sigma2_mat[i, ])
      S  <- SigmaX + Di
      Ki <- SigmaX %*% .safe_inv(S)
      Vi <- SigmaX - Ki %*% SigmaX
      Li[[i]] <- .chol_psd(Vi)
      mi[i, ] <- as.numeric(Ki %*% W[i, ])
    }
    w_names <- if (is.null(colnames(W))) paste0("W", 1:p) else colnames(W)
    list(mi = mi, Li = Li, w_names = w_names)
  }
}
if (!exists(".posterior_scores_fast")) {
  .posterior_scores_fast <- function(beta, Y, post, B = 150, clamp_eta = 30) {
    n <- length(Y); p <- ncol(post$mi)
    X_stack <- matrix(NA_real_, n * B, p)
    for (i in seq_len(n)) {
      Zi <- matrix(rnorm(B * p), B, p)
      Xi <- Zi %*% t(post$Li[[i]])
      Xi <- sweep(Xi, 2, post$mi[i, ], `+`)
      idx <- ((i - 1L) * B + 1L):(i * B)
      X_stack[idx, ] <- Xi
    }
    X_aug   <- cbind(1, X_stack)
    eta     <- as.vector(X_aug %*% beta)
    mu_aug  <- .logistic(eta, clamp = clamp_eta)
    
    mu_mat <- matrix(mu_aug, nrow = n, ncol = B, byrow = TRUE)
    Y_mat  <- matrix(Y, nrow = n, ncol = B)                   # 每欄同一個 Y 向量
    w_mat  <- ifelse(Y_mat == 1, mu_mat, 1 - mu_mat)
    w_mat  <- pmax(w_mat, 1e-20); w_mat <- w_mat / rowSums(w_mat)
    
    w_vec  <- as.vector(t(w_mat))         #  用 row-wise 展平成 [i1..i1_B, i2..i2_B, ...]
    y_stack <- rep(Y, each = B)
    resid   <- y_stack - as.vector(t(mu_mat))   # 同上 row-wise 對齊
    idx_i   <- rep(seq_len(n), each = B)        # 與 X_stack 的建構一致
    
    score_post <- matrix(0, n, p + 1)
    for (j in seq_len(p + 1)) {
      contrib <- w_vec * resid * X_aug[, j]
      score_post[, j] <- rowsum(contrib, idx_i, reorder = FALSE)
    }
    colnames(score_post) <- c("(Intercept)", post$w_names)
    score_post
  }
}

# 向量化 pseudo‑R²
.pseudoR2_vec <- function(y, pi, type = c("CS","MF","NK"), eps = 1e-12) {
  type <- match.arg(type)
  y  <- as.numeric(y)
  if (!all(y %in% c(0,1))) stop("y 必須是 0/1")
  n  <- length(y)
  if (length(pi) != n) stop("y 與 pi 長度需相同")
  pi <- pmin(pmax(pi, eps), 1 - eps)
  mu <- pmin(pmax(mean(y), eps), 1 - eps)
  ll_full <- sum(y * log(pi) + (1 - y) * log(1 - pi))
  ll_null <- sum(y * log(mu) + (1 - y) * log(1 - mu))
  if (type == "CS") {
    return(1 - exp((2 / n) * (ll_null - ll_full)))
  } else if (type == "MF") {
    return(1 - (ll_full / ll_null))
  } else {
    cs <- 1 - exp((2 / n) * (ll_null - ll_full))
    denom <- 1 - exp((2 / n) * ll_null)
    if (!is.finite(denom) || denom <= 0) return(NA_real_)
    return(cs / denom)
  }
}

# === 一階近似影響點診斷（使用 EM 物件；快速後驗分數） ===
diag_all_1stepEM <- function(Y, W, sigma2_mat,
                             G.ind = NULL,
                             cut_vec = NULL,
                             em_out = NULL,
                             em_B = 300, em_maxit = 500, em_tol = 1e-6,
                             dbg = FALSE) {
  stopifnot(is.matrix(W), is.matrix(sigma2_mat),
            all(dim(W) == dim(sigma2_mat)), length(Y) == nrow(W))
  
  # cut_vec 安全處理
  .fallback_cut <- function(cut_vec) {
    if (!is.null(cut_vec)) return(cut_vec)
    if (exists("cut_vec", envir = .GlobalEnv, inherits = FALSE))
      return(get("cut_vec", envir = .GlobalEnv))
    stop("cut_vec 未提供；請先用 make_cut_vec(...) 產生並傳入。")
  }
  cut_vec <- .fallback_cut(cut_vec)
  if (is.na(cut_vec["GCD_GSPR"])) cut_vec["GCD_GSPR"] <- 1
  
  n <- nrow(W);  p <- ncol(W); k <- p + 1
  if (is.null(G.ind)) G.ind <- seq_len(n)
  B.ind <- setdiff(seq_len(n), G.ind);  b <- length(B.ind)
  
  # (1) 全資料 EM 結果
  em_fit <- if (is.null(em_out)) stop("請提供 em_out（hetero_em_logit 的輸出）") else em_out
  full_beta <- em_fit$beta
  if (length(full_beta) < (p + 1)) {
    miss <- setdiff(colnames(W), names(full_beta))
    full_beta <- c("(Intercept)" = 0, full_beta)
    if (length(miss)) full_beta[miss] <- 0
    full_beta <- full_beta[c("(Intercept)", colnames(W))]
  }
  if (any(!is.finite(full_beta))) {
    warning("EM β 有 NA/Inf，改用 Naïve β")
    full_beta <- tryCatch(
      coef(glm(Y ~ ., family = binomial(), data = data.frame(Y, W))),
      error = function(e) { mu0 <- pmin(pmax(mean(Y), 1e-6), 1 - 1e-6); c("(Intercept)" = qlogis(mu0), setNames(rep(0, ncol(W)), colnames(W))) }
    )
  }
  
  # (2) Good-only 參考 β 與資訊矩陣
  if (length(B.ind) > 0 && length(B.ind) < n) {
    em_clean <- hetero_em_logit(
      Y = Y[G.ind],
      W = W[G.ind, , drop = FALSE],
      sigma2_mat = sigma2_mat[G.ind, , drop = FALSE],
      init_beta = full_beta,
      B0 = 120,            # 小 B 起步
      Bmax = 300,          # 參考子樣本就不用太大
      growth = 1.5,
      B_vcov = 600,
      maxit = min(300, em_maxit),
      tol = max(1e-4, em_tol),
      verbose = FALSE,
      vcov = "bread",
      SigmaX = em_fit$SigmaX,
      use_sobol = TRUE, antithetic = TRUE
    )
    beta_ref     <- em_clean$beta
    J_clean_inv  <- em_clean$vcov
  } else {
    beta_ref     <- full_beta
    Xg   <- cbind(1, W[G.ind, , drop = FALSE])
    pig  <- .logistic(as.vector(Xg %*% beta_ref), clamp = 30)
    vg   <- pig * (1 - pig)
    XtVX_g     <- crossprod(Xg, Xg * vg)
    Sig_delta_g<- diag(c(0, colMeans(sigma2_mat[G.ind, , drop = FALSE])), k, k)
    J_clean    <- XtVX_g - length(G.ind) * Sig_delta_g
    J_clean_inv<- .safe_inv(J_clean)
  }
  J_clean <- .safe_inv(J_clean_inv)
  
  # (2.5) 預算一次後驗（用好點 ΣX ）
  SigmaX_good <- .est_SigmaX_full(
    W[G.ind, , drop = FALSE],
    sigma2_mat[G.ind, , drop = FALSE]
  )
  post <- .precompute_post_gauss(SigmaX_good, W, sigma2_mat)
  
  # (3) 在 β_ref 下：用 ME-aware 的邊際機率 μ̄_i = E[logit^{-1}(β'X_i) | W_i]
  X    <- cbind(1, W)
  
  B_mu   <- 400L                    # 可 300~600
  mu_bar <- numeric(n)
  for (i in seq_len(n)) {
    Z  <- matrix(rnorm(B_mu * p), B_mu, p)
    Xi <- Z %*% t(post$Li[[i]])
    Xi <- sweep(Xi, 2, post$mi[i, ], `+`)
    mu_bar[i] <- mean(.logistic(as.vector(cbind(1, Xi) %*% beta_ref), clamp = 30))
  }
  pi_ref <- mu_bar
  v_all  <- pmax(mu_bar * (1 - mu_bar), 1e-8)
  
  # generalized / modified 帽子值
  tmp         <- X %*% J_clean_inv
  quad_xJx    <- rowSums(tmp * X)
  hii_G_all   <- pmin(pmax(v_all * quad_xJx, 1e-8), 1 - 1e-6)
  
  X_hat   <- cbind(1, W)
  Vref_s  <- sqrt(v_all)
  Zhat    <- Vref_s * X_hat
  XtVtX_hat <- crossprod(Zhat)
  M_hat_inv  <- .safe_inv(XtVtX_hat)
  quad_xhat  <- rowSums((X_hat %*% M_hat_inv) * X_hat)
  hii_M_all  <- pmin(pmax(v_all * quad_xhat, 1e-8), 1 - 1e-6)
  
  # (4) 基準 pseudo‑R²（Good 子樣本）
  R_CS_all <- .pseudoR2_vec(Y[G.ind],  pi_ref[G.ind],  type = "CS")
  R_MF_all <- .pseudoR2_vec(Y[G.ind],  pi_ref[G.ind],  type = "MF")
  R_NK_all <- .pseudoR2_vec(Y[G.ind],  pi_ref[G.ind],  type = "NK")
  
  term_vec <- quad_xJx
  
  out <- data.frame(id = seq_len(n),
                    GD = NA, MD = NA, GDF = NA, MDF = NA,
                    GDB = NA, MDB = NA, GCD_GSPR = NA,
                    mCDstar = NA, R_CS = NA, R_MF = NA, R_NK = NA)
  
  Rdiff_CS <- Rdiff_MF <- Rdiff_NK <- numeric(n)
  
  # 傳統 CookD/DFFITS 作為參考
  glm_fit <- try(glm(Y ~ ., data = data.frame(Y, W), family = binomial(),
                     control = glm.control(maxit = 200, epsilon = 1e-8)), silent = TRUE)
  if (inherits(glm_fit, "try-error") || any(!is.finite(coef(glm_fit)))) {
    if (requireNamespace("brglm2", quietly = TRUE)) {
      glm_fit <- try(brglm2::brglm(Y ~ ., data = data.frame(Y, W), family = binomial("logit"), type = "AS_mixed"), silent = TRUE)
    }
  }
  CookD_trad  <- DFFITS_trad <- rep(NA_real_, nrow(W))
  if (!inherits(glm_fit, "try-error") && is.list(glm_fit)) {
    CookD_trad  <- suppressWarnings(cooks.distance(glm_fit))
    DFFITS_trad <- suppressWarnings(dffits(glm_fit))
  }
  out <- cbind(out, CookD = CookD_trad, DFFITS = DFFITS_trad)
  
  # 先一次算好每一筆的後驗分數（B=150 可調）
  set.seed(123)  # 可重現
  scores_all <- .posterior_scores_fast(
    beta = beta_ref, Y = Y, post = post, B = 400, clamp_eta = 30
  )
  
  # (5) 逐點一階近似
  for (i in seq_len(n)) {
    vi     <- v_all[i]
    hii_G  <- hii_G_all[i]
    hii_M  <- hii_M_all[i]
    xi     <- X[i, ]
    
    # 快速版後驗分數（沿用預算 post）
    score_i <- scores_all[i, ]
    if (dbg && i == 1L) {
      score_old <- as.vector(X[i, ] * (Y[i] - pi_ref[i]))
      cat(sprintf("[diag] ||post - plug-in|| = %.3e\n", sqrt(sum((score_i - score_old)^2))))
    }
    
    den   <- if (i %in% B.ind) (1 + hii_G) else (1 - hii_G)
    d_beta <- as.vector(J_clean_inv %*% score_i) / den
    if (!(i %in% B.ind)) d_beta <- -d_beta
    
    beta_new_i <- as.vector(beta_ref + d_beta)
    idx_R      <- if (i %in% B.ind) c(G.ind, i) else setdiff(G.ind, i)
    X_clean    <- X[idx_R, , drop = FALSE]
    # 邊際化：對 idx_R 上的每個 j，算 μ̄_j(β_new_i)
    B_mu2 <- 200L
    pi_approx_i <- numeric(length(idx_R))
    for (t in seq_along(idx_R)) {
      j  <- idx_R[t]
      Z  <- matrix(rnorm(B_mu2 * p), B_mu2, p)
      Xi <- Z %*% t(post$Li[[j]])
      Xi <- sweep(Xi, 2, post$mi[j, ], `+`)
      pi_approx_i[t] <- mean(.logistic(as.vector(cbind(1, Xi) %*% beta_new_i), clamp = 30))
    }
    
    GD_i <- as.numeric(t(d_beta) %*% J_clean %*% d_beta / k)
    MD_i <- as.numeric(t(d_beta) %*% XtVtX_hat %*% d_beta / k)
    
    R_CS_i <- .pseudoR2_vec(Y[idx_R],  pi_approx_i, type = "CS")
    R_MF_i <- .pseudoR2_vec(Y[idx_R],  pi_approx_i, type = "MF")
    R_NK_i <- .pseudoR2_vec(Y[idx_R],  pi_approx_i, type = "NK")
    Rdiff_CS[i] <- abs(R_CS_all - R_CS_i)
    Rdiff_MF[i] <- abs(R_MF_all - R_MF_i)
    Rdiff_NK[i] <- abs(R_NK_all - R_NK_i)
    
    if (i %in% B.ind) {
      ri_G <- (Y[i] - pi_ref[i]) / sqrt(vi * (1 + hii_G));  hiiS_G <- hii_G / (1 + hii_G)
      ri_M <- (Y[i] - pi_ref[i]) / sqrt(vi * (1 + hii_M));  hiiS_M <- hii_M / (1 + hii_M)
    } else {
      ri_G <- (Y[i] - pi_ref[i]) / sqrt(vi * (1 - hii_G));  hiiS_G <- hii_G / (1 - hii_G)
      ri_M <- (Y[i] - pi_ref[i]) / sqrt(vi * (1 - hii_M));  hiiS_M <- hii_M / (1 - hii_M)
    }
    GDF_i <- ri_G * sqrt(hiiS_G)
    MDF_i <- ri_M * sqrt(hiiS_M)
    mCD_i <- abs(MDF_i) * sqrt((n - b - p) / p)
    GDB_i <- hiiS_G * ri_G^2
    MDB_i <- hiiS_M * ri_M^2
    
    GCD_i <- (1 / (p + 1)) * (ri_G^2) * term_vec[i]
    
    out[i, c("GD","MD","GDF","MDF","GDB","MDB","GCD_GSPR","mCDstar","R_CS","R_MF","R_NK")] <-
      c(GD_i, MD_i, GDF_i, MDF_i, GDB_i, MDB_i, GCD_i, mCD_i, Rdiff_CS[i], Rdiff_MF[i], Rdiff_NK[i])
  }
  
  # quick check：尺度訊息 + 採用之共同尺度（CS/MF/NK 都印）
  sd_safe  <- function(x) { x <- x[is.finite(x)]; if (length(x) >= 2) sd(x) else NA_real_ }
  mad_safe <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) return(NA_real_);
  stats::mad(x, center = stats::median(x), constant = 1.4826, na.rm = TRUE) }
  
  # 小工具：印出 bad/good 尺度摘要，並回傳共同尺度（good MAD → overall SD → 1e-6）
  .report_and_get_scale <- function(name, Rdiff, G.ind, B.ind) {
    sdb <- sd_safe(Rdiff[B.ind]); sdg <- sd_safe(Rdiff[G.ind])
    if (exists("summ_sep", mode = "function")) {
      print(summ_sep(abs(Rdiff), G.ind, B.ind))
    } else {
      message(sprintf("[check] R_%s scales: sd_bad=%.4g, sd_good=%.4g, ratio(good/bad)=%.2f",
                      name, sdb, sdg, ifelse(is.finite(sdb) && sdb > 0, sdg / sdb, NA_real_)))
    }
    s <- mad_safe(Rdiff[G.ind])
    if (!is.finite(s) || s <= 0) s <- sd_safe(Rdiff)
    if (!is.finite(s) || s <= 0) s <- 1e-6
    message(sprintf("[scale] R_%s common scale used = %.6g (good MAD → overall SD → 1e-6)", name, s))
    s
  }
  
  # --- CS：共同尺度 + 輸出 scale 訊息 ---
  s_pool_cs <- .report_and_get_scale("CS", Rdiff_CS, G.ind, B.ind)
  out$R_CS_raw <- abs(Rdiff_CS)
  out$R_CS     <- out$R_CS_raw / s_pool_cs
  
  # --- MF：改成共同尺度 + 輸出 scale 訊息（與 CS/NK 一致）---
  s_pool_mf <- .report_and_get_scale("MF", Rdiff_MF, G.ind, B.ind)
  out$R_MF_raw <- abs(Rdiff_MF)
  out$R_MF     <- out$R_MF_raw / s_pool_mf
  
  # --- NK：共同尺度 + 輸出 scale 訊息 ---
  s_pool_nk <- .report_and_get_scale("NK", Rdiff_NK, G.ind, B.ind)
  out$R_NK_raw <- abs(Rdiff_NK)
  out$R_NK     <- out$R_NK_raw / s_pool_nk
  
  
  # 旗標
  GDFFITS_cut  <- cut_vec["GDFFITS"]; GSDFBETA_cut <- cut_vec["GSDFBETA"]
  out <- transform(out,
                   flag_CookD    = CookD   > cut_vec["CookD"],
                   flag_DFFITS   = abs(DFFITS) > cut_vec["DFFITS"],
                   flag_GD       = GD  > cut_vec["GD"],
                   flag_MD       = MD  > cut_vec["MD"],
                   flag_GDF      = abs(GDF) > GDFFITS_cut,
                   flag_MDF      = abs(MDF) > GDFFITS_cut,
                   flag_GDB      = abs(GDB) > GSDFBETA_cut,
                   flag_MDB      = abs(MDB) > GSDFBETA_cut,
                   flag_mCDstar  = mCDstar > cut_vec["mCDstar"],
                   flag_GCD_GSPR = abs(GCD_GSPR) > cut_vec["GCD_GSPR"],
                   flag_R_CS     = R_CS > cut_vec["StdR2_CS"],
                   flag_R_MF     = R_MF > cut_vec["StdR2_MF"],
                   flag_R_NK     = R_NK > cut_vec["StdR2_NK"])
  return(out)
}
