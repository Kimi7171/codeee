
## =========
##   MCEM
## =========

if (!requireNamespace("MASS", quietly = TRUE))
  stop("請先安裝 MASS 以使用 ginv()。")

# -------- 小工具：穩定 logistic、與安全反矩陣 --------
.logistic <- function(z, clamp = 30) {
  z <- pmin(pmax(z, -clamp), clamp)
  1/(1 + exp(-z))
}
.safe_inv <- function(A) tryCatch(solve(A), error = function(e) MASS::ginv(A))

# ========== Sobol 檢查（無套件就退回 rnorm） ==========
.use_sobol <- function() requireNamespace("randtoolbox", quietly = TRUE)

# --------- Cholesky（PSD 守備）---------
.chol_psd <- function(S, jitter = 1e-8) {
  S <- (S + t(S)) / 2
  tryCatch(chol(S), error = function(e) {
    eg  <- eigen(S, symmetric = TRUE)
    lam <- pmax(eg$values, jitter)
    chol(eg$vectors %*% diag(lam) %*% t(eg$vectors))
  })
}

# --------- 估 Σ_X（full-cov；可從外部傳入 SigmaX）---------
.est_SigmaX_full <- function(W, sigma2_mat, jitter = 1e-6) {
  SW   <- stats::cov(W)
  Dbar <- diag(colMeans(sigma2_mat))
  S    <- (SW - Dbar + t(SW - Dbar)) / 2
  eg   <- eigen(S, symmetric = TRUE)
  eg$values[eg$values < jitter] <- jitter
  eg$vectors %*% diag(eg$values) %*% t(eg$vectors)
}

# --------- 預算後驗常態參數（m_i, L_i）---------
.precompute_post_gauss <- function(SigmaX, W, sigma2_mat) {
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

# --------- 一次抽好：每個 i 產生最多 Bmax 個樣本（支援 Sobol/反對稱）---------
.pre_sample_Xi <- function(post, Bmax, use_sobol = TRUE, antithetic = TRUE, seed = NULL) {
  n <- nrow(post$mi); p <- ncol(post$mi)
  if (!is.null(seed)) set.seed(seed)
  if (antithetic && (Bmax %% 2 == 1L)) Bmax <- Bmax + 1L  # 反對稱需要偶數
  Xi_arr <- array(NA_real_, dim = c(n, Bmax, p))          # [i, b, j]
  make_Z <- function(m, d) {
    if (use_sobol && .use_sobol()) {
      U <- randtoolbox::sobol(n = m, dim = d, scrambling = 1, normal = FALSE)
      qnorm(pmin(pmax(U, 1e-12), 1 - 1e-12))
    } else {
      matrix(rnorm(m * d), m, d)
    }
  }
  for (i in seq_len(n)) {
    Z <- make_Z(Bmax, p)
    if (antithetic) {
      half <- Bmax / 2L
      Z[(half + 1L):Bmax, ] <- -Z[seq_len(half), ]
    }
    Xi <- Z %*% t(post$Li[[i]])
    Xi <- sweep(Xi, 2, post$mi[i, ], `+`)
    Xi_arr[i, , ] <- Xi
  }
  Xi_arr
}

# --------- 從 cache 取前 B_t 個，組成堆疊，並回傳對應 w_mat ---------
.do_E_from_cache <- function(beta, Y, Xi_arr, B_t, clamp_eta = 30) {
  n <- dim(Xi_arr)[1]; Bmax <- dim(Xi_arr)[2]; p <- dim(Xi_arr)[3]
  stopifnot(B_t <= Bmax)
  X_stack <- matrix(NA_real_, n * B_t, p)
  offset <- 0L
  for (i in seq_len(n)) {
    Xi <- Xi_arr[i, seq_len(B_t), , drop = FALSE]
    Xi <- matrix(Xi, B_t, p)
    idx <- (offset + 1L):(offset + B_t)
    X_stack[idx, ] <- Xi
    offset <- offset + B_t
  }
  y_stack <- rep(Y, each = B_t)
  eta     <- as.vector(cbind(1, X_stack) %*% beta)
  mu_aug  <- .logistic(eta, clamp = clamp_eta)
  
  mu_mat <- matrix(mu_aug, nrow = n, ncol = B_t, byrow = TRUE)
  w_mat  <- ifelse(matrix(Y, n, B_t), mu_mat, 1 - mu_mat)
  w_mat  <- pmax(w_mat, 1e-20); w_mat <- w_mat / rowSums(w_mat)
  
  list(X_stack = X_stack, y_stack = y_stack, w_mat = w_mat)
}

# ===================== 主要：MCEM（一次抽好 + 逐步加大 B） =====================
hetero_em_logit <- function(
    Y, W, sigma2_mat,
    init_beta = NULL,
    B0 = 100,        # 起始 B（小）
    Bmax = 800,      # 最大 B（大）
    growth = 1.6,    # B 的成長倍率：100 -> 160 -> 256 -> ...
    B_vcov = NULL,   # 變異數用 B；預設 2*Bmax
    maxit = 300, tol = 1e-6,
    verbose = TRUE, clamp_eta = 30,
    vcov = c("sandwich", "bread"),
    seed = NULL,
    SigmaX = NULL,           # 可外部指定；否則自動用 Cov(W)-diag(E delta^2)
    use_sobol = TRUE,        # 有裝 randtoolbox 就會用
    antithetic = TRUE        # 反對稱抽樣
) {
  vcov <- match.arg(vcov)
  W <- as.matrix(W); sigma2_mat <- as.matrix(sigma2_mat)
  stopifnot(nrow(W) == length(Y), all(dim(W) == dim(sigma2_mat)))
  if (is.null(B_vcov)) B_vcov <- 2L * Bmax
  
  Y <- as.numeric(Y); if (!all(Y %in% c(0,1))) stop("Y 必須是 0/1")
  n <- nrow(W); p <- ncol(W)
  x_names <- c("(Intercept)", if (is.null(colnames(W))) paste0("W", 1:p) else colnames(W))
  
  # Σ_X
  if (is.null(SigmaX)) SigmaX <- .est_SigmaX_full(W, sigma2_mat)
  
  # 後驗參數 + 一次抽好
  post   <- .precompute_post_gauss(SigmaX, W, sigma2_mat)
  Xi_arr <- .pre_sample_Xi(post, Bmax = Bmax, use_sobol = use_sobol, antithetic = antithetic, seed = seed)
  
  # 初值
  if (is.null(init_beta)) {
    beta_old <- tryCatch(
      coef(glm.fit(x = cbind(1, W), y = Y, family = binomial())),
      error = function(e) { mu0 <- min(max(mean(Y), 1e-6), 1 - 1e-6); c(qlogis(mu0), rep(0, p)) }
    )
  } else {
    stopifnot(length(init_beta) == p + 1)
    beta_old <- init_beta
  }
  names(beta_old) <- x_names
  
  # EM
  converged <- FALSE; Q_trace <- numeric(0); beta_hat <- beta_old
  B_t <- B0
  for (it in seq_len(maxit)) {
    Eobj  <- .do_E_from_cache(beta_old, Y, Xi_arr, B_t, clamp_eta)
    Xs    <- Eobj$X_stack
    ys    <- Eobj$y_stack
    w_mat <- Eobj$w_mat
    w_vec <- as.vector(t(w_mat))   # 注意 byrow=TRUE 的排列
    
    # ---- M-step（brglm2 優先；失敗才退 glm）----
    fit <- NULL
    if (requireNamespace("brglm2", quietly = TRUE)) {
      fit <- try(
        brglm2::brglm.fit(
          x = cbind(1, Xs), y = ys, weights = w_vec,
          family = binomial("logit"), type = "AS_mixed"
        ),
        silent = TRUE
      )
    }
    if (inherits(fit, "try-error") || is.null(fit) || any(!is.finite(coef(fit)))) {
      fit <- try(glm.fit(x = cbind(1, Xs), y = ys, weights = w_vec, family = binomial()),
                 silent = TRUE)
    }
    
    beta_m <- if (inherits(fit, "try-error") || is.null(fit) || any(!is.finite(coef(fit))))
      beta_old else as.numeric(coef(fit))
    names(beta_m) <- x_names
    
    # ---- 阻尼 + 步長上限（trust-region）----
    gamma <- 0.35                             # 阻尼（0.25~0.5 都可）
    beta_new_raw <- beta_old + gamma * (beta_m - beta_old)
    
    d <- beta_new_raw - beta_old
    step_max <- 2.0                           # 每次 ||Δβ|| 不超過 2
    if (sqrt(sum(d^2)) > step_max) {
      beta_new_raw <- beta_old + d * (step_max / sqrt(sum(d^2)))
    }
    
    # 絕對夾限，避免數值爆到非常大
    beta_cap <- 25
    beta_new <- pmin(pmax(beta_new_raw, -beta_cap), beta_cap)
    
    # 監控 Q（byrow=TRUE）
    eta_vec <- as.vector(cbind(1, Xs) %*% beta_new)
    eta_mat <- matrix(pmin(pmax(eta_vec, -clamp_eta), clamp_eta), n, B_t, byrow = TRUE)
    Q_new   <- sum(w_mat * (matrix(Y, n, B_t) * eta_mat - log1p(exp(eta_mat))))
    Q_trace <- c(Q_trace, Q_new)
    
    diff_norm <- sqrt(sum((beta_new - beta_old)^2))
    if (verbose && (it %% 10 == 0 || diff_norm < tol)) {
      cat(sprintf(" iter %d | B=%d | ||Δβ||=%.3e | Q=%.6f\n", it, B_t, diff_norm, Q_new))
    }
    if (diff_norm < tol) { converged <- TRUE; beta_hat <- beta_new; break }
    beta_old <- beta_new
    
    # 成長排程（最多到 Bmax）
    B_t <- min(Bmax, as.integer(ceiling(B_t * growth)))
  }
  if (!converged) { warning("EM 未在 maxit 收斂；回傳最後一輪 β"); beta_hat <- beta_old }
  
  # 變異數：直接用預抽樣，取前 B_vcov
  B_v <- min(B_vcov, dim(Xi_arr)[2])
  Eobj_v <- .do_E_from_cache(beta_hat, Y, Xi_arr, B_v, clamp_eta)
  Xs_v    <- Eobj_v$X_stack
  ys_v    <- Eobj_v$y_stack
  w_mat_v <- Eobj_v$w_mat
  w_vec_v <- as.vector(t(w_mat_v))
  
  X_aug  <- cbind(1, Xs_v)
  mu_aug <- .logistic(as.vector(X_aug %*% beta_hat), clamp = clamp_eta)
  idx_i  <- rep(seq_len(n), each = B_v)
  
  # 後驗期望分數（依個體聚合）
  resid_aug  <- (ys_v - mu_aug)
  score_post <- matrix(0, n, p + 1)
  for (j in seq_len(p + 1)) {
    contrib <- w_vec_v * resid_aug * X_aug[, j]
    score_post[, j] <- rowsum(contrib, idx_i, reorder = FALSE)
  }
  colnames(score_post) <- x_names
  
  # bread / sandwich
  wV    <- w_vec_v * (mu_aug * (1 - mu_aug))
  bread <- crossprod(X_aug, X_aug * wV)
  Hinv  <- .safe_inv(bread)
  vcov_hat <- if (vcov == "bread") Hinv else { meat <- crossprod(score_post); Hinv %*% meat %*% Hinv }
  # --- 安全指派 dimnames（避免 dimnames 長度不符） ---
  vcov_hat <- as.matrix(vcov_hat)
  nx <- length(x_names)
  
  if (nrow(vcov_hat) != nx || ncol(vcov_hat) != nx) {
    stop(sprintf("[hetero_em_logit] 內部維度不符：vcov_hat 是 %dx%d，但預期是 %dx%d (p+1)。",
                 nrow(vcov_hat), ncol(vcov_hat), nx, nx))
  }
  dimnames(vcov_hat) <- list(x_names, x_names)
  
  
  list(
    beta = structure(beta_hat, names = x_names),
    vcov = vcov_hat,
    score_mat_post = score_post,
    converged = converged,
    Q_trace = Q_trace,
    B = B_t, B_vcov = B_v,
    SigmaX = SigmaX,                 # 給診斷程式使用
    call = match.call()
  )
}
