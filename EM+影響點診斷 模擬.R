
## ==========
##   主程式
## ==========

options(error = function() {traceback(3)})

# ---- 套件 ----
library(simex)
library(brglm2)
library(dplyr)
library(data.table)
library(purrr)
library(furrr)
library(future)
library(future.apply)
library(writexl)

# ---- 載入函式檔 ----
source("C:/Users/kimi1/OneDrive/文件/論文/diagnostic_funcs.R")
source("C:/Users/kimi1/OneDrive/文件/論文/hetero_EM_logistic.R")

# ---- 平行策略（預設單執行緒；要平行再打開） ----
# plan(sequential)
# 若需要平行，改成：
plan(multisession, workers = max(1, parallel::detectCores() - 1L))


# ---------- 工具函式 ----------
# (1) 收斂旗標安全轉型
safe_converged <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(FALSE)
  isTRUE(as.logical(x))
}

# (2) Monte‑Carlo SE —— 取最後 k 步 β trace 標準差
calc_mcse <- function(beta_trace, last_k = 50L) {
  if (is.null(beta_trace) || !is.matrix(beta_trace)) return(numeric(0))
  nrow_tr <- nrow(beta_trace)
  if (nrow_tr == 0L) return(rep(0, ncol(beta_trace)))
  k <- min(last_k, nrow_tr)
  apply(beta_trace[(nrow_tr - k + 1):nrow_tr, , drop = FALSE], 2, sd)
}

# (3) AUC（ROC 面積）
AUC_fast <- function(y, p) {
  y <- as.integer(y); stopifnot(all(y %in% 0:1))
  pos <- p[y == 1]; neg <- p[y == 0]
  n1 <- length(pos); n0 <- length(neg)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# --- 用 RC 做 EM 起始值 ---
get_rc_init <- function(Y, W, sigma2_mat,
                        use_brglm2 = TRUE,
                        max_abs_beta = 50,
                        verbose = TRUE) {
  stopifnot(is.matrix(W), is.matrix(sigma2_mat),
            all(dim(W) == dim(sigma2_mat)), length(Y) == nrow(W))
  nmW <- colnames(W)
  
  # 0) 去掉近零變異欄位（避免 κ 計算不穩）
  nzv <- which(apply(W, 2, function(z) sd(z, na.rm=TRUE) < 1e-8))
  if (length(nzv)) {
    if (verbose) message("[RC-init] drop near-zero var: ", paste(nmW[nzv], collapse=", "))
    W <- W[, -nzv, drop=FALSE]
    sigma2_mat <- sigma2_mat[, -nzv, drop=FALSE]
    nmW <- colnames(W)
    if (ncol(W) == 0L) {
      mu0 <- pmin(pmax(mean(Y), 1e-6), 1 - 1e-6)
      return(list(beta = c("(Intercept)" = qlogis(mu0)), type = "SAFE"))
    }
  }
  
  df  <- data.frame(Y, W)
  frm <- as.formula(paste("Y ~", paste(nmW, collapse=" + ")))
  
  # 1) 先配一個 naive 基礎模型（prefer brglm2 → glm）
  base_fit <- NULL
  if (use_brglm2 && requireNamespace("brglm2", quietly = TRUE)) {
    base_fit <- try(brglm2::brglm(frm, family = binomial("logit"), data = df, type = "AS_mixed"), silent = TRUE)
  }
  if (inherits(base_fit, "try-error") || is.null(base_fit) || any(!is.finite(coef(base_fit)))) {
    base_fit <- try(glm(frm, family = binomial(), data = df, control = glm.control(maxit = 200, epsilon = 1e-8)), silent = TRUE)
  }
  if (inherits(base_fit, "try-error") || any(!is.finite(coef(base_fit)))) {
    if (verbose) message("[RC-init] SAFE (base fit failed)")
    mu0 <- pmin(pmax(mean(Y), 1e-6), 1 - 1e-6)
    return(list(beta = c("(Intercept)" = qlogis(mu0), setNames(rep(0, ncol(W)), nmW)), type = "SAFE"))
  }
  
  # 2) RC-like 放縮：β1_RC = β1_naive / κ̂
  b0 <- coef(base_fit)[c("(Intercept)", nmW)]
  varW <- pmax(apply(W, 2, var), 1e-12)
  ed   <- pmax(colMeans(sigma2_mat), 0)
  kappa_hat <- pmax((varW - ed) / varW, 1e-3)
  b0[nmW] <- b0[nmW] / kappa_hat[nmW]
  
  # 3) 夾限（避免極端初值）
  if (any(abs(b0) > max_abs_beta)) {
    if (verbose) message("[RC-init] clip |beta| to ", max_abs_beta)
    b0 <- pmin(pmax(b0, -max_abs_beta), max_abs_beta)
  }
  
  list(beta = b0, type = "RC_LIKE")
}

# ---------- 2. 全域設定 ----------
p_list <- c(3)           # 之後會增加 5
n_list <- c(200)         # 100之後會增加 50, 200, 500
r_target <- 0.5
kappa_set <- c(0.75, 0.85)     # 之後會增加 0.85
contam_rate_list <- c(0.1, 0.2)  # 0.1之後會增加 0.3
sigma_fac_list <- c(1.5, 2)    # 4之後會增加 7
R <- 3                   # 之後會從100次逐漸改為200次、500次，最後1000次
diag_mode_list <- c("em_step")

# 影響點門檻
make_cut_vec <- function(p, n_main, contam_n) {
  k <- p + 1
  c(
    CookD = 4 / (n_main - contam_n),
    DFFITS = 2 * sqrt(k / (n_main - contam_n)),
    GDFFITS = 3 * sqrt(k / (n_main - contam_n)),
    GSDFBETA = 9 * k / (n_main - contam_n - 3*p),
    mCDstar = 3 / sqrt((n_main - p) / (n_main - contam_n)),
    GCD_GSPR = 1, StdR2_CS = 3, StdR2_MF = 3, StdR2_NK = 3,
    GD = 4 / (n_main - p), MD = 4 / (n_main - p)
  )
}

# ---------- 3. 單回合 ----------
one_run <- function(kappa_tar, n_total, n_val,
                    p, beta_true, contam_n,
                    sigma_fac = 7, cut_vec, diag_func = diag_all_1stepEM) {
  message(sprintf("── one_run() | κ=%.2f p=%d n=%d contam=%d σ=%.2f",
                  kappa_tar, p, n_total, contam_n, sigma_fac))
  
  ## 3-1 生成 X, W, Δ²
  rho <- sqrt(r_target)
  Z   <- matrix(rnorm(n_total * (p + 1)), n_total, p + 1)
  X   <- sqrt(1 - rho^2) * Z[, 1:p] + rho * Z[, p + 1]
  colnames(X) <- paste0("X", 1:p)
  
  sigma_X2   <- apply(X, 2, var)
  bar_delta2 <- sigma_X2 * (1 / kappa_tar - 1)
  
  delta2_mat <- matrix(0, n_total, p)
  for (j in 1:p) delta2_mat[, j] <- runif(n_total, 0.5 * bar_delta2[j], 1.5 * bar_delta2[j])
  
  W <- X + matrix(rnorm(n_total * p, sd = sqrt(as.vector(delta2_mat))), n_total, p)
  colnames(W) <- paste0("W", 1:p)
  
  ## 3-2 產生 Y
  eta <- beta_true[1] + X %*% beta_true[-1]
  Y   <- rbinom(n_total, 1, plogis(eta))
  
  ## 3-3 切主樣本
  full_dat <- data.frame(Y, W)
  main_dat <- full_dat[(n_val + 1):n_total, ]
  n_main   <- n_total - n_val
  
  ## 3-5 汙染（固定偏移 k*s_eta + 遠邊界 + 只動部分座標）
  contam_idx <- sample(seq_len(n_main), contam_n)
  
  Wmat    <- as.matrix(main_dat[, paste0("W", 1:p)])
  beta_s  <- beta_true[-1]
  b0_true <- beta_true[1]
  
  eta0  <- as.numeric(b0_true + Wmat %*% beta_s)  # 汙染前 η
  s_eta <- sd(eta0)
  
  k <- sigma_fac         # 4 或 7
  bound_fac <- 7        # 遠邊界倍率（可調大一點，確保多半不碰邊）
  B <- bound_fac * s_eta # 讓 η ∈ [-B, +B]
  
  q_prop <- 0.2          # ★ 只動 40% 的座標（可改 0.3~0.6）
  tau_dir <- 0.2        # ★ 方向微擾強度（0=不加噪；例如 0.1 會更不整齊）
  mix_prob  <- 0.2      # 新增：30% 的點不用「與 y 相反」方向
  hit <- 0L              # 邊界命中計數器
  
  for (i in contam_idx) {
    # 名目要推動的 Δη（在 logits 空間），方向跟 Y 相反
    sgn_base <- if (main_dat$Y[i] == 0) +1 else -1
    sgn <- if (runif(1) < mix_prob) sample(c(-1, +1), 1) else sgn_base
    d_move_nom <- sgn * (k * s_eta)
    
    # 邊界限制：最多只能走到 ±B
    eta_i <- eta0[i]
    d_max <- sgn * (B - sgn * eta_i)
    d_move <- sgn * min(abs(d_move_nom), abs(d_max))
    if (abs(d_move) >= abs(d_max) - 1e-12) hit <- hit + 1L
    
    # ★ 只動部分座標
    q <- max(1, floor(q_prop * p))
    j_sub <- sample.int(p, q, replace = FALSE)
    
    # 子向量方向（必要時加入微擾，避免方向過於整齊）
    b_sub <- beta_s[j_sub]
    b_norm_sub <- sqrt(sum(b_sub^2))
    if (b_norm_sub < 1e-12) next  # β 在這些座標太小就跳過（很罕見）
    
    u_sub <- b_sub / b_norm_sub
    if (tau_dir > 0) {
      z <- rnorm(length(j_sub))
      u_sub <- u_sub + tau_dir * z
      u_sub <- u_sub / sqrt(sum(u_sub^2))
    }
    
    # ★ 步長校正：用子向量長度，保證真正的 Δη = d_move
    t_step <- d_move / b_norm_sub
    Wmat[i, j_sub] <- Wmat[i, j_sub] + t_step * u_sub
  }
  
  cat(sprintf("boundary hits: %d/%d (%.1f%%)\n",
              hit, length(contam_idx), 100*hit/length(contam_idx)))
  
  # 寫回
  for (j in seq_len(p)) main_dat[[paste0("W", j)]] <- Wmat[, j]
  
  
  # 快速檢查污染是否推開
  eta_ref <- as.numeric(b0_true + Wmat %*% beta_s)
  mu_ref  <- plogis(eta_ref)
  cat(sprintf("[check] mu_ref bad: med=%.3f, good: med=%.3f\n",
              median(mu_ref[contam_idx]), median(mu_ref[-contam_idx])))
  
  # 好/壞組與類別檢查
  good_idx <- setdiff(seq_len(n_main), contam_idx)
  p1 <- mean(main_dat$Y[good_idx] == 1)
  if (p1 <= 0.1 || p1 >= 0.9) {
    warning("skip: good-sample class imbalance ", round(p1, 2))
    return(NULL)
  }
  G.ind <- good_idx
  if (length(unique(main_dat$Y)) < 2L) {
    warning("skip: overall Y is all the same after contamination")
    return(NULL)
  }
  
  ## 3-7 EM（RC 起始）
  message("   · EM (RC init) ...")
  Wmat <- as.matrix(main_dat[, paste0("W", 1:p)])  # 用污染後的 W
  ini  <- get_rc_init(main_dat$Y, Wmat,
                      sigma2_mat = delta2_mat[(n_val + 1):n_total, ],
                      verbose = TRUE)
  beta_init <- ini$beta
  init_type <- ini$type
  message("   · start = ", init_type)
  
  em_out <- hetero_em_logit(
    Y = main_dat$Y, W = Wmat, sigma2_mat = delta2_mat[(n_val + 1):n_total, ],
    init_beta = beta_init,
    B0 = 120, Bmax = 800, growth = 1.6, B_vcov = 1200,
    maxit = 300, tol = 1e-4, verbose = TRUE, vcov = "sandwich",
    seed = 123, use_sobol = TRUE, antithetic = TRUE
  )
  
  print(em_out$beta)
  
  # 影響點診斷
  infl <- diag_func(
    Y = main_dat$Y,
    W = Wmat,
    sigma2_mat = delta2_mat[(n_val + 1):n_total, ],
    cut_vec = cut_vec,
    G.ind = G.ind,
    em_out = em_out,
    dbg = FALSE
  )
  
  # ============ 多輪小步修剪（blind）開始 ============
  # 小工具：robust z-score（用 median/MAD，比 scale 穩）
  .z <- function(x) {
    m <- stats::median(x, na.rm = TRUE)
    s <- stats::mad(x, center = m, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
    (x - m) / s
  }
  
  keep <- seq_len(n_main)                 # 目前保留的索引（原始列編號）
  dropped <- integer(0)                   # 累計刪除
  trim_budget_frac <- 0.10                # 最多刪 10%（實務上保守）
  trim_budget <- floor(trim_budget_frac * n_main)
  max_round <- 3                          # 1～3 輪通常足夠
  em_curr <- em_out                       # 先用第一次 EM 當起點
  
  for (iter in seq_len(max_round)) {
    message(sprintf("[trim] round %d | keep=%d/%d", iter, length(keep), n_main))
    
    # 針對「目前 keep 集」重算診斷（門檻沿用原 cut_vec）
    infl_sub <- diag_func(
      Y = main_dat$Y[keep],
      W = as.matrix(main_dat[keep, paste0("W", 1:p), drop = FALSE]),
      sigma2_mat = delta2_mat[(n_val + 1):n_total, ][keep, , drop = FALSE],
      cut_vec = cut_vec,
      G.ind = match(intersect(G.ind, keep), keep),    # 轉成子樣本位置
      em_out = em_curr,
      dbg = FALSE
    )
    
    # infl_sub$id 是 1..length(keep)，把它對應回原始列編號
    idx_keep <- keep[infl_sub$id]
    
    # ---- 共識旗標（≥ 2 項命中才算可疑；門檻照舊）----
    flag_rcs <- infl_sub$R_CS       > cut_vec["StdR2_CS"]
    flag_mdf <- abs(infl_sub$MDF)   > cut_vec["GDFFITS"]
    flag_mcd <- infl_sub$mCDstar    > cut_vec["mCDstar"]
    hit_cnt  <- flag_rcs + flag_mdf + flag_mcd
    
    cand_local <- which(hit_cnt >= 2L)                     # 子樣本位置
    if (!length(cand_local)) cand_local <- which(hit_cnt > 0)  # 保底
    if (!length(cand_local)) { message("[trim] no candidates; stop."); break }
    
    # 綜合分數排名（越大越像壞點）
    score_local <- .z(infl_sub$R_CS) + .z(abs(infl_sub$MDF)) + .z(infl_sub$mCDstar)
    ord_local   <- order(score_local[cand_local], decreasing = TRUE)
    
    # ---- 決定這一輪刪幾個：只刪「可疑的一半」、但總量 ≤ 10% n ----
    remain <- trim_budget - length(dropped)
    if (remain <= 0L) { message("[trim] budget used up; stop."); break }
    drop_n <- min(remain, max(0L, floor(0.5 * length(cand_local))))
    if (drop_n <= 0L) { message("[trim] drop_n=0; stop."); break }
    
    # 轉回原始列編號後，選出要刪的 index
    cand_global <- idx_keep[cand_local]
    to_drop     <- cand_global[ord_local][seq_len(drop_n)]
    
    # 更新 keep / dropped
    keep    <- setdiff(keep, to_drop)
    dropped <- c(dropped, to_drop)
    
    # 類別防呆：兩類都要有、且不極端
    yk  <- main_dat$Y[keep]
    p1k <- mean(yk == 1)
    if (length(unique(yk)) < 2L || p1k <= 0.05 || p1k >= 0.95) {
      message("[trim] class collapsed after trimming; revert and stop.")
      keep    <- setdiff(keep, integer(0))       # 保留現狀
      dropped <- dropped[seq_len(length(dropped) - drop_n)]  # 回退剛剛的刪除
      break
    }
    
    # 用 keep 集重估 ΣX、重新 RC 起始，再配 EM（小步）
    SigmaX_keep <- .est_SigmaX_full(
      W = as.matrix(main_dat[keep, paste0("W", 1:p), drop = FALSE]),
      sigma2_mat = delta2_mat[(n_val + 1):n_total, ][keep, , drop = FALSE]
    )
    ini2 <- get_rc_init(
      Y = main_dat$Y[keep],
      W = as.matrix(main_dat[keep, paste0("W", 1:p), drop = FALSE]),
      sigma2_mat = delta2_mat[(n_val + 1):n_total, ][keep, , drop = FALSE],
      verbose = FALSE
    )
    
    em_curr <- hetero_em_logit(
      Y = main_dat$Y[keep],
      W = as.matrix(main_dat[keep, paste0("W", 1:p)]),
      sigma2_mat = delta2_mat[(n_val + 1):n_total, ][keep, , drop = FALSE],
      init_beta = ini2$beta,
      B0 = 150, Bmax = 600, growth = 1.5, B_vcov = 1000,
      maxit = 250, tol = 1e-4, vcov = "sandwich", verbose = FALSE,
      SigmaX = SigmaX_keep, use_sobol = TRUE, antithetic = TRUE
    )
    
    # 若這輪沒有新的可疑或改善幅度很小，可以提前停（可選）
    if (iter == max_round) break
  }
  # 最終 EM 結果
  em_out <- em_curr
  # ============ 多輪小步修剪（blind）結束 ============
  
  
  # 只取斜率（依名稱對齊）
  nm <- colnames(Wmat)  # 原始欄名順序
  b <- tryCatch(em_out$beta, error = function(e) NULL)
  if (is.null(b) || any(!is.finite(b)) || max(abs(b)) > 15) {  # 15 比 25 更嚴
    warning("[perf] 本回合 β 極端/非有限：不計入 perf（設 NA）")
    beta_EM <- rep(NA_real_, length(nm))
  } else {
    slope_names <- setdiff(names(b), "(Intercept)")
    beta_EM <- as.numeric(b[slope_names][match(nm, slope_names)])
  }
  
  # 效能表
  get_bias_mse <- function(est, beta_true) {
    d <- est - beta_true
    safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    c(bias = safe_mean(d), mse = safe_mean(d^2))
  }
  
  perf_df <- data.frame(method = "EM_trim",
                        t(get_bias_mse(beta_EM, beta_true[-1])))
  perf_df$kappa <- kappa_tar
  perf_df$init_type <- init_type
  perf_df$estimator <- "TrimEM"
  
  
  # 指標匯總
  bad_RCS  <- infl$R_CS[infl$id %in% contam_idx]
  good_RCS <- infl$R_CS[infl$id %in% setdiff(seq_len(n_main), contam_idx)]
  cs_cut <- cut_vec["StdR2_CS"]
  qs_bad  <- quantile(bad_RCS,  probs = c(.5,.9,.95,.99), na.rm = TRUE)
  qs_good <- quantile(good_RCS, probs = c(.5,.9,.95,.99), na.rm = TRUE)
  message(sprintf("[R_CS] cut=%.3g | bad med=%.3g p95=%.3g p99=%.3g | good med=%.3g p95=%.3g",
                  cs_cut, qs_bad[1], qs_bad["95%"], qs_bad["99%"], qs_good[1], qs_good["95%"]))
  
  if (!identical(infl$id, seq_len(n_main))) infl <- infl[order(infl$id), , drop = FALSE]
  flag_mat <- as.matrix(infl[, grepl("^flag_", names(infl)), drop = FALSE]); flag_mat[is.na(flag_mat)] <- FALSE
  is_bad <- seq_len(n_main) %in% contam_idx
  TP <- colSums(flag_mat[ is_bad, , drop = FALSE])
  FP <- colSums(flag_mat[!is_bad, , drop = FALSE])
  CIR <- TP / sum(is_bad)
  SR  <- FP / (n_main - sum(is_bad))
  
  labels <- as.integer(is_bad)
  scores <- list(
    CookD = infl$CookD, DFFITS = abs(infl$DFFITS),
    GDF = abs(infl$GDF), MDF = abs(infl$MDF),
    GDB = abs(infl$GDB), MDB = abs(infl$MDB),
    GD = infl$GD, MD = infl$MD,
    GCD_GSPR = abs(infl$GCD_GSPR), mCDstar = infl$mCDstar,
    R_CS = infl$R_CS_raw, R_MF = infl$R_MF_raw, R_NK = infl$R_NK_raw
  )
  AUC_fast <- function(y, p) {
    y <- as.integer(y); pos <- p[y == 1]; neg <- p[y == 0]
    n1 <- length(pos); n0 <- length(neg); if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- rank(c(pos, neg)); (sum(r[1:n1]) - n1*(n1+1)/2) / (n1*n0)
  }
  auc_by_method <- vapply(scores, function(s) {
    ok <- is.finite(s) & !is.na(labels)
    if (sum(ok) < 2 || length(unique(labels[ok])) < 2 || length(unique(s[ok])) < 2) return(NA_real_)
    AUC_fast(labels[ok], s[ok])
  }, numeric(1))
  
  diag_df <- data.frame(method = sub("^flag_", "", colnames(flag_mat)), CIR = CIR, SR = SR, kappa = kappa_tar)
  diag_df$init_type <- init_type
  diag_df$estimator <- "TrimEM"
  diag_df$AUC <- as.numeric(auc_by_method[match(diag_df$method, names(auc_by_method))])
  
  out <- list(perf = perf_df, diag = diag_df)
  attr(out, "infl") <- infl
  return(out)
}

# ---------- 4. 主迴圈 ----------
set.seed(2025)  # Monte‑Carlo 共用種子
result_perf <- list(); result_diag <- list(); id <- 1L

for (diag_mode in diag_mode_list) {
  # 這裡把 ... 直接 forward 到 diag_all_1stepEM，
  # 以便 one_run() 可以傳 dbg 等參數進去
  diag_func <- switch(
    diag_mode,
    em_step = function(Y, W, sigma2_mat, cut_vec, G.ind = NULL, em_out = NULL, ...) {
      if (is.null(em_out)) stop("em_out 應由 one_run() 提供")
      diag_all_1stepEM(
        Y = Y, W = W, sigma2_mat = sigma2_mat,
        G.ind = G.ind, cut_vec = cut_vec, em_out = em_out, ...
      )
    }
  )
  
  for (p in p_list) {
    for (n_total in n_list) {
      for (contam_rate in contam_rate_list) {
        for (sigma_fac in sigma_fac_list) {
          message(sprintf("▶ run %d | mode=%s p=%d n=%d contam=%.2g σ=%.2f",
                          id, diag_mode, p, n_total, contam_rate, sigma_fac))
          contam_n <- floor(n_total * contam_rate)
          n_val <- 0L
          n_main   <- n_total - n_val
          beta_true <- c(1, rep(2, p))
          cut_vec   <- make_cut_vec(p, n_main, contam_n)
          
          setting_res <- future_lapply(seq_len(R), function(r) {
            lapply(kappa_set, one_run,
                   n_total   = n_total,
                   n_val     = n_val,
                   p         = p,
                   beta_true = beta_true,
                   contam_n  = contam_n,
                   sigma_fac = sigma_fac,
                   cut_vec   = cut_vec,
                   diag_func = diag_func)
          }, future.seed = TRUE, future.scheduling = 1)
          
          # 把所有 one_run() 結果攤平成同一層，並過濾掉 NULL
          flat_res <- Filter(Negate(is.null), unlist(setting_res, recursive = FALSE))
          
          # perf
          perf_part <- rbindlist(lapply(flat_res, `[[`, "perf"),  fill = TRUE)
          # diag
          diag_part <- rbindlist(lapply(flat_res, `[[`, "diag"),  fill = TRUE)
          
          perf_part[, `:=`(diag_mode = diag_mode, p = p, n_total = n_total, contam_rate = contam_rate, sigma_fac = sigma_fac)]
          diag_part[, `:=`(diag_mode = diag_mode, p = p, n_total = n_total, contam_rate = contam_rate, sigma_fac = sigma_fac)]
          
          result_perf[[id]] <- perf_part
          result_diag[[id]] <- diag_part
          id <- id + 1
        }
      }
    }
  }
}

## 合併所有模式的結果
final_perf <- rbindlist(result_perf)
final_diag <- rbindlist(result_diag)

## ================= 5. 摘要輸出 ==================
summary_tbl <- final_perf %>%
  group_by(method, kappa, p, n_total, contam_rate, sigma_fac) %>%
  summarise(
    mean_bias = if (any(is.finite(bias))) mean(bias[is.finite(bias)]) else NA_real_,
    sd_bias   = if (sum(is.finite(bias)) > 1) sd(bias[is.finite(bias)]) else NA_real_,
    mean_mse  = if (any(is.finite(mse)))  mean(mse[is.finite(mse)])   else NA_real_,
    .groups   = "drop"
  ) %>% arrange(method)
print(summary_tbl, n = Inf)

## ── 影響點診斷成功率 ───────────────────────────
summary_all <- final_diag %>%
  group_by(diag_mode, method, kappa, p, n_total, contam_rate, sigma_fac) %>%
  summarise(
    CIR_mean = mean(CIR, na.rm = TRUE),
    CIR_sd   = sd(CIR,  na.rm = TRUE),
    SR_mean  = mean(SR,  na.rm = TRUE),
    SR_sd    = sd(SR,   na.rm = TRUE),
    BA_mean  = mean((CIR + (1 - SR))/2, na.rm = TRUE),
    BA_sd    = sd((CIR + (1 - SR))/2,   na.rm = TRUE),
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd   = sd(AUC,  na.rm = TRUE),
    .groups = "drop"
  )

# 原本的 diag_rate（只要平均）
diag_rate <- summary_all %>%
  dplyr::select(diag_mode, method, kappa, p, n_total, contam_rate, sigma_fac,
                CIR_mean, SR_mean, AUC_mean)
print(diag_rate, n = Inf, width = Inf)

# 原本的 compare_multi（有 BA 與 sd）
methods_keep <- c("R_CS","R_MF","R_NK","GDF","GCD_GSPR","MDB","MDF","GDB","mCDstar","CookD","DFFITS")
compare_multi <- summary_all %>% dplyr::filter(method %in% methods_keep)
print(compare_multi, n = Inf, width = Inf)

write_xlsx(
  x    = compare_multi,
  path = file.path("C:/Users/kimi1/OneDrive/文件/論文/run200 .xlsx")
)