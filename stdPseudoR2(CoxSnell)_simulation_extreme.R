# Std Pseudo R-square (CoxSnell): cut-off (3)
# 高度影響點 (+7σ)
rm(list = ls()) # 清空 Environment
cat("\014") # 清空 Console

library(DescTools) # PseudoR2()
library(openxlsx) # 匯出.xlsx檔(綜合分析結果)

# setting for simulation -------------------------------------------------------------------------------------
n_rep <- 100                                     # 重複數
n_list <- c(20, 30, 50, 100, 1000)    # 樣本數
# n_list <- c(20, 30, 50, 100, 150, 250, 500, 1000)    # 樣本數
IDvariable_list <- c(1, 2, 5)                   # 自變數(X)個數
ContaminR_list <- c((0.1), (0.2), (0.3))        # 汙染率
c <- 3

dataAnalysisPath <- "D:/OL"

folder_nameSave <- vector()
result_column_names <- c(
  "Cook's_Distance", "IDENTIFY_CD", "DFFITS", "IDENTIFY_DFFITS",
  "GDFFITS", "IDENTIFY_GDFFITS", "GSDFBETA", "IDENTIFY_GSDFBETA",
  "GCD_GSPR", "IDENTIFY_GCD_GSPR", "mCD", "IDENTIFY_mCD",
  "Std_Pseudo_Rsqu (CoxSnell)", "IDENTIFY_Std_Pseudo_Rsqu (CoxSnell)"
)

# CIR、SR 欄位名稱
CS_column_names <- c(
  "Cook's_Distance", "DFFITS", "GDFFITS",
  "GSDFBETA", "GCD_GSPR", "mCD*", "Std_Pseudo_Rsqu (CoxSnell)"
)
# 總報表格式設定
Result_column_names <- c(
  "# Variable", "Sample Size", "Contamination Rates", "Evaluation",
  "Cook's_Distance", "DFFITS", "GDFFITS", "GSDFBETA",
  "GCD_GSPR", "mCD*", "Std_Pseudo_Rsqu (CoxSnell)"
)
extrResult <- matrix(NA,
                     nrow = length(IDvariable_list)*length(n_list)*length(ContaminR_list)*2,
                     ncol = length(Result_column_names))
colnames(extrResult) <- Result_column_names
rec <- 1 # 計數此次設定下，不同方法的 CIR 與 SR

# Helper Functions -------------------------------------------------------------------------------------------
detect_outliers <- function(values, threshold) {
  abs_values <- abs(values)
  detected <- abs_values > threshold
  if (sum(na.omit(detected)) == 0) { # ignore values with "NaN" to sum up
    # some values is "NaN", so identification of its denote as "NA"
    return(ifelse(is.na(values), NA, 0))
  } else {
    result <- rep(0, length(values))
    result[detected] <- 1
    # some values is "NaN", so identification of its denote as "NA"
    return(ifelse(is.na(values), NA, result))
  }
}

folder_times <- 0
# Generate Data & Compute ------------------------------------------------------------------------------------
for (IDvariable in IDvariable_list) { # 自變數數量
  # beta <- matrix(c(1, rep(2, IDvariable)), (IDvariable + 1), 1)   # (IDvariable + 1)*1 matrix
  if (IDvariable == 1) {
    beta <- matrix(c(1, 2), (IDvariable + 1), 1)
  } else if (IDvariable == 2) {
    beta <- matrix(c(1, 2, 2), (IDvariable + 1), 1)
  } else if (IDvariable == 5) {
    beta <- matrix(c(1, 2, 2, 2, 2, 2), (IDvariable + 1), 1)
  }
  
  # generate parameter (std normal distribution) of each independent variable
  mean_ <- rep(0, IDvariable)
  sd_ <- rep(1, IDvariable)
  
  data_names <- character(0)
  data_names <- c(
    data_names,
    paste("X", 1:IDvariable, sep = ""),
    "y",
    "Contamin",
    paste("extrX", 1:IDvariable, sep = ""),
    "NEWy"
  )
  
  beta_column_names <- character(0)
  beta_column_names <- c(
    beta_column_names,
    paste("Estimate_beta", 0:IDvariable, sep = ""),
    paste("Std.Error_beta", 0:IDvariable, sep = ""),
    paste("Z-value_beta", 0:IDvariable, sep = ""),
    paste("Pr(>|z|)_beta", 0:IDvariable, sep = "")
  )
  
  for (n in n_list) { # 樣本數數量
    for (ContaminR in ContaminR_list) { # 汙染點比例
      ContaminS <- n * ContaminR
      
      # Cut-off values setting -------------------------------------------------------------------------
      d <- ContaminS                                              # number of points that are suspicious
      cutoff_CD <- 1
      cutoff_DFFITS <- 3 * sqrt((IDvariable + 1) / n)
      cutoff_GDFFITS <- 3 * sqrt((IDvariable + 1) / (n - d))
      cutOff_GDFBETA <- (9 * (IDvariable + 1)) / (n - d - 3 * IDvariable)
      cutoff_GCD_GSPR <- 1
      cutoff_mCD <- c / sqrt((n - IDvariable)/(n - d))
      cutoff_stdpseudoR2 <- 3
      
      # Put all cut-off values into a vector
      cutoff <- c(
        cutoff_CD, cutoff_DFFITS,
        cutoff_GDFFITS, cutOff_GDFBETA, 
        cutoff_GCD_GSPR, cutoff_mCD, 
        cutoff_stdpseudoR2
      )
      
      # evaluation setting -----------------------------------------------------------------------------
      extrCIR <- matrix(nrow = n_rep, ncol = length(cutoff))        # correct identification rate (extreme)
      extrSR <- matrix(nrow = n_rep, ncol = length(cutoff))         # swamping rate (extreme)
      
      # record Estimate, Std. Error, z value, Pr(>|z|) of beta in MODEL_All & MODEL_R
      extrbetaA <- matrix(nrow = n_rep, ncol = 4 * (IDvariable + 1))
      extrbetaR <- matrix(nrow = n_rep, ncol = 4 * (IDvariable + 1))
      
      cat("Extreme | # Variable:", IDvariable, "| Sample Size:", n,  "| Contamination Rate:", ContaminR,
          "| Contaminated Sample:", ContaminS, "\n")
      
      # 建立存放「模擬資料」與「分析結果」的資料夾
      folder_times <- folder_times + 1
      folder_name <- paste0(dataAnalysisPath, "/Extreme/V", IDvariable, " n=", n,
                            " CR=", ContaminR)
      folder_nameSave <- rbind(folder_nameSave, folder_name) # 紀錄存放的資料夾名稱
      if (!file.exists(folder_name)) { # 如果資料夾不存在則需建立
        dir.create(folder_name, recursive = TRUE) # 使用recursive = TRUE來建立多層文件夾
        cat(paste0("\033[36m資料夾建立成功！路徑：", folder_name, "\033[0m\n"))
      } else { # 如果資料夾存在，跳出覆蓋檔案的訊息
        cat(paste0("\033[36m", folder_name, "，此資料夾已存在。注意！舊檔案將被覆蓋。\033[0m\n"))
      }
      r <- 1
      while (r <= n_rep) {
        cat(r, "\n")
        # Results data frame setting -------------------------------------------------------------------
        extrDatRES <- as.data.frame(matrix(ncol = length(cutoff) * 2, nrow = n)) # 7 methods * 2
        colnames(extrDatRES) <- result_column_names
        
        Contamin <- sample(n, ContaminS) # draw the sample to be contaminated
        
        # independent variable (X)
        x <- matrix(NA, nrow = n, ncol = IDvariable)
        for (i in 1:IDvariable) {
          x[, i] <- rnorm(n, mean = mean_[i], sd = sd_[i]) # Generate one independent variable
        }
        
        # add a column with value 1 as x0
        X <- cbind(1, x)                                          # n*(IDvariable + 1) matrix
        
        linpred <- X %*% beta
        prob <- 1/(1 + exp(-linpred))
        y <- matrix(NA, nrow = n, ncol = 1)
        y = rbinom(n = n, size = 1, prob = prob)
        
        # add a column (contaminated) to record which point is contaminated
        contaminated <- matrix(ncol = 1, nrow = n)
        contaminated[Contamin] = 1
        contaminated[-Contamin] = 0
        # number of variable*2(original、extreme)+3(y、Contamin、NEWy)
        dat <- matrix(nrow = n, ncol = IDvariable * 2 + 3)
        dat[, 1:IDvariable] <- x
        dat[, IDvariable + 1] <- y
        dat[, IDvariable + 2] <- contaminated
        # 製造汙染點
        for (idv in 1:IDvariable){
          for (j in 1:n) {
            if (dat[j, (IDvariable + 2)] == 1){ # 汙染點的樣本處理
              dat[j, idv + (IDvariable + 2)] = dat[j, idv] + 7 * sd_[idv]
              dat[j, (2 * IDvariable + 3)] = 0 # contaminated point its y is set to 0
            } else{ # 非汙染點的樣本
              dat[j, idv + (IDvariable + 2)] = dat[j, idv]
              dat[j, (2 * IDvariable + 3)] = dat[j, IDvariable + 1] # original y
            }
          }
        }
        
        if (sum(dat[, (2 * IDvariable + 3)][dat[, IDvariable + 2] == 0]) / (n-ContaminS) <= 0.1 ||
            sum(dat[, (2 * IDvariable + 3)][dat[, IDvariable + 2] == 0]) / (n-ContaminS) >= 0.9) {
          r <- r+1 # 不平衡資料就跳過(newY)
        } else {
          colnames(dat) <- data_names
          # 儲存生成的資料
          yourfilenamecsv = paste0(folder_name,
                                   "/V", IDvariable, "_n", n, "_CR", ContaminR,
                                   "(", r, ").csv"
          )
          write.csv(dat, file = yourfilenamecsv)
          
          # Extract independent variable which is contaminated with EXTREME suspicious observations
          extrDatO <- dat[, c((IDvariable + 2):(IDvariable * 2 + 3))]
          extrDatO <- as.data.frame(extrDatO)
          extrDatO$NEWy <- as.factor(extrDatO$NEWy)
          # Logistic regression model using the original data
          MODEL_All <- glm(NEWy ~ . - Contamin - NEWy, data = extrDatO,
                           family = "binomial", maxit = 100)
          # summary(MODEL_All)
          
          # Keep remaining data set
          extrDatR <- extrDatO[-Contamin, ]
          # Logistic regression model using the remaining data
          MODEL_R <- glm(NEWy ~ . - Contamin - NEWy, data = extrDatR, family = "binomial",
                         maxit = 100)
          # summary(MODEL_R)
          
          # record information about beta
          coefA <- numeric()
          coefR <- numeric()
          for (e in 1:4) { # Estimate, Std. Error, z value, Pr(>|z|)
            for(v in 1:(IDvariable + 1)){
              coefA <- c(coefA, coef(summary(MODEL_All))[v, c(e)])
              coefR <- c(coefR, coef(summary(MODEL_R))[v, c(e)])
            }
          }
          
          extrbetaA[r, ] <- coefA
          extrbetaR[r, ] <- coefR
          
          ## Cook's Distance #############################################################################
          extrDatRES[, 1] <- cooks.distance(MODEL_All)
          extrDatRES[, 2] <- detect_outliers(extrDatRES[, 1], cutoff[1])
          
          ## DFFITS ######################################################################################
          extrDatRES[, 3] <- dffits(MODEL_All)
          extrDatRES[, 4] <- detect_outliers(extrDatRES[, 3], cutoff[2])
          
          # Data Preparation -----------------------------------------------------------------------------
          X <- matrix(NA, n, (IDvariable + 1))                  # n*(IDvariable+1) matrix
          pi_i_R <- numeric(n)                                  # n vector
          vi_R <- numeric(n)                                    # n vector
          X_R <- matrix(NA, nrow(extrDatR), (IDvariable + 1))   # (n-Contamin)*(IDvariable+1) matrix
          PI_i_R <- numeric(nrow(extrDatR))                     # (n-Contamin) vector
          Vi_R <- numeric(nrow(extrDatR))                       # (n-Contamin) vector
          V_R <- matrix(NA, nrow(extrDatR), nrow(extrDatR))     # (n-Contamin)*(n-Contamin) diagonal matrix
          Mid <- matrix(NA, (IDvariable + 1), (IDvariable + 1)) # (IDvariable + 1)*(IDvariable + 1) matrix
          term <- numeric(n)                                    # n vector
          hii_R <- numeric(n)                                   # n vector
          y_i <- numeric(n)                                     # n vector
          
          # *** vi_R ***
          pi_i_R <- predict(MODEL_R, extrDatO, type = 'response') # 移除可疑點後，計算完整資料的pi_hat
          vi_R <- pi_i_R * (1 - pi_i_R)
          # *** x_i ***
          X <- model.matrix(MODEL_All) # 完整資料自變數(包含截距項)
          # *** X_R ***
          X_R <- model.matrix(MODEL_R) # 未汙染資料自變數(包含截距項)
          # *** V_R ***
          PI_i_R <- predict(MODEL_R, extrDatR, type = 'response') # 移除可疑點後，計算未汙染資料的pi_hat
          Vi_R <- PI_i_R * (1 - PI_i_R)
          V_R <- diag(Vi_R)
          # *** hii_R ***
          Mid <- solve(t(X_R) %*% V_R %*% X_R)
          term <- X %*% Mid %*% t(X)
          hii_R <- diag(vi_R * term) # 取出 leverage matrix 對角線的元素，存成 hii_R（完整資料個數）
          
          # compute the generalized standardized Pearson residual (GSPR) and
          #         the generalized weights (GWs)
          r_si_R <- numeric(n)                                  # n vector (GSPR)
          hiiS_R <- numeric(n)                                  # n vector (GWs)
          y_i <- as.numeric(extrDatO$NEWy)-1
          for (k in 1:n) {
            if (k %in% Contamin){
              r_si_R[k] <- (y_i[k] - pi_i_R[k]) / sqrt(vi_R[k] * (1 + hii_R[k]))
              hiiS_R[k] <- hii_R[k] / (1 + hii_R[k])
            } else{
              r_si_R[k] <- (y_i[k] - pi_i_R[k]) / sqrt(vi_R[k] * (1 - hii_R[k]))
              hiiS_R[k] <- hii_R[k] / (1 - hii_R[k])
            }
          }
          
          ## GDFFITS #####################################################################################
          GDFFITS <- numeric(n)
          GDFFITS <- r_si_R * sqrt(hiiS_R)
          
          extrDatRES[, 5] <- GDFFITS
          extrDatRES[, 6] <- detect_outliers(extrDatRES[, 5], cutoff[3])
          
          ## GDFBETAS ####################################################################################
          GSDFBETAS <- numeric(n)
          for (m in 1:n) {
            if (m %in% Contamin){
              GSDFBETAS[m] <- (hiiS_R[m] * (r_si_R[m] ^ 2)) / (1 + hii_R[m])
            } else{
              GSDFBETAS[m] <- (hiiS_R[m] * (r_si_R[m] ^ 2)) / (1 - hii_R[m])
            }
          }
          
          extrDatRES[, 7] <- GSDFBETAS
          extrDatRES[, 8] <- detect_outliers(extrDatRES[, 7], cutoff[4])
          
          ## GCD.GSPR ####################################################################################
          GCD_GSPR <- numeric(n)
          GCD_GSPR <- (1/(IDvariable + 1)) * (r_si_R ^ 2) * diag(term)
          
          extrDatRES[, 9] <- GCD_GSPR
          extrDatRES[, 10] <- detect_outliers(extrDatRES[, 9], cutoff[5])
          
          ## mCD #########################################################################################
          mCD <- numeric(n)
          mCD <- abs(GDFFITS) * sqrt((n - d - (IDvariable + 1)) / (IDvariable + 1))
          
          extrDatRES[, 11] <- mCD
          extrDatRES[, 12] <- detect_outliers(extrDatRES[, 11], cutoff[6])
          
          ## Standard pseudo R^2: "CoxSnell" #############################################################
          # https://search.r-project.org/CRAN/refmans/DescTools/html/PseudoR2.html
          # PseudoR2(MODEL_All, which = "all")
          # PseudoR2(MODEL_R, which = "all")
          # pseudo_name <- c("McFadden", "McFaddenAdj", "CoxSnell",
          #                  "Nagelkerke", "AldrichNelson", "VeallZimmermann",
          #                  "Efron", "McKelveyZavoina", "Tjur")
          
          pseudoR2 <- vector("numeric")
          stdpseudoR2 <- vector("numeric")
          
          # 分組計算 pseudoR2
          for (p in 1:n) {
            if (p %in% Contamin) { # 如果是壞點組
              datSTP <- rbind(extrDatR, extrDatO[p, ]) # 好點組加入一個壞的點
              MODEL <- glm(NEWy ~ . - Contamin - NEWy, data = datSTP,
                           family = "binomial", maxit = 100)
              pseudoR2[p] <- PseudoR2(MODEL, which = "CoxSnell")
            } else { # 好點組
              datSTP <- extrDatO[-c(p, Contamin), ] # 移除好點組的其中一點
              MODEL <- glm(NEWy ~ . - Contamin - NEWy, data = datSTP,
                           family = "binomial", maxit = 100)
              pseudoR2[p] <- PseudoR2(MODEL, which = "CoxSnell")
            }
          }
          pseudoR2 # 分組計算結果
          # plot(pseudoR2)
          
          # 「pseudoR2」與「所有好點」的 PseudoR2 差異
          pseudoR2_diff <- PseudoR2(MODEL_R, which = "CoxSnell") - pseudoR2
          # plot(pseudoR2_diff)
          
          sdB <- sd(pseudoR2_diff[Contamin])   # 壞點組差異的標準差
          sdG <- sd(pseudoR2_diff[-Contamin])  # 好點組差異的標準差
          
          # 分組計算 stdpseudoR2
          for (b in 1:n) {
            if (b %in% Contamin) { # 如果是壞點組
              stdpseudoR2[b] <- abs(pseudoR2_diff[b]) / sdB
            } else { # 好點組
              stdpseudoR2[b] <- abs(pseudoR2_diff[b]) / sdG
            }
          }
          stdpseudoR2
          # plot(stdpseudoR2)
          
          extrDatRES[, 13] <- stdpseudoR2
          extrDatRES[, 14] <- detect_outliers(extrDatRES[, 13], cutoff[7])
          
          extrDatRES <- cbind(extrDatRES, Contamin = extrDatO$Contamin)
          # write.csv(extrDatRES, "C:/Users/jenny/Downloads/extrDatRES.csv")
          
          # 先計算 CIR 和 SR
          for (cs in 1:length(cutoff)) {
            NArow <- sum(is.na(extrDatRES[, (cs * 2 - 1)]))
            # CIR：Contamin = 1 (影響點) & IDENTIFY = 1（正確找到影響點）
            extrCIR[r, cs] <- nrow(extrDatRES[extrDatRES[, (length(cutoff) * 2 + 1)] == 1 &
                                                extrDatRES[, (cs * 2)] == 1, ]) / ContaminS - NArow
            # SR：Contamin = 0 (影響點) & IDENTIFY = 1（不是影響點卻被誤判成影響點）
            extrSR[r, cs] <- nrow(extrDatRES[extrDatRES[, (length(cutoff) * 2 + 1)] == 0 &
                                               extrDatRES[, (cs * 2)] == 1, ]) / (n - ContaminS - NArow)
          }
          cutoff <- c(
            cutoff_CD, cutoff_DFFITS,
            cutoff_GDFFITS, cutOff_GDFBETA,
            cutoff_GCD_GSPR, cutoff_mCD,
            cutoff_stdpseudoR2
          )
          cutoffArray <- rep(NA, length(cutoff)*2+1)
          for (cut in 1:length(cutoff)) {
            cutoffArray[(cut*2)-1] <- cutoff[cut]
          }
          extrDatRES <- rbind(cutoffArray, extrDatRES) # 將臨界值與分析結果合併
          rownames(extrDatRES) <- c("cut-off", 1:n)
          
          # 匯出每筆資料不同方法的異常點檢測結果
          yourfilenamecsv2 = paste0(folder_name,
                                    "/V", IDvariable, "_extrDatResult(", r, ").csv")
          write.csv(extrDatRES, file = yourfilenamecsv2)
          r <- r + 1
        }
      }
      
      # 計算 CIR 平均 (extr) ---------------------------------------------------------------------------
      CIR_column_means <- colMeans(extrCIR)
      extrCIR <- rbind(extrCIR, CIR_column_means) # 合併平均
      colnames(extrCIR) <- CS_column_names
      rownames(extrCIR) <- c(1:n_rep, "mean")
      extrCIR <- cbind("Rep" = c(1:n_rep, "Mean"), extrCIR) # 為了匯出EXCEL
      
      # 計算 SR 平均 (extr) ----------------------------------------------------------------------------
      SR_column_means <- colMeans(extrSR)
      extrSR <- rbind(extrSR, SR_column_means)
      colnames(extrSR) <- CS_column_names
      rownames(extrSR) <- c(1:n_rep, "mean")
      extrSR <- cbind("Rep" = c(1:n_rep, "Mean"), extrSR) # 為了匯出EXCEL
      
      # BETA of MODEL_All (extr) -----------------------------------------------------------------------
      betaA_column_means <- colMeans(extrbetaA)
      extrbetaA <- rbind(extrbetaA, betaA_column_means)
      colnames(extrbetaA) <- beta_column_names
      rownames(extrbetaA) <- c(1:n_rep, "mean")
      extrbetaA <- cbind("Rep" = c(1:n_rep, "Mean"), extrbetaA) # 為了匯出EXCEL
      
      # BETA of MODEL_R (extr) -------------------------------------------------------------------------
      betaR_column_means <- colMeans(extrbetaR)
      extrbetaR <- rbind(extrbetaR, betaR_column_means)
      colnames(extrbetaR) <- beta_column_names
      rownames(extrbetaR) <- c(1:n_rep, "mean")
      extrbetaR <- cbind("Rep" = c(1:n_rep, "Mean"), extrbetaR) # 為了匯出EXCEL
      
      # 匯出某設定下重複完成後計算分析結果的 EXCEL 檔案 ------------------------------------------------
      wb <- createWorkbook() # 建立活頁簿
      # 將分析結果依序放入EXCEL中不同的Sheet
      addWorksheet(wb, "1 CIR")
      writeData(wb, "1 CIR", extrCIR)
      addWorksheet(wb, "2 SR")
      writeData(wb, "2 SR", extrSR)
      addWorksheet(wb, "3 beta MODEL_All")
      writeData(wb, "3 beta MODEL_All", extrbetaA)
      addWorksheet(wb, "4 beta MODEL_R")
      writeData(wb, "4 beta MODEL_R", extrbetaR)
      # 儲存
      current_time <- format(Sys.time(), "%Y-%m-%d %H-%M-%S") # 獲取當下的時間
      file_name <- paste0(dataAnalysisPath,
                          "/Extreme/V", IDvariable, "_n", n, "_CR", ContaminR,
                          "_Output (", current_time, ").xlsx")
      saveWorkbook(wb, file_name, overwrite = TRUE)
      
      # 記錄此次設定下，不同方法 CIR 與 SR 的平均 ------------------------------------------------------
      extrResult[rec, ] <- c(IDvariable, n, ContaminR, "CIR", round(CIR_column_means, 4))
      rec <- rec + 1
      extrResult[rec, ] <- c(IDvariable, n, ContaminR, "SR", round(SR_column_means, 4))
      rec <- rec + 1
    }
  }
}

# 綜合結果報表 -----------------------------------------------------------------------------------------
current_time <- format(Sys.time(), "%Y-%m-%d %H-%M-%S") # 獲取當下的時間
write.csv(extrResult,
          file = paste0(dataAnalysisPath,
                        "/Extreme/# extremeResult(", current_time, ").csv"),
          row.names = F)
