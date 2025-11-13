library("gmnl")
library("mlogit")
library(R.matlab)
library(MASS)
library(readxl)
set.seed(2025)

eff_code <- function(X, nlevels) {
  nf <- length(nlevels)  
  df <- sum(nlevels) - nf  
  xrow <- numeric(df)  
  
  startidx <- 1  
  
  for (i in seq_len(nf)) {
    nl <- nlevels[i]  
    xtmp <- numeric(nl - 1) 
    
    if (X[i] < nl) {
      xtmp[X[i]] <- 1 
    } else {
      xtmp[] <- -1  
    }
    
    xrow[startidx:(startidx + nl - 2)] <- xtmp
    startidx <- startidx + nl - 1  
  }
  
  return(xrow)
}

transform_X <- function(X, nlevels) {
  numRows <- nrow(X) 
  totalLevels <- sum(nlevels) - length(nlevels) 
  
  X_code <- matrix(0, nrow = numRows, ncol = totalLevels)
  
  for (row in 1:numRows) {
    X_code[row, ] <- eff_code(X[row, ], nlevels)
  }
  
  return(X_code)
}

transform_X_inter <- function(X, nlevels, interactions) {
  
  colIndex <- c(1, cumsum(nlevels - 1) + 1)
  colIndex <- colIndex[-length(colIndex)]
  
  interactionMatrix <- NULL  
  
  for (interaction in interactions) {
    col1 <- interaction[[1]]
    col2 <- interaction[[2]]
    startCol1 <- colIndex[col1]
    startCol2 <- colIndex[col2]
    nLevelsCol1 <- nlevels[col1] - 1
    nLevelsCol2 <- nlevels[col2] - 1
    
    for (level1 in 0:(nLevelsCol1 - 1)) {
      for (level2 in 0:(nLevelsCol2 - 1)) {
        interactionColumn <- X[, startCol1 + level1] * X[, startCol2 + level2]
        interactionMatrix <- cbind(interactionMatrix, interactionColumn)
      }
    }
  }
  
  interactionModel <- cbind(X, interactionMatrix)
  return(interactionModel)
}

data_robust = readMat("D:/博士项目/partial profile/2025.1.26/case_study/design/design_robust.mat")
data_true = readMat("D:/博士项目/partial profile/2025.1.26/case_study/design/design_true.mat")
data_main = readMat("D:/博士项目/partial profile/2025.1.26/case_study/design/design_main.mat")
data_org = read_excel("D:/博士项目/partial profile/2025.1.26/case_study/design/X_org.xlsx", col_names = FALSE)


nlevels = c(2,3,3,3,3,3,5)
interactions <- list(c(1, 4), c(1, 7))
nresps <- 100
ndraws <- 500
#14cs 3alt
nalts <- 2
nsets <- 42
nruns <- nsets*nalts
nruns_group <- nsets*nalts/3
betavec <- c(-0.4,-0.5,0,-0.4,0.1,-0.8,0,-0.5,0,-0.5,0.2,-0.5,-0.25,0,0.25,-0.0431, 0.0345, 0.012,-0.0676,-0.048,0.1103)

X_main_robust = transform_X(data_robust$design.robust,nlevels)
X_robust = transform_X_inter(X_main_robust, nlevels, interactions)
X_robust_1 = X_robust[1:28,]
X_robust_2 = X_robust[29:56,]
X_robust_3 = X_robust[57:84,]
X_main_main = transform_X(data_main$design.main,nlevels)
X_main = transform_X_inter(X_main_main, nlevels, interactions)
X_main_1 = X_main[1:28,]
X_main_2 = X_main[29:56,]
X_main_3 = X_main[57:84,]
X_main_true = transform_X(data_true$design.true,nlevels)
X_true = transform_X_inter(X_main_true, nlevels, interactions)
X_true_1 = X_true[1:28,]
X_true_2 = X_true[29:56,]
X_true_3 = X_true[57:84,]
x_main_org = transform_X(as.matrix(data_org),nlevels)
X_org = transform_X_inter(x_main_org, nlevels, interactions)
X_org_1 = X_org[1:28,]
X_org_2 = X_org[29:56,]
X_org_3 = X_org[57:84,]
nr <- nresps*nruns
nrcs <- nr/nalts
nparam = length(betavec)


utility_robust <- matrix(0, nrow = nr, ncol = 1)
utility_true <- matrix(0, nrow = nr, ncol = 1)
utility_org <- matrix(0, nrow = nr, ncol = 1)
utility_main <- matrix(0, nrow = nr, ncol = 1)
for(i in 1:100) {
  utility_robust[((i-1)*nruns_group + 1):(i*nruns_group)] <- X_robust_1%*%betavec
  utility_main[((i-1)*nruns_group + 1):(i*nruns_group)] <- X_main_1%*%betavec
  utility_true[((i-1)*nruns_group + 1):(i*nruns_group)] <- X_true_1%*%betavec
  utility_org[((i-1)*nruns_group + 1):(i*nruns_group)] <- X_org_1%*%betavec
}
for(i in 101:200) {
  utility_robust[((i-1)*nruns_group + 1):(i*nruns_group)] <- X_robust_2%*%betavec
  utility_main[((i-1)*nruns_group + 1):(i*nruns_group)] <- X_main_2%*%betavec
  utility_true[((i-1)*nruns_group + 1):(i*nruns_group)] <- X_true_2%*%betavec
  utility_org[((i-1)*nruns_group + 1):(i*nruns_group)] <- X_org_2%*%betavec
}
for(i in 201:300) {
  utility_robust[((i-1)*nruns_group + 1):(i*nruns_group)] <- X_robust_3%*%betavec
  utility_main[((i-1)*nruns_group + 1):(i*nruns_group)] <- X_main_3%*%betavec
  utility_true[((i-1)*nruns_group + 1):(i*nruns_group)] <- X_true_3%*%betavec
  utility_org[((i-1)*nruns_group + 1):(i*nruns_group)] <- X_org_3%*%betavec
}
exputility_robust <- exp(utility_robust)
exputility_true <- exp(utility_true)
exputility_org <- exp(utility_org)
exputility_main <- exp(utility_main)
prob_main = rep(0, nr)
for (i in 1:nrcs) {
  u1 <- exputility_main[2*i-1,]  
  u2 <- exputility_main[2*i,]    
  
  p1 <- u1 / (u1 + u2)  
  prob_main[2*i-1] = p1
  prob_main[2*i] = 1-p1
}
prob_robust = rep(0, nr)

for (i in 1:nrcs) {
  u1 <- exputility_robust[2*i-1,]  
  u2 <- exputility_robust[2*i,]    
  
  p1 <- u1 / (u1 + u2)  
  prob_robust[2*i-1] = p1
  prob_robust[2*i] = 1-p1
}

prob_true = rep(0, nr)

for (i in 1:nrcs) {
  u1 <- exputility_true[2*i-1,]  
  u2 <- exputility_true[2*i,]    
  
  p1 <- u1 / (u1 + u2)  
  prob_true[2*i-1] = p1
  prob_true[2*i] = 1-p1
}

prob_org = rep(0, nr)

for (i in 1:nrcs) {
  u1 <- exputility_org[2*i-1,]  
  u2 <- exputility_org[2*i,]    
  
  p1 <- u1 / (u1 + u2)  
  prob_org[2*i-1] = p1
  prob_org[2*i] = 1-p1
}

#all data
nminresp <- nresps - 1
Xtotal_robust <- X_robust_1
Xtotal_main <- X_main_1
Xtotal_true <- X_true_1
Xtotal_org <- X_org_1
for(i in 1:nminresp){
  Xtotal_robust <- rbind(Xtotal_robust, X_robust_1) 
  Xtotal_main <- rbind(Xtotal_main, X_main_1)
  Xtotal_true <- rbind(Xtotal_true, X_true_1)
  Xtotal_org <- rbind(Xtotal_org, X_org_1)
}
for(i in 1:nresps ){
  Xtotal_robust <- rbind(Xtotal_robust, X_robust_2) 
  Xtotal_main <- rbind(Xtotal_main, X_main_2)
  Xtotal_true <- rbind(Xtotal_true, X_true_2)
  Xtotal_org <- rbind(Xtotal_org, X_org_2)
}
for(i in 1:nresps ){
  Xtotal_robust <- rbind(Xtotal_robust, X_robust_3) 
  Xtotal_main <- rbind(Xtotal_main, X_main_3)
  Xtotal_true <- rbind(Xtotal_true, X_true_3)
  Xtotal_org <- rbind(Xtotal_org, X_org_3)
}
#wide table
Xtotal_wide_robust <- matrix(nrow = nsets*nresps, ncol = 42)
for (i in 1:(nsets*nresps)) {
  Xtotal_wide_robust[i, ] <- c(Xtotal_robust[2*i-1, ], Xtotal_robust[2*i, ])
}

Xtotal_wide_true <- matrix(nrow = nsets*nresps, ncol = 42)
for (i in 1:(nsets*nresps)) {
  Xtotal_wide_true[i, ] <- c(Xtotal_true[2*i-1, ], Xtotal_true[2*i, ])
}

Xtotal_wide_org <- matrix(nrow = nsets*nresps, ncol = 42)
for (i in 1:(nsets*nresps)) {
  Xtotal_wide_org[i, ] <- c(Xtotal_org[2*i-1, ], Xtotal_org[2*i, ])
}

Xtotal_wide_main <- matrix(nrow = nsets*nresps, ncol = 42)
for (i in 1:(nsets*nresps)) {
  Xtotal_wide_main[i, ] <- c(Xtotal_main[2*i-1, ], Xtotal_main[2*i, ])
}

colnames(Xtotal_wide_robust) = c("main1_A","main2_A","main3_A","main4_A","main5_A", "main6_A","main7_A","main8_A","main9_A","main10_A","main11_A","main12_A","main13_A","main14_A","main15_A","int1_A","int2_A","int3_A","int4_A","int5_A","int6_A","main1_B","main2_B","main3_B","main4_B","main5_B",         "main6_B","main7_B","main8_B","main9_B","main10_B","main11_B","main12_B","main13_B","main14_B","main15_B","int1_B","int2_B","int3_B","int4_B","int5_B","int6_B")
Xtotal_wide_robust = as.data.frame(Xtotal_wide_robust)
Xtotal_wide_robust$subject_id = rep(1:300,each = 14)
Xtotal_wide_robust$cs_id = c(rep(c(1:14),time = 100),rep(c(15:28),time = 100),rep(c(29:42),time = 100))


colnames(Xtotal_wide_main) = c("main1_A","main2_A","main3_A","main4_A","main5_A", "main6_A","main7_A","main8_A","main9_A","main10_A","main11_A","main12_A","main13_A","main14_A","main15_A","int1_A","int2_A","int3_A","int4_A","int5_A","int6_A","main1_B","main2_B","main3_B","main4_B","main5_B",         "main6_B","main7_B","main8_B","main9_B","main10_B","main11_B","main12_B","main13_B","main14_B","main15_B","int1_B","int2_B","int3_B","int4_B","int5_B","int6_B")
Xtotal_wide_main = as.data.frame(Xtotal_wide_main)
Xtotal_wide_main$subject_id = rep(1:300,each = 14)
Xtotal_wide_main$cs_id = c(rep(c(1:14),time = 100),rep(c(15:28),time = 100),rep(c(29:42),time = 100))

colnames(Xtotal_wide_true) = c("main1_A","main2_A","main3_A","main4_A","main5_A", "main6_A","main7_A","main8_A","main9_A","main10_A","main11_A","main12_A","main13_A","main14_A","main15_A","int1_A","int2_A","int3_A","int4_A","int5_A","int6_A","main1_B","main2_B","main3_B","main4_B","main5_B",         "main6_B","main7_B","main8_B","main9_B","main10_B","main11_B","main12_B","main13_B","main14_B","main15_B","int1_B","int2_B","int3_B","int4_B","int5_B","int6_B")
Xtotal_wide_true = as.data.frame(Xtotal_wide_true)
Xtotal_wide_true$subject_id = rep(1:300,each = 14)
Xtotal_wide_true$cs_id = c(rep(c(1:14),time = 100),rep(c(15:28),time = 100),rep(c(29:42),time = 100))

colnames(Xtotal_wide_org) = c("main1_A","main2_A","main3_A","main4_A","main5_A",  "main6_A","main7_A","main8_A","main9_A","main10_A","main11_A","main12_A","main13_A","main14_A","main15_A","int1_A","int2_A","int3_A","int4_A","int5_A","int6_A","main1_B","main2_B","main3_B","main4_B","main5_B",         "main6_B","main7_B","main8_B","main9_B","main10_B","main11_B","main12_B","main13_B","main14_B","main15_B","int1_B","int2_B","int3_B","int4_B","int5_B","int6_B")
Xtotal_wide_org = as.data.frame(Xtotal_wide_org)
Xtotal_wide_org$subject_id = rep(1:300,each = 14)
Xtotal_wide_org$cs_id = c(rep(c(1:14),time = 100),rep(c(15:28),time = 100),rep(c(29:42),time = 100))


#robust
data_robust <- vector("list", ndraws)

for (b in 1:ndraws) {
  start_time <- Sys.time()
  Choice_robust <- rep(0, nrcs)  
  for (i in 1:nrcs) {
    p1 <- prob_robust[2*i - 1]  
    rnd <- runif(1)
    if (p1 >= rnd) {
      Choice_robust[i] <- 'A'
    } else {
      Choice_robust[i] <- 'B'
    }
  }
  
  # Combine the design matrix with the simulated choice
  data_choice <- cbind(Xtotal_wide_robust, Choice_robust)
  # Convert to mlogit data
  data_z <- mlogit.data(data_choice, 
                        choice = "Choice_robust", 
                        shape  = "wide", 
                        varying = 1:42, 
                        sep    = "_", 
                        id.var = "subject_id")
  
  
  # Store this data frame (or matrix) in the list
  data_robust[[b]] <- data_z
  iteration_time <- Sys.time() - start_time
  cat("draws:", b, "Time:", iteration_time, "\n")
}


#true
data_true <- vector("list", ndraws)

for (b in 1:ndraws) {
  start_time <- Sys.time()
  Choice_true <- rep(0, nrcs)  
  for (i in 1:nrcs) {
    p1 <- prob_true[2*i - 1]  
    rnd <- runif(1)
    if (p1 >= rnd) {
      Choice_true[i] <- 'A'
    } else {
      Choice_true[i] <- 'B'
    }
  }
  
  # Combine the design matrix with the simulated choice
  data_choice <- cbind(Xtotal_wide_true, Choice_true)
  # Convert to mlogit data
  data_z <- mlogit.data(data_choice, 
                        choice = "Choice_true", 
                        shape  = "wide", 
                        varying = 1:42, 
                        sep    = "_", 
                        id.var = "subject_id")
  
  
  # Store this data frame (or matrix) in the list
  data_true[[b]] <- data_z
  iteration_time <- Sys.time() - start_time
  cat("draws:", b, "Time:", iteration_time, "\n")
}




#org
data_org <- vector("list", ndraws)

for (b in 1:ndraws) {
  start_time <- Sys.time()
  Choice_org <- rep(0, nrcs)  
  for (i in 1:nrcs) {
    p1 <- prob_org[2*i - 1]  
    rnd <- runif(1)
    if (p1 >= rnd) {
      Choice_org[i] <- 'A'
    } else {
      Choice_org[i] <- 'B'
    }
  }
  
  # Combine the design matrix with the simulated choice
  data_choice <- cbind(Xtotal_wide_org, Choice_org)
  # Convert to mlogit data
  data_z <- mlogit.data(data_choice, 
                        choice = "Choice_org", 
                        shape  = "wide", 
                        varying = 1:42, 
                        sep    = "_", 
                        id.var = "subject_id")
  
  
  # Store this data frame (or matrix) in the list
  data_org[[b]] <- data_z
  iteration_time <- Sys.time() - start_time
  cat("draws:", b, "Time:", iteration_time, "\n")
}


#main
data_main <- vector("list", ndraws)

for (b in 1:ndraws) {
  start_time <- Sys.time()
  Choice_main <- rep(0, nrcs)  
  for (i in 1:nrcs) {
    p1 <- prob_main[2*i - 1]  
    rnd <- runif(1)
    if (p1 >= rnd) {
      Choice_main[i] <- 'A'
    } else {
      Choice_main[i] <- 'B'
    }
  }
  
  # Combine the design matrix with the simulated choice
  data_choice <- cbind(Xtotal_wide_main, Choice_main)
  # Convert to mlogit data
  data_z <- mlogit.data(data_choice, 
                        choice = "Choice_main", 
                        shape  = "wide", 
                        varying = 1:42, 
                        sep    = "_", 
                        id.var = "subject_id")
  
  
  # Store this data frame (or matrix) in the list
  data_main[[b]] <- data_z
  iteration_time <- Sys.time() - start_time
  cat("draws:", b, "Time:", iteration_time, "\n")
}


bhatCov_robust     <- matrix(0, nrow = 21, ncol = ndraws)
for (b in 1:ndraws) {
  start_time <- Sys.time()
  data <- data_robust[[b]]
  fit <- mlogit(Choice_robust ~ main1 + main2 + main3 + main4 + main5 +main6 + main7 + main8 + main9 + main10 + 
                  main11 + main12 + main13 + main14 + main15 +int1 + int2 + int3 + int4 + int5 + int6 | 0, 
                data  = data)
  
  # Extract the coefficients
  bhatCov_robust[, b]       <- fit[["coefficients"]][1:21]
  iteration_time <- Sys.time() - start_time
  cat("draws:", b, "Time:", iteration_time, "\n")
}
write.csv(bhatCov_robust, "MNL_beta_robust.csv", row.names = FALSE)

bhatCov_true     <- matrix(0, nrow = 21, ncol = ndraws)
for (b in 1:ndraws) {
  start_time <- Sys.time()
  data <- data_true[[b]]
  fit <- mlogit(Choice_true ~ main1 + main2 + main3 + main4 + main5 +main6 + main7 + main8 + main9 + main10 + 
                  main11 + main12 + main13 + main14 + main15 +int1 + int2 + int3 + int4 + int5 + int6 | 0, 
                data  = data)
  
  # Extract the coefficients
  bhatCov_true[, b]       <- fit[["coefficients"]][1:21]
  iteration_time <- Sys.time() - start_time
  cat("draws:", b, "Time:", iteration_time, "\n")
}
write.csv(bhatCov_true, "MNL_beta_true.csv", row.names = FALSE)

bhatCov_org     <- matrix(0, nrow = 21, ncol = ndraws)
for (b in 1:ndraws) {
  start_time <- Sys.time()
  data <- data_org[[b]]
  fit <- mlogit(Choice_org ~ main1 + main2 + main3 + main4 + main5 +main6 + main7 + main8 + main9 + main10 + 
                  main11 + main12 + main13 + main14 + main15 +int1 + int2 + int3 + int4 + int5 + int6 | 0, 
                data  = data)
  
  # Extract the coefficients
  bhatCov_org[, b]       <- fit[["coefficients"]][1:21]
  iteration_time <- Sys.time() - start_time
  cat("draws:", b, "Time:", iteration_time, "\n")
}
write.csv(bhatCov_org, "MNL_beta_org.csv", row.names = FALSE)

bhatCov_main     <- matrix(0, nrow = 21, ncol = ndraws)
for (b in 1:ndraws) {
  start_time <- Sys.time()
  data <- data_main[[b]]
  fit <- mlogit(Choice_main ~ main1 + main2 + main3 + main4 + main5 +main6 + main7 + main8 + main9 + main10 + 
                  main11 + main12 + main13 + main14 + main15 +int1 + int2 + int3 + int4 + int5 + int6 | 0, 
                data  = data)
  
  # Extract the coefficients
  bhatCov_main[, b]       <- fit[["coefficients"]][1:21]
  iteration_time <- Sys.time() - start_time
  cat("draws:", b, "Time:", iteration_time, "\n")
}
write.csv(bhatCov_main, "MNL_beta_main.csv", row.names = FALSE)
