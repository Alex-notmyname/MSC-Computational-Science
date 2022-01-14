data <- read.csv('Cachexia.csv')
head(data[,1:3])

# Make a plot of this data table
intensity <- as.matrix(data)
par(mfrow=c(1,2))
matplot(1:77, intensity, type='l', lwd=1, xlab='individuls', ylab='intensity')
matplot(1:63, t(intensity), type='l', lwd=1, xlab='metabolites', ylab='intensity')

# Make transformation on metabolites data
logdata <- log(intensity)
sqrdata <- sqrt(intensity)

par(mfrow=c(1,3))
hist(intensity[,1])
hist(logdata[,1])
hist(sqrdata[,1])

# Perform PCA on both raw data and transformed data
m = colMeans(intensity)                # calculate mean value of columns
Xm = sweep(intensity, 2, m, FUN='-')   # Substract mean on each columns(centering)
ssqtotal <- sum(Xm*Xm)                 # compute sum of squares
USV <- svd(Xm)                         # Perform SVD on centered data sets

# Note that SVD() in R returns a list with length 3: d(n): single values, u(p*n):
# Unitary matrix, v(n*n): orthogonal matrix

# With SVD, it's super easy to perform PCA: T = US, P = V

T <- USV$u %*% diag(USV$d)             # compute scores (principal components)
P <- USV$v                             # compute loadings (principal directions)

npc <- 10                              # compute 10 principal components
ssq <- 0 * (1:npc)                     # initialize an empty ssq vector

for (i in 1:npc){
  Xtest <- T[,i] %*% t(P[,i])          # calculte ssq_n / ssqtotal for each PC
  ssq[i] <- 100*sum(Xtest*Xtest)/ssqtotal
}
ssqcum <- cumsum(ssq)                  # calculate cumulative ssq

data.frame(ssq=ssq, ssqto=ssqcum)      # Put ssq and ssq cumulative value together

# Now do the same procedure for logdata
mlog= colMeans(logdata)
Xlogm=sweep(logdata,2,mlog,FUN="-")
ssqlogtotal <- sum(Xlogm*Xlogm)
USVlog <- svd(Xlogm)
T_log <- USVlog$u %*% diag(USVlog$d)
P_log <- USVlog$v
ssqlog <- 0 * (1:npc)
for (i in 1:npc){
  Xlog_est  <- T_log[,i] %*% t(P_log[,i])
  ssqlog[i] <- 100 * sum(Xlog_est*Xlog_est)/ssqlogtotal
}
ssqlogcum = cumsum(ssqlog)
data.frame(ssqlog=ssqlog,ssqlogcum=ssqlogcum)

# Now let's take a closer look at scores and residuals of both sets
par(mfrow=c(2,2))
plot(T[1:47,1], T[1:47,2], pch=3, col='blue', xlab='T1', ylab='T2', main=
     'No transformation scores', cex.main=0.8)
points(T[48:77,1], T[48:77,2], pch=1, col='red')
res <- Xm - T[,1:2] %*% t(P[,1:2])          # calculate res
matplot(1:77, res, type='l', lwd=1, xlab='', ylab='', main='No transformation residuals')

plot(T_log[1:47,1], T_log[1:47,2], pch=3, col='blue', xlab='T1_log', ylab=
       'T2_log', main='Transformation scores', cex.main=0.8)
points(T_log[48:77,1], T_log[48:77,2], pch=1, col='red')
res_log <- Xlogm - T_log[,1:2] %*% t(P_log[,1:2])          # calculate res
matplot(1:77, res_log, type='l', lwd=1, xlab='', ylab='', main='Transformation residuals')

# Discuss the obvious outlier in No transformation data (individual 12)
par(mfrow=c(1,1))
matplot(1:63, t(Xm), type='l', col='red', xlab='', ylab='')
lines(1:63, Xm[12,], col='black', lwd=2)

# Perform Cross-validation on log_data
source('RMSECV.R')
RMSE <- RMSECV(Xlogm,10,5)
plot(RMSE, xlab='#PC', ylab='RMSECV')

# Now let's continue on log_data models with only 3 PCs
XPCA <- T_log[,1:3] %*% t(P_log[,1:3])
E <- Xlogm - XPCA
Res_IND <- rowSums(E*E)
Res_VARS <- as.numeric(colSums(E*E))
par(mfrow=c(1,2))
barplot(Res_IND, names.arg=c(1:77), main='residuals per individual', cex.main=0.8)
barplot(Res_VARS, names.arg=c(1:63), main='residuals per variables', cex.main=0.8)

# Try to interpret why there are outlying individuals and variables
par(mfrow=c(1,2))
matplot(1:63, t(E), type='l', lwd=1, col='red', xlab='Individual 5', ylab='Residuals')
lines(1:63,t(E[5,]), lwd=2, col='black')

matplot(1:77, E, type='l', lwd=1, col='red', xlab='Variable 19', ylab='Residuals')
lines(1:77,E[,19], lwd=2, col='black')
