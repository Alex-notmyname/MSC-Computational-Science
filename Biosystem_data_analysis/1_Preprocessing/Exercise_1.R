counts = read.csv('counts.csv')
head(counts)

totalcounts = colSums(counts)

ratespermillion = as.data.frame(sweep(10^6*counts,2,totalcounts,FUN="/"))

ct.names = paste('ct',1:4,sep = '.')
kd.names = paste('kd',1:3,sep = '.')

control.mean = apply(ratespermillion[ct.names],1,mean)
control.sd =  apply(ratespermillion[ct.names],1,sd)
knockdown.mean = apply(ratespermillion[kd.names],1,mean)
knockdown.sd =  apply(ratespermillion[kd.names],1,sd)

plot(control.mean,control.sd,log = 'xy',pch = 20,cex=.5,col='blue',
     xlab='Mean Rate',ylab='Stand dev')
points(knockdown.mean, knockdown.sd, pch=20, cex=.5, col='red')
title('Standard deviation vs mean rate')

all.fit = lm(sd~mean,
          data.frame(mean=log(c(control.mean, knockdown.mean)), sd=
                     log(c(control.sd, knockdown.sd))))

coef(all.fit)

t.test(ratespermillion['FBgn0000043', ct.names], ratespermillion['FBgn0000043',
       kd.names], var.equal = FALSE)

pvalues.orig = apply(ratespermillion, 1, function(row){t.test(x=row[ct.names], 
                     y=row[kd.names], var.equal=FALSE)$p.value})


# Applying the variance-stabilizing transformation
# g(x) = x^(1-b)
transformed = as.data.frame(ratespermillion^(1-coef(all.fit)['mean']))
head(transformed)

#calculate the mean and std for variance-stabilized data

control.t.mean = apply(transformed[ct.names], 1, mean)
control.t.sd = apply(transformed[ct.names], 1, sd)
knockdown.t.mean = apply(transformed[kd.names], 1, mean)
knockdown.t.sd = apply(transformed[kd.names], 1, sd)

# Plot std vs mean of transformed data
plot(x=control.t.mean, y=control.t.sd, pch=20, cex=.5, col='blue',
     ylim=c(-0.03, 0.3), xlab='Mean transformed value', ylab='transformed std')
points(x=knockdown.t.mean, y=knockdown.t.sd, pch=20, cex=.5, col='red')
title('STD versus Mean of transformed value')

# t-test on tranformed data
t.test(x=transformed['FBgn0000043', ct.names], y=transformed['FBgn0000043',
       kd.names], var.equal=TRUE)

# T-test on all rows, and extracting the p-values from the test results
pvalues.t <- apply(transformed, 1, function(row){t.test(x=row[ct.names], 
                   y=row[kd.names], var.equal=TRUE)$p.value})

# Add pvalues.t to ratespermillion table
ratespermillion$pvalue.t <- pvalues.t
head(ratespermillion)

# Compare the p-values before and after transformation
plot(x=-log10(pvalues.orig), y=-log10(pvalues.t), xlim=c(0,8), ylim=c(0,8),
     xlab='-log10(p-values) before transformation', ylab=
       '-log10(p-values) after transformation')
abline(coef=c(0,1), col='blue')
legend(0.2, 7.5, legend=c('x=y'), col=c('blue'), lty=1, box.lty=0)
title('p-values before and after transformation')

# Calculate the log2ratio
log2ratio <- apply(counts, 1, function(x){log(sum(x[kd.names])/
                   sum(totalcounts[kd.names]), 2) - log(sum(x[ct.names])/
                   sum(totalcounts[ct.names]), 2)})
head(log2ratio)

# Combining log2ratios and p-values to create volcano plot
par(mfrow=c(1,1))
p05 = -log10(0.05)
plot(y=-log10(pvalues.t),x=log2ratio)
lines(x=c(-1,-1),y=c(p05,7),col='gray')
lines(x=c(1,1),y=c(p05,7),col='gray')
lines(x=c(-6,-1),y=c(p05,p05),col='gray')
lines(x=c(1,4),y=c(p05,p05),col='gray')
text( x = 3, y = p05+0.2, label = "P = 0.05")
text (x=1.15, y = 6, label="FC = 2",srt=90)  
title('volcano plot')

# Plot a histogram of p-values
nbins <- 100
bins <- hist(pvalues.t, nclass=nbins, col='grey50', main='Histogram p-values')

# P-value cut-off
null.bins <- bins$mids > 0.4
null.level <- mean(bins$counts[null.bins])
plot(bins, col='grey50')
abline(a=null.level, b=0, col='red', lwd=2)

# Calculate False discovery rate (FDR)
fdr.05 <- null.level*nbins*0.05 / sum(pvalues.t <= 0.05) 
fdr.05

# Plot the relationship between p-value cut-off and FDR
pvalue.cutoff <- seq(1/nbins, 1, by=1/nbins)
fdr <- vector(mode="numeric",length=length(pvalue.cutoff))
for (i in seq_along(pvalue.cutoff)) {
  fdr[i] <- null.level*(pvalue.cutoff[i]*nbins)/sum(pvalues.t <= pvalue.cutoff[i])
}
plot(pvalue.cutoff, fdr, xlab='p-value cut-off', ylab='FDR', ylim=c(0,1))
title('FDR and p-value relationship')
