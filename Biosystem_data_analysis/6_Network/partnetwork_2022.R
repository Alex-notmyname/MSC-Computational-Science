# Partial Correlation Network calculation

#install.packages("igraph")
library(igraph)

realadj <- read.csv(file='realadj.csv', header= FALSE)
names <- c("GLCi","G6P","F6P","F16bP","TRIO","BPG","3PGA","2PGA","PEP","PYR","ACE","P","NADH")
names_P <- c("GLCi","G6P","F6P","F16bP","TRIO","BPG","3PGA","2PGA","PEP","PYR","ACE","NADH")
# build the graph
realadj_P <- as.matrix(realadj[-12,-12])
colnames(realadj_P) <- rownames(realadj_P) <- names_P
real_net_P <- graph_from_adjacency_matrix(realadj_P, mode = 'undirected')
# calculate a layout to keep
real_layout_P <- layout_nicely(real_net_P)
# plot the graph
plot(real_net_P,layout=real_layout_P)

## 
load("networkdata.Rdata") # load data
dataset <- data_n10  # Here you can fill in the corresponding data frame

dims <- dim(dataset)
I <- dims[1]
J <- dims[2]

### Calculate partial correlation via regression SLOW APPROACH
dataset <- as.matrix(dataset)

partcor2 <- sapply(1:J,function(i){
  sapply(1:J,function(j){
    set <- 1:J
    s <- set[-c(i,j)]
    fit_i <- lm(dataset[,i] ~ as.matrix(dataset[,s]))
    fit_j <- lm(dataset[,j] ~ as.matrix(dataset[,s]))
    cor(residuals(fit_i),residuals(fit_j))
  })
})

# calculate Partial Correlation via inverse FAST APPROACH
cormat <- cor(dataset, method="pearson")
w <- solve(cormat)
denom <- sqrt(diag(w)) %*% t(sqrt(diag(w))) 
partcor <- -w/denom

## We will write this into a function
pcc <- function(dat){
  w <- solve(cor(dat, method="pearson"))
  -w/(sqrt(diag(w)) %*% t(sqrt(diag(w))))
}


## Start permutation of data and calculate Partial correlation of permuted data
counting <- matrix(0,J,J)
partcor_perm <- matrix(0,J,J)
Conf_threshold <- 0.01
Nperm <- 1000
for (p in 1:Nperm) {
  print(p)
  data_perm <- matrix(0,I,J)
  for (j in 1:J) { 
    data_perm[,j] <- dataset[sample(1:I,size=I),j] 
  } # data_perm is permuted dataset
  
  partcor_perm <- pcc(data_perm)    
  counting[abs(partcor_perm)>abs(partcor)] <- counting[abs(partcor_perm)>abs(partcor)]+1
}
adj_part_pert <- counting < Conf_threshold*Nperm
diag(adj_part_pert) <- 0

adj_part_pert_P = adj_part_pert[-12,-12]
colnames(adj_part_pert_P) <- rownames(adj_part_pert_P) <- names_P
net_part_pert_P <- graph_from_adjacency_matrix(adj_part_pert_P, mode = 'undirected')

# plot the graph
plot(net_part_pert_P, main = 'Partial correlation p<0.01',layout=real_layout_P)


# below the tp, tn, fp, fn, tpr, tnr and g-score are calculated
source('quality.R')
Q = quality(realadj,adj_part_pert)
Q
