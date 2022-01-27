# INSTALL PACKAGES IF NECESSARY
install.packages("igraph")
library(igraph)

# SET WORKING DIRECTORY HERE using setwd() or in the Files panel

# Exercise 1

# calculates and plot the correlation network 
realadj <- read.csv(file='realadj.csv', header= FALSE)
names <- c("GLCi","G6P","F6P","F16bP","TRIO","BPG","3PGA","2PGA","PEP","PYR","ACE","P","NADH")
colnames(realadj) <- rownames(realadj) <- names
# build the graph
realadj <- as.matrix(realadj)
real_net <- graph_from_adjacency_matrix(realadj, mode = 'undirected')
# plot the graph
# calculate a layout to keep
real_layout <- layout_nicely(real_net)
plot(real_net, main = 'Network including P',layout=real_layout)


# Without P (var 12)
realadj_P <- realadj[-12,-12]
names_P <- names[-12]
colnames(realadj_P) <- rownames(realadj_P) <- names_P 
real_net_P <-graph_from_adjacency_matrix(realadj_P, mode = 'undirected')
# calculate the layout to keep
real_layout_P <- layout_nicely(real_net_P)
# plot the graph
plot(real_net_P, main='Network excluding P',layout=real_layout_P)


## Exercise 2
load("networkdata.Rdata") # load data

# part 1
dataset <- data999  # Here you can fill the specific data frame you want
cormat <- cor(dataset, method="pearson")    #Pearson correlation is calculated.

#part 2
par(mfrow=c(1,3))
plot(dataset[,4],dataset[,3],xlab="F16b",ylab = "F6P")
plot(dataset[,4],dataset[,5],xlab="F16b",ylab = "TRIO")
plot(dataset[,4],dataset[,6],xlab="F16b",ylab = "BPG")

#part 3
# Set threshold
Corr_threshold <- 0.1

# calculate associations if abs(correlation) > threshold
cormat_threshold <- (abs(cormat) > Corr_threshold) * 1

# make diagonal 0 as we do not want each metabolite be connected to itself
diag(cormat_threshold) <- 0

# Make graph without P
cormat_threshold_P <- cormat_threshold[-12,-12]
colnames(cormat_threshold_P) <- rownames(cormat_threshold_P) <- names_P
net_P <- graph_from_adjacency_matrix(cormat_threshold_P, mode = 'undirected')
# plot the graph
plot(net_P, main='Network cor > threshold;excluding P',layout=real_layout_P)

# part 4
## Use permutations to calculate significance
# below the data is permuted and the correlation of the true data is compared with the permuted data.
# if the value of the permuted data is larger than the correlation of the true data, it is counted.
# if the fraction of counts is lower than the cut, it is inferred as a connection
dims <- dim(dataset)
I <- dims[1]
J <- dims[2]
Conf_threshold <- 0.01         # Confidence threshold = 0.01
Nperm <- 1000         # number of permutations
counting <- matrix(0,J,J)
for (p in 1:Nperm) {
  print(p)
  data_perm <- matrix(0,I,J)
  for (j in 1:J) { 
    data_perm[,j] <- dataset[sample(1:I,size=I),j] 
  } # data contains permuted values from dataset
  cor_perm <- cor(data_perm, method="pearson") # calculate correlations in permuted data
  counting[abs(cor_perm)>abs(cormat)] <- counting[abs(cor_perm)>abs(cormat)]+1
}
adjpert <- counting < Conf_threshold*Nperm
diag(adjpert) <- 0

# Make graph without P
adjpert_P <- adjpert[-12,-12]
colnames(adjpert_P) <- rownames(adjpert_P) <- names_P
adj_net_P <- graph_from_adjacency_matrix(adjpert_P, mode = 'undirected')
# plot the graph
plot(adj_net_P, main='Network p < cut ;excluding P',layout=real_layout_P)

# Part 5
# below the true positives, true negatives, false positives, false negatives, true positive rate, true negative rate and g-score are calculated
source('quality.R')
Q <- quality(realadj,adjpert)
Q
