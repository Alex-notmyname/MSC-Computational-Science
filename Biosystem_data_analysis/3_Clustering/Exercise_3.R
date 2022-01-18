# Load packages
library(fossil)
library(kohonen) # used for SOM
library(cluster) # used for pam
library(datasets) # For the Iris data

# Load Iris data
data(iris)
summary(iris)

X = as.matrix(iris[,-5])
Y = as.numeric(iris$Species)

# Make some plots to have an initial intuition on the data
par(mfrow=c(2,2))
plot(iris$Sepal.Length, iris$Sepal.Width, col=iris$Species)
legend(0,7.1, c('Setosa', 'Versicolor', 'Virginica'), col=1:3, pch=1)
plot(iris$Sepal.Length, iris$Petal.Width, col=iris$Species)
plot(iris$Sepal.Width, iris$Petal.Width, col=iris$Species)
plot(iris$Petal.Length, iris$Petal.Width, col=iris$Species)

# Start with k-medoids clustering method
pam.Euc <- pam(X, 3, diss=FALSE)    # Euclidean distance is the default setting
# diss=FALSE means X is a matrix of observations
table(pam.Euc$cluster, iris$Species)

# summary(pam.Euc)

rand.index(pam.Euc$cluster, Y)
adj.rand.index(pam.Euc$cluster, Y)

# Make cosine distance function
CosDist = function(X){
  X = as.matrix(X)
  X.norm = X / sqrt(rowSums(X*X))   # normalization, lead to A/||A||
  Cos.Sim = X.norm %*% t(X.norm)    # calculate angles between all rows in X
  return(Cos.Dist = 1-Cos.Sim)
}

# Compute the cosine distance
Cos_iris = CosDist(X)
pam.Cos <- pam(X, 3, diss=TRUE)
table(pam.Cos$cluster, iris$Species)
rand.index(pam.Cos$cluster,Y)
adj.rand.index(pam.Cos$cluster,Y)
pam.Cos$id.med

"Hierarchical clustering"
par(mfrow=c(1,1))
d_iris <- dist(X, method='euclidean')
hc_iris <- hclust(d_iris, method='average')
plot(hc_iris)

clusterCut <- cutree(hc_iris, 3)
table(clusterCut, iris$Species)

rand.index(clusterCut,Y)
adj.rand.index(clusterCut,Y)

# Hierarchical Clustering Cosine distance
hc_iris_cos <- hclust(as.dist(Cos_iris), method = "average")
plot(hc_iris_cos)
clusterCut_cos <- cutree(hc_iris_cos, 3)
table(clusterCut_cos, iris$Species)
adj.rand.index(clusterCut_cos,Y)

# Let's try clustering with 4 clusters
clusterCut_cos2 <- cutree(hc_iris_cos, 4)
table(clusterCut_cos2, iris$Species)
adj.rand.index(clusterCut_cos2, Y)

"Now test the hierarchical clustering method on random data"
par(mfrow=c(1,1))
Gene_data = read.table('random_genes')
plot(Gene_data)

d_Gene <- dist(Gene_data, method='euclidean')
hc_Gene <- hclust(d_Gene, method='average')
plot(hc_Gene)

"Self Organising Map (SOM)"

# Start with a 5 by 5 grid
c1 = 5
c2 = 5
kohmap <- som(as.matrix(X), grid=somgrid(c1, c2, 'hexagonal'), rlen=100)
par(mfrow=c(1,1))
plot(kohmap, type='changes')


par(mfrow=c(2,2))
counts <- plot(kohmap, type='counts', shape='straight')
plot(kohmap,type='dist.neighbours')
similarities <- plot(kohmap, type="quality", palette.name = terrain.colors)
plot(kohmap, type = "mapping",labels = as.integer(Y), col = as.integer(Y),
     main = "mapping plot")

par(mfrow = c(1,1))
plot(kohmap, type="codes", shape="straight",main="profiles")

"SOM on gene data"
Genes = read.csv('Genes.csv', header = FALSE)
G = Genes[,1:13]
G2 = log2(G)
matplot(t(G2), type='l')
Labels = rep(1:10, times=c(100,162,38,43,7,39,19,65,19,25))

aG2 = aggregate(G2, list(Labels), mean)
matplot(t(aG2[,c(-1,-2)]), type="l", xlim = c(0,14),col=1:10, lty = 1:10)
legend(12,2.5, legend = 1:10, col=1:10, lty = 1:10)

#Now train the SOM
c1 = 10
c2 = 7
Genemap <- som(as.matrix(G2[,2:13]), grid=somgrid(c1, c2, 'hexagonal'), rlen=200)
plot(Genemap, type='change')
counts <- plot(Genemap, type='counts', shape='straight')
similarities <- plot(Genemap, type='quality', palette.name=terrain.colors)

plot(Genemap, type="mapping",
     labels = as.integer(Labels), col = as.integer(Labels),
     main = "mapping plot")
plot(Genemap, type="codes", codeRendering="lines", shape="straight",main="profiles")


"DO IT YOURSELF"
# Make a hierarchical clustering using dist and hclust functions
d_gene <- dist(G2, method='euclidean')
hclust_gene <- hclust(d_gene, method='average')
par(mfrow=c(1,1))
plot(hclust_gene, ylab='Height')
