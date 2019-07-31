# Installing needed packages

install.packages("BiocManager")
library("BiocManager")
BiocManager::install(c("limma", "affy","ecoli2.db","annotate","IRanges"))
library("limma")
library("affy")
library("annotate")
library ("ecoli2.db")
library("readxl")
library("data.table")
library("factoextra")
library("dplyr")
library("igraph")
library("tidygraph")
library("ggraph")




setwd("C:/Users/mexicore/Dropbox/Universidad_pasada/Máster/Asignaturas/TFM/Tablas")

# Data reading

celpath = "C:/Users/mexicore/Dropbox/Universidad_pasada/Máster/Asignaturas/TFM/Datos/Archivos CEL/Q&R/Q&R-X-A/"

data = ReadAffy(celfile.path=celpath)

ph <- data@phenoData

# Normalization

data.rma = rma(data)
data.matrix = exprs(data.rma)
plotMDS(data.rma)

# Quality control

## Images

for (i in 1:2)
{
  name = paste("Q-X-GL",i,".jpg",sep="")
  jpeg(name)
  image(data[,i],main="Q-X-GL")
  dev.off()
}

for (i in 3:4)
{
  name = paste("R-X-GL",i,".jpg",sep="")
  jpeg(name)
  image(data[,i],main="R-X-GL")
  dev.off()
}


## Dendograms

par(mfrow=c(1,2))

### Before normalization

eset <-exprs(data)
distance <-dist(t(eset), method = "maximum")
clusters <-hclust(distance)
plot(clusters, main="Q&R-X-Gl - Before normalization",labels = c("Q1","Q2","R1","R2"))

### After normalization

eset_rma <-exprs(data.rma)
distance_rma <-dist(t(eset_rma), method = "maximum")
clusters_rma <-hclust(distance_rma)
plot(clusters_rma,     main="Q&R-X-Gl -After normalization", labels = c("Q1","Q2","R1","R2"))

## Boxplot

par(mfrow=c(1,2))

### Before normalization
boxplot(data,which='pm', names = c("Q1","Q2","R1","R2"),
        target='core', main="Q&R-X-Gl - Before normalization")

### Después de normalización
boxplot(data.matrix,name = c("Q1","Q2","R1","R2"),main="Q&R-X-Gl - After normalization")

## Histogram

### Antes de normalización
hist(data,lwd=2,which='pm',ylab='Density',xlab='Log2 intensities', main='Q&R-X-Gl - Before normalization', 
     target='core')
legend('topright',  legend=c("Q1","Q2","R1","R2"), col = c("green","blue","red","black") ,  lty=1, 
       bty='n', cex=.75)

### After normalization
plotDensities(data.rma,main="Q&R-X-Gl - After normalization", legend= F)
legend('topright',  legend=c("Q1","Q2","R1","R2"), col = c("green","blue","red","black"),lty=1, bty='n', cex=.75)


# Annotation

ID <- featureNames(data.rma)

Symbol <- getSYMBOL(ID, "ecoli2.db")
Name <- as.character(lookUp(ID, "ecoli2.db", "GENENAME"))
Ensembl <- as.character(lookUp(ID, "ecoli2.db", "ENSEMBL"))

# Experimental design

## We asign experimental groups to the measures and we build the experimental design

ph@data[ ,2] = c("Q","Q","R","R")
colnames(ph@data)[2]="source"
ph@data[ ,3] = c("uno","dos","uno","dos")
colnames(ph@data)[3]="experiment"

ph@data


groups = ph@data$source

experiment = ph@data$experiment

f = factor(groups,levels=c("Q","R"))
e = factor(experiment, levels =c("uno","dos"))

design = model.matrix(~0+f+e)

colnames(design) = c("Q","R","Exp")

# We fit the experimental design to a 

data.fit = lmFit(data.matrix,design)
data.fit$coefficients[1:4,]

x <-c("R-Q", "Exp")
contrast.matrix = makeContrasts(contrasts = x, levels=design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)
names(data.fit.eb)

results <- decideTests(data.fit.eb)

summary(results)

## Obtener genes diferencialmente expresados 

options(digits=2)
tab = topTable(data.fit.eb,coef=1,number=20000,adjust.method="BH")
topgenes = tab[tab[, "P.Value"] < 0.01, ]
dim(topgenes)
topups = topgenes[topgenes[, "logFC"] > 0.01, ]
dim(topups)
topdowns = topgenes[topgenes[, "logFC"] < -0.01, ]
dim(topdowns)

IDs.up = rownames(topups)
IDs.down = rownames(topdowns)
IDs.up
IDs.down

#

DB<-read_excel("C:/Users/mexicore/Downloads/jbc.M115.655688-2.xlsx")
DB<-data.frame(DB[, 1:3])

IDs.up <- data.frame(IDs.up)
IDs.down <- data.frame(IDs.down)

DB.IDsup <- na.omit(DB[match(IDs.up[,1],DB[,1]),])
DB.IDsdown <- na.omit(DB[match(IDs.down[,1],DB[,1]),])

write.csv(DB.IDsup, file="IDSUPDEF1")
write.csv(DB.IDsdown,file="IDSDOWNDEF1")

data.matrix.up = data.frame(data.matrix[(rownames(topups)),])
data.matrix.down = data.matrix[(rownames(topdowns)),]

matrix <- data.frame(cbind(rownames(data.matrix),data.matrix))

up.a.x<-read_excel("C:/Users/mexicore/Dropbox/Universidad_pasada/Máster/Asignaturas/TFM/Tablas/Genes_TFM.xlsx", 
                    sheet="UP-X-A-CRP")

Dt1 <- data.table(up.a.x, key ="Column.ID")
Dt2 <- data.table (matrix, key = "V1")

cluster.up <- data.frame(Dt1[Dt2, nomatch = 0])
rownames(cluster.up)<- cluster.up[,2]
cluster.up<- data.frame(cluster.up[,-c(1,2,3,4)])

cluster.up1 <- lapply(cluster.up, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
cluster.up1 <- data.frame(cluster.up1)
rownames(cluster.up1)<-rownames(cluster.up)


# Non hierarchical clustering

## K-means

### Number of centers

fviz_nbclust(x = cluster.up1, FUNcluster = kmeans, method = "wss", k.max = 10, 
             diss = get_dist(cluster.up1, method = "euclidean"), nstart = 50)

### Plotting

set.seed(123)
km_clusters <- kmeans(x = cluster.up1, centers = 3, nstart = 50)


fviz_cluster(object = km_clusters, data = cluster.up1, show.clust.cent = TRUE,
             ellipse.type = "euclid", star.plot = TRUE, repel = TRUE) +
  labs(title = "Clustering K-means - Regulados positivamente X-A") +
  theme_bw() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))



## K-medoids clustering (PAM)

### Number of clusters

library(cluster)
fviz_nbclust(x = cluster.up1, FUNcluster = pam, method = "wss", k.max = 10,
             diss = dist(cluster.up1, method = "manhattan"))

set.seed(123)
pam_clusters <- pam(x = cluster.up1, k = 3, metric = "manhattan")

fviz_cluster(object = pam_clusters, data = cluster.up1, ellipse.type = "t",
             repel = TRUE) +
  theme_bw() +
  labs(title = "Clustering PAM - Regulados positivamente X-A") +
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))


# Hierarchical clustering 

## Hierarchical K-means

# Number of centers

set.seed(123)
hc_euclidea_completo <- hclust(d = dist(x = cluster.up1, method = "euclidean"),
                               method = "complete")
fviz_dend(x = hc_euclidea_completo, cex = 0.5, main = "Clustering jerárquico por linkage completo - Regulados positivamente X-A",
          sub = "Distancia euclídea", k= 4, rect = TRUE) +
  theme(plot.title =  element_text(hjust = 0.5, size = 15))

# Plotting

hkmeans_cluster <- hkmeans(x = cluster.up1, hc.metric = "euclidean",
                           hc.method = "complete", k = 4)

fviz_cluster(object = hkmeans_cluster, pallete = "jco", repel = TRUE) +
  theme_bw() + labs(title = "K-means Jerarquico - Regulados positivamente X-GL")+  theme(plot.title =  element_text(hjust = 0.5, size = 15))

# Hierarchical clustering

mat_dist <- dist(x = cluster.up1, method = "euclidean")

# Dendrograms with linkage complete y average

hc_euclidea_complete <- hclust(d = mat_dist, method = "complete")
hc_euclidea_average  <- hclust(d = mat_dist, method = "average")
cor(x = mat_dist, cophenetic(hc_euclidea_complete))
cor(x = mat_dist, cophenetic(hc_euclidea_average))

plot(hc_euclidea_average)

plot(hc_euclidea_complete)


# Data joining

out.up <- data.frame(cbind(rownames(cluster.up1), km_clusters$cluster,pam_clusters$cluster, cutree(hc_euclidea_completo, 4)))

DB<-read_excel("C:/Users/mexicore/Downloads/jbc.M115.655688-2.xlsx")
DB<-data.frame(DB[, 1:3])

Dt1 <- data.table(DB, key ="Gene.Symbol")
Dt2 <- data.table (out.up, key = "X1")

cluster.def <- data.frame(Dt1[Dt2, nomatch = 0])
colnames(cluster.def)<-colnames<-c("Transcript.ID","Gene.Symbol","Gene.Title", "Kmeans cluster", "PAM cluster", "Hierarchical cluster")

write.csv(cluster.def, file="clusters-up-x-a-crp")




