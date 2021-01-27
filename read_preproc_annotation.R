# Data reading

# celpath = "C:/Users/mexicore/Dropbox/Universidad_pasada/M?ster/Asignaturas/TFM/Datos/Archivos CEL/Q&R/Q&R-X-A/"

celpath = file.path(root, "cells")
data = ReadAffy(celfile.path=celpath)
ph <- data@phenoData

data_name1 <- "Q-X-GL" 
data_name2 <- "R-X-GL" 

# Normalization

data.rma = rma(data)
data.matrix = exprs(data.rma)
plotMDS(data.rma)

# Quality control

## Images

for (i in 1:2)
{
  name = paste(data_name1,i,".jpg",sep="")
  jpeg(name)
  image(data[,i],main=data_name1)
  dev.off()
}

for (i in 3:4)
{
  name = paste(data_name2,i,".jpg",sep="")
  jpeg(name)
  image(data[,i],main=data_name2)
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

### After normalization
boxplot(data.matrix,name = c("Q1","Q2","R1","R2"),main="Q&R-X-Gl - After normalization")

## Histogram

### Antes de normalizaci?n
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