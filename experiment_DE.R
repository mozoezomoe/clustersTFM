

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

up.a.x<-read_excel("C:/Users/mexicore/Dropbox/Universidad_pasada/M?ster/Asignaturas/TFM/Tablas/Genes_TFM.xlsx", 
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
