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
fviz_dend(x = hc_euclidea_completo, cex = 0.5, main = "Clustering jer?rquico por linkage completo - Regulados positivamente X-A",
          sub = "Distancia eucl?dea", k= 4, rect = TRUE) +
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



