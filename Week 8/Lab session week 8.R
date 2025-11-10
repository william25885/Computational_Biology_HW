# ==============================
# 0) Packages
# ==============================
library(haven)
library(dplyr)
library(cluster)
library(factoextra)
library(janitor)

# ==============================
# 1) Load and merge NHANES data
# ==============================
data_dir <- "data_raw"

demo <- read_xpt(file.path(data_dir,"DEMO_L.XPT")) %>% clean_names()
bpx  <- read_xpt(file.path(data_dir,"BPXO_L.XPT")) %>% clean_names()
bmx  <- read_xpt(file.path(data_dir,"BMX_L.XPT"))  %>% clean_names()

dat <- demo %>%
  inner_join(bmx, by="seqn") %>%
  inner_join(bpx, by="seqn")

# ==============================
# 2) Select variables
# ==============================
set.seed(123)
vars <- dat %>%
  dplyr::select(bmxbmi, bmxwaist, ridageyr, bpxosy1, bpxodi1) %>%
  na.omit() %>%
  sample_n(15)


# ==============================
# 3) Distance / Similarity
# ==============================
d <- dist(scale(vars), method="euclidean")
dist(scale(vars), method="manhattan")

dist(scale(vars), diag = FALSE, upper = FALSE)  # 預設，只顯示下三角 (By default, only the lower triangle is displayed.)
dist(scale(vars), diag = TRUE)                # 顯示對角線 (Show Diagonal)
dist(scale(vars), upper = TRUE)               # 顯示上三角 (Show upper triangle)
dist(scale(vars), diag=TRUE, upper=TRUE)    # 完整矩陣 (Complete Matrix)


library(pheatmap)
pheatmap(dist(scale(vars), method = "euclidean"), 
         cluster_rows = T, cluster_cols = T, 
         main = "Distance Heatmap")

# ==============================
# 4) Mahalanobis distance
# ==============================
mu <- colMeans(vars)
S  <- cov(vars)
 
  #每個人到全體平均得距離 (The average distance from each person to the whole group)
  #用途：常見在統計檢測「異常值 (outlier detection)」，因為離群點會和平均距離特別大。(Usage: Commonly used in statistical detection of outliers, because outliers are very far from the mean.)
mahal <- mahalanobis(vars, center = mu, cov = S) 
mahal

  #把「每個人」都依次當作中心，算其他人和他的距離。(Treat "each person" as the center in turn and calculate the distance between each person and him.)
  #用途：更接近 dist() 的概念，因為可以得到「樣本與樣本之間」的相似度。(Usage: Closer to the concept of dist(), because it can get the similarity between samples.)
mah_mat <- apply(vars, 1, function(x)
  mahalanobis(vars, center = x, cov = cov(vars)) #逐一把每一列（x）當作「中心點」，然後計算所有樣本到這個中心的馬氏距離。(# Treat each column (x) as the "center point" one by one, and then calculate the Mahalanobis distance of all samples to this center.)
)


# ==============================
# 5) Hierarchical Clustering
# ==============================
# Complete linkage
hc_complete <- hclust(d, method = "complete")
plot(hc_complete, labels = rownames(vars), main = "Cluster Dendrogram (Complete Linkage)")
rect.hclust(hc_complete, k = 3, border = 2)

# Single linkage
hc_single <- hclust(d, method="single")
plot(hc_single, labels=rownames(vars), main="Cluster Dendrogram (Single Linkage)")
rect.hclust(hc_single, k = 3, border = 2)

# Ward.D2
hc_ward <- hclust(d, method="ward.D2")
plot(hc_ward, labels=rownames(vars), main="Cluster Dendrogram (Ward.D2)")
rect.hclust(hc_ward, k = 3, border = 2)

# ==============================
# 6) K-Means Clustering
# ==============================
set.seed(123)
km3 <- kmeans(scale(vars), centers = 3, nstart = 20, iter.max = 1000)
km3

# 繪圖 (以 BMI vs SBP 為例)
plot(vars$bmxbmi, vars$bpxosy1,
     col = km3$cluster, pch = 19,
     xlab = "BMI", ylab = "SBP",
     main = "K-means Clustering (k=3)")

# 查看分群
km3$cluster

# ==============================
# 7) Choosing the number of clusters
# ==============================

## (A) Elbow Method
wss <- sapply(1:10, function(k){
  kmeans(scale(vars), centers = k, nstart = 20)$tot.withinss
})

plot(1:10, wss, type="b", pch=19, frame=FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares",
     main="Elbow Method")

## (B) Gap Statistic
set.seed(123)
gap <- clusGap(scale(vars), FUN = kmeans, nstart = 20, K.max = 12, B = 50)
plot(gap, main="Gap Statistic", ylim = c(-0.5,0.15))

## (C) Silhouette Statistic
fviz_nbclust(scale(vars), kmeans, method="silhouette")

