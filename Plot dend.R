coefData <- read.csv("C:/Users/luisr/Documentos/Documents/Documentos - Copia/NanoIA/Topicos_em_Nanotecnologia/Tabela de coeficientes NanoIA.csv")

head(coefData)

coefDataPca <- coefData[, c("b1", "b2", "b3", "b12", "b13", "b23")]

rownames(coefDataPca) <- coefData$Combination

distPair <- dist(
  scale(coefDataPca)
)

coefHc <- hclust(
  distPair,
  method = "ward.D2"
)

library(dendextend)

coefDend <- as.dendrogram(coefHc)

coefDend <- color_branches(
  coefDend,
  k = 4
)

png(
  "cluster.png",
  width = 18,
  height = 14,
  units = "cm",
  res = 600
)

par(
  mar = c(10, 4, 4, 2)
)

plot(
  coefDend,
  main = "Cluster Dendrogram",
  xlab = "",
  sub = "",
  cex = 1.2
)


dev.off()

#Now I'll plot the dendogram of the sole lipids

coefSing <- read.csv("coeficientesSingular.csv")

coefSing <- coefSing[rowSums(is.na(coefSing)) < ncol(coefSing), ]

coefSing <- coefSing[, colSums(is.na(coefSing)) < nrow(coefSing)]

coefSing <- coefSing[!is.na(coefSing$Lipid), ]

coefSing <- coefSing[coefSing$Lipid != "", ]

nrow(coefSing)

str(coefSing)

head(coefSing)

pcaData <- coefSing[, c(
  "b1",
  "b2",
  "b3",
  "b12",
  "b13",
  "b23"
)]

rownames(pcaData) <- coefSing$Lipid

distMat <- dist(scale(pcaData))

hc <- hclust(
  distMat,
  method = "ward.D2"
)

library(dendextend)

dend <- as.dendrogram(hc)

dend <- color_branches(
  dend,
  k = 4,
  col = c("#FA6F55", "#418DBA", "#CCC941", "#339E2C" )
)

dend <- set(dend, "branches_lwd", 2)

png(
  "cluster.png",
  width = 18,
  height = 14,
  units = "cm",
  res = 600
)

plot(
  dend,
  main = "Cluster Dendrogram of Lipids",
  xlab = "",
  sub = "",
  cex = 1.2
)

abline(
  h = 4,
  lty = 2,
  lwd = 2
)

dev.off()

#Now i'll try to save the heatmap

library(pheatmap)

pcaDataScale <- scale(pcaData)

heatMap <- pheatmap(
  pcaDataScale,
  clustering_method = "ward.D2",
  border_color = NA,
  fontsize_row = 12,
  fontsize_col = 12,
  main = "Lipid clustering by Scheffé's Coefficients"
)

class(heatMap)

heatMap

png(
  "pcaHeat.png",
  width = 18,
  height = 14,
  units = "cm",
  res = 600
)

heatMap

dev.off()
