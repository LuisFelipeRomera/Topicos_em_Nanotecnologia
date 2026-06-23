setwd("C:/Users/luisr/Documentos/Documents/Documentos - Copia/NanoIA/planilhas")

bacuriAcai <- read.csv("C:/Users/luisr/Documentos/Documents/Documentos - Copia/NanoIA/planilhas/NanoBacuriAcai.csv")

head(bacuriAcai)

library(mixexp)

bacuriAcaiVar <- c("Bacuri", "Acai", "Brijo10")

bacuriAcaiMod1 <- MixModel(bacuriAcai, "PDI", bacuriAcaiVar, 1)

coef(bacuriAcaiMod1)

pcaData <- read.csv("C:/Users/luisr/Documentos/Documents/Documentos - Copia/NanoIA/Topicos_em_Nanotecnologia/Tabela de coeficientes NanoIA.csv")

head(pcaData)
str(pcaData)

pcaCol <- pcaData[, c("b1", "b2", "b3", "b12", "b13", "b23")]

pca <- prcomp(
  pcaCol,
  center = TRUE,
  scale. = TRUE
)

pca

summary(pca)

round(pca$rotation, 3)

round(pca$x, 3)

rownames(pcaCol) <- pcaData$Combination

pca <- prcomp(
  pcaCol,
  center = TRUE,
  scale. = TRUE
)

install.packages("factoextra")
library(factoextra)

fviz_pca_biplot(
  pca,
  repel = TRUE
)

pcaEig <- fviz_eig(pca)

pcaPlot <- fviz_pca_biplot(
  pca,
  repel = TRUE,
  pointsize = 4,
  labelsize = 4,
  arrowsize = 0.8
) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(size = 12)
  ) +
  labs(
    title = "PCA of Scheffé model's coefficients"
  )

pcaPlot

fviz_pca_ind(
  pca,
  habillage = pcaData$Model,
  addEllipses = TRUE,
  ellipse.type = "confidence",
  ellipse.level = 0.95,
  repel = TRUE
)

pcaColor <- fviz_pca_var(
  pca,
  col.var = "contrib",
  repel = TRUE,
  gradient.cols = c("blue", "orange", "red")
) +
  theme_classic(base_size = 14)

pcaElipse <- fviz_pca_ind(
  pca,
  geom = "point",
  habillage = pcaData$Model,
  addEllipses = TRUE,
  ellipse.level = 0.95,
  pointsize = 4,
  repel = TRUE
) +
  theme_classic(base_size = 14)

getwd()

ggsave(
  "pcaPlot.png",
  plot = pcaPlot,
  width = 18,
  height = 14,
  units = "cm",
  dpi = 600
)

setwd("C:/Users/luisr/Documentos/Documents/Documentos - Copia/NanoIA/Topicos_em_Nanotecnologia")

ggsave(
  "pcaEig.png",
  plot = pcaEig,
  width = 18,
  height = 14,
  units = "cm",
  dpi = 600
)

ggsave(
  "pcaColor.png",
  plot = pcaColor,
  width = 18,
  height = 14,
  units = "cm",
  dpi = 600
)

ggsave(
  "pcaElipse.png",
  plot = pcaElipse,
  width = 18,
  height = 14,
  units = "cm",
  dpi = 600
)

head(pcaCol)

pca

install.packages("pheatmap")
library(pheatmap)

pcaHeat <- pheatmap(
  pcaCol,
  clustering_method = "ward.D2"
)

ggsave(
  "pcaHeat.png",
  plot = pcaHeat,
  width = 18,
  height = 14,
  units = "cm",
  dpi = 600
)

distMat <- dist(pcaCol)

distMat

pcaHc <- hclust(
  distMat,
  method = "ward.D2"
)

plot(pcaHc)

png(
  "cluster.png",
  width = 18,
  height = 14,
  units = "cm",
  res = 600
)

plot(
  pcaHc,
  main = "Cluster Dendrogram",
  xlab = "",
  sub = "",
  cex = 1.2
)

dev.off()

install.packages("ggdendro")

library(ggdendro)
library(ggplot2)

dend <- dendro_data(pcaHc)

dendPlot <- ggplot(
  segment(dend)
) +
  geom_segment(
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    )
  ) +
  theme_classic()

dendPlot

clusterFviz <- fviz_dend(
  pcaHc,
  k = 3,
  rect = TRUE,
  cex = 0.9,
  lwd = 0.1,
  rect_lty = 1
)

ggsave(
  "clusterFviz.png",
  plot = clusterFviz,
  width = 18,
  height = 14,
  units = "cm",
  dpi = 300
)

fviz_dend(
  pcaHc,
  k = 3,
  rect = TRUE,
  rect_fill = TRUE,
  rect_border = "black",
  cex = 0.9,
  lwd = 1
) +
  theme_classic()

install.packages("dendextend")
library(dendextend)

dend <- as.dendrogram(pcaHc)

dend <- color_branches(
  dend,
  k = 3
)

plot(dend)

png(
  "cluster.png",
  width = 18,
  height = 14,
  units = "cm",
  res = 600
)

plot(
  dend,
  main = "Cluster Dendrogram",
  xlab = "",
  sub = "",
  cex = 1.2
)

dev.off()

t.test(
  b12 ~ Murumuru,
  data = pcaData
)

lm(
  b12 ~ Murumuru,
  data = pcaData
)

t.test(
  b13 ~ Pequi,
  data = pcaData
)
