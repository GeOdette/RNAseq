library(DT)
library(plotly)
library(gt)
install.packages("DT")
targets
group <- factor(targets$group)
#group <- factor(group)

distance <- dist(t(log2.cpm.filtered.norm), method = "maximum")
clusters <- hclust(distance, method = "average")
plot(clusters, labels = sampleLabels)
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale. = F, retx = T)
ls(pca.res)
summary(pca.res)
pca.res$rotationf
pca.res$x
screeplot(pca.res)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var/sum(pca.var)*100, 1)
pc.per
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  #geom_label() +
  stat_ellipse() +
  xlab(paste0("PC1 (", pc.per[1], "%",")")) +
  ylab(paste0("PC2 (", pc.per[2], "%",")")) +
  labs(title="PCA plot", 
        caption=paste0("produced on ", Sys.Date())) +
  coord_fixed() +
  theme_bw()
pca.res.df <- pca.res$x[,1:4] %>%
  as_tibble() %>%
  add_column(sample = sampleLabels, group = group)
pca.pivot <- pca.res.df |> pivot_longer(cols = PC1:PC4,
                                        names_to = "PC",
                                        values_to = "loadings")
pca.pivot |> ggplot() +
  aes(x=sample, y=loadings, fill = group) +
  geom_bar(stat = "identity") +
  facet_wrap(~PC) +
  labs(title = "PCA 'small multiples' plot",
       caption = paste0("produced on ", Sys.Date())) +
  theme_bw() +
  coord_flip()
