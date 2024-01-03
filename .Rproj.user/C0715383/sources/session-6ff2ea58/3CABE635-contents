#Load packages
packages_to_load <- c("tidyverse", "edgeR", "matrixStats", "cowplot")

# Assuming packages_to_load is a vector of package names
for (package in packages_to_load) {
  if (!(package %in% installed.packages()[, "Package"])) {
    tryCatch(
      {
        BiocManager::install(package)
      },
      error = function(e) {
        cat("Error occurred while installing", package, ":", conditionMessage(e), "\n")
      },
      finally = {
      }
    )
  }
  
  # Load the package
  library(package, character.only = TRUE)
}
myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts
myTranscriptTPM <- Txi_transcript$abundance
colSums(myTPM)
colSums(myTranscriptTPM)
colSums(myCounts)

sampleLabels <- targets$sample
myTPM.stats <- transform(myTPM,
                         SD=rowSds(myTPM),
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM)
  
)
head(myTPM.stats)
# ggplot ----
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) + 
  geom_point(shape=16, size=2) +
  geom_hex(show.legend = FALSE) +
  geom_smooth(method = lm) +
  labs(
    y="Median", x="Standard deviation",
    title="Transcripts per million (TPM)",
    subtitle = "unfiltered, non-normalized data",
    caption = "Ge Odette 2023"
  )
  theme_classic() +
  theme_dark() +
  theme_bw()

  #DGE LIST ----
myDGEList <- DGEList(myCounts)
save(myDGEList, file = "myDGEList")
load(file = myDGEList)
cpm <- cpm(myDGEList)
colSums(cpm)
log2.cpm <- cpm(myDGEList, log = TRUE)
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df,
                                  cols = HS01:CL13,
                                  names_to = "samples",
                                  values_to = "expression")
p1 <- ggplot(log2.cpm.df.pivot) + aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = 'median', 
               geom = 'point',
               shape = 95, 
               size = 10,
               color = 'black') +
  labs(y="log2 expression", x= "sample",
       title = "log2 Counts per Million (cpm)", 
       subtitle = "unfiltered, non-normalized",
       caption = paste0("prodiced on ", Sys.Date())) +
  theme_bw()

table(rowSums(myDGEList$counts==0)==10)
keepers <- rowSums(cpm>1)>=5
myDGEList.filtered <- myDGEList[keepers,]
dim(myDGEList.filtered)
log2.cpm.filtered <- cpm(myDGEList.filtered, log = TRUE)
log2.cpm.filtered.df = as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df,
                                  cols = HS01:CL13,
                                  names_to = "samples",
                                  values_to = "expression")
p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = 'median', 
               geom = 'point',
               shape = 95, 
               size = 10,
               color = 'black') +
  labs(y="log2 expression", x= "sample",
       title = "log2 Counts per Million (cpm)", 
       subtitle = "filtered, non-normalized",
       caption = paste0("prodiced on ", Sys.Date())) +
  theme_bw() 
  
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log =TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df,
                                                cols = HS01:CL13,
                                                names_to = "samples",
                                                values_to = "expression")
p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = 'median', 
               geom = 'point',
               shape = 95, 
               size = 10,
               color = 'black') +
  labs(y="log2 expression", x= "sample",
       title = "log2 Counts per Million (cpm)", 
       subtitle = "filtered, TMM normalized",
       caption = paste0("prodiced on ", Sys.Date())) +
  theme_bw() 
plot_grid(p1, p2, p3, labels = c("A", "B", "C"),label_size = 12)
p2

