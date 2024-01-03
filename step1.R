#A vector of needed packages ----
packages_to_install <- c("rhdf5", "tximport", "ensembldb", "EnsDb.Hsapiens.v86",
                         "biomaRt", "beepr")
install.packages("datapasta")
# Install needed packages ----
for (package in packages_to_install){
  BiocManager::install(package)
}

#Load needed packages ----
for (package in packages_to_install){
  library(package, character.only = TRUE)
}
library(datapasta)
library(tidyverse)


#Reading the study design ----
targets <- read_tsv("studydesign.txt")

path <- file.path(targets$sra_accession, "abundance.tsv")
all(file.exists(path))
#Get annotations ----
Tx <- transcripts(EnsDb.Hsapiens.v86, columns = c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)

#Rename tx_id to target-id ----
Tx <- dplyr::rename(Tx, target_id = tx_id) 
#Select target_id and gene_name columns ----
Tx <- dplyr::select(Tx, "target_id", "gene_name")

#Using biomat ----
listMarts()

myMart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
available.datasets <- listDatasets(myMart)
dog.anno <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "clfamiliaris_gene_ensembl")
dog.filters <- listFilters(dog.anno)
Tx.dog <- getBM(attributes = c("ensembl_transcript_id_version", "external_gene_name"), mart = dog.anno)
Tx.dog <- as_tibble(Tx.dog)

Tx.dog <- dplyr::rename(Tx.dog, target_id = ensembl_transcript_id_version, gene_name = external_gene_name )

#Read kallisto data ----
Txi_gene <- tximport(path, 
                     type = "kallisto",
                     tx2gene = Tx,
                     txOut = FALSE,
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE
                     )
names(Txi_gene)
