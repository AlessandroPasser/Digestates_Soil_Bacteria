#############################################################
#
# ARTICLE "Evolution of biogenic nitrogen from digestates for lettuce fertilization and the effect on the bacterial community"
# 
# Alessandro Passera
# alessandro.passera@unimi.it
# 
# Revision May 2024
# 
# script to reproduce calculations and figures presented in the manuscript
# 
#
#############################################################

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#1st time installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")


#required packages 
library("phyloseq")
library("ggplot2")


#retrieve R and package versions and compare to the uploaded file in gitHub for the reproducibility of the code
sessionInfo()

#set the correct working directory (insert correct folder path)
setwd("")


#############################################################
#import the count matrix and the design file
#############################################################


#OTU table 
dat_info <- read.delim("ont2illumina_OTUTable_10More.csv", sep = ",", header=T, row.names=1, blank.lines.skip = FALSE)

#inspect the file 
dim(dat_info)
colnames(dat_info)

#design files
designU <- read.delim("Mapping_File_DIMITRI_Unplanted.txt", sep = "\t", header=TRUE, row.names=1)
designP <- read.delim("Mapping_File_DIMITRI_Planted.txt", sep = "\t", header=TRUE, row.names=1)

#create files with count data
dat_count_U <-dat_info[, rownames(designU)] 
dat_count_P <-dat_info[, rownames(designP)] 

#inspect the files 
dim(dat_count_U)
colnames(dat_count_U)
dim(dat_count_P)
colnames(dat_count_P)

#check for chloroplast and mitochondria presence
Chloroplast <- dat_info[grepl("Chloroplast", dat_info$taxonomy), ]
dim(Chloroplast)
Chl_list <- c(rownames(Chloroplast))

mitochondria <- dat_info[grepl("mitochondria", dat_info$taxonomy), ]
dim(mitochondria)
Mit_list <- c(rownames(mitochondria))

#############################################################
#Genererate the phyloseq object
#############################################################

#The OTU Table counts
DIMITRI_OTU_U <- otu_table(dat_count_U, taxa_are_rows=TRUE)
DIMITRI_OTU_P <- otu_table(dat_count_P, taxa_are_rows=TRUE)

#The taxonomy information
DIMITRI_taxa_ordered <- read.delim ("ont2illumina_TaxTable_10More.csv", sep = ";", row.names=1, header=T, blank.lines.skip = FALSE)
DIMITRI_taxa <- tax_table(as.matrix(DIMITRI_taxa_ordered))
dim(DIMITRI_taxa)

#The mapping file 
DIMITRI_map_U <- sample_data(designU)
DIMITRI_map_P <- sample_data(designP)

#merge the files and create the phyloseq object
DIMITRI_phyloseq_U <- merge_phyloseq(DIMITRI_OTU_U, DIMITRI_taxa, DIMITRI_map_U)
DIMITRI_phyloseq_U
DIMITRI_phyloseq_P <- merge_phyloseq(DIMITRI_OTU_P, DIMITRI_taxa, DIMITRI_map_P)
DIMITRI_phyloseq_P

#Rarefy to even depth
sort(sample_sums(DIMITRI_phyloseq_U))
DIMITRI_rar_U <- rarefy_even_depth(DIMITRI_phyloseq_U, sample.size = min(sample_sums(DIMITRI_phyloseq_U)), rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

sort(sample_sums(DIMITRI_phyloseq_P))
DIMITRI_rar_P <- rarefy_even_depth(DIMITRI_phyloseq_P, sample.size = min(sample_sums(DIMITRI_phyloseq_P)), rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#############################
#Histograms
#############################

#Phylum level_Plant
DIMITRI_P_data_Phylum <- tax_glom(DIMITRI_rar_P, taxrank = "Phylum")
Phylum_color <- c("#1D91C0", "#67001F", "#F7FCFD", "#CB181D", "#78C679", "#F46D43", "#A6CEE3", "#FD8D3C", "#A6D854", "#6A51A3", "#7F0000", "#FFF7BC", "#000000", "#F0F0F0", "#C7EAE5", "#003C30", "#F16913", "#FFF7FB", "#8C6BB1", "#C7E9B4", "#762A83", "#FC9272", "#AE017E", "#F7F7F7", "#DF65B0", "#EF3B2C", "#74C476")
p <- plot_bar(DIMITRI_P_data_Phylum, fill = "Phylum", x = "Treatment", facet_grid = "Time")
p = p + scale_fill_manual(values = Phylum_color)
p

#Family level_Plant
DIMITRI_P_data_Family <- tax_glom(DIMITRI_rar_P, taxrank = "Family")
Top_N_Fams <- names(sort(taxa_sums(DIMITRI_P_data_Family), TRUE)[1:20])
DIMITRI_P_data_Family_Top <- prune_taxa(Top_N_Fams, DIMITRI_P_data_Family)
DIMITRI_P_data_Family_prune <- filter_taxa(DIMITRI_P_data_Family_Top, function(x) mean(x) >0.1, TRUE)
Fam_color <- c("#1D91C0", "#67001F", "#F7FCFD", "#CB181D", "#78C679", "#F46D43", "#A6CEE3", "#FD8D3C", "#A6D854", "#6A51A3", "#7F0000", "#FFF7BC", "#000000", "#F0F0F0", "#C7EAE5", "#003C30", "#F16913", "#FFF7FB", "#8C6BB1", "#C7E9B4", "#762A83", "#FC9272", "#AE017E", "#F7F7F7", "#DF65B0", "#EF3B2C", "#74C476")
p <- plot_bar(DIMITRI_P_data_Family_prune, fill = "Family", x = "Treatment", facet_grid = "Time")
p = p + scale_fill_manual(values = Fam_color)
p

#Phylum level_Unplanted
DIMITRI_U_data_Phylum <- tax_glom(DIMITRI_rar_U, taxrank = "Phylum")
Phylum_color <- c("#1D91C0", "#67001F", "#F7FCFD", "#CB181D", "#78C679", "#F46D43", "#A6CEE3", "#FD8D3C", "#A6D854", "#6A51A3", "#C7EAE5", "#7F0000", "#FFF7BC", "#000000", "#F0F0F0", "#003C30", "#F16913", "#FFF7FB", "#8C6BB1", "#C7E9B4", "#762A83", "#FC9272", "#AE017E", "#F7F7F7", "#DF65B0", "#EF3B2C", "#74C476")
p <- plot_bar(DIMITRI_U_data_Phylum, fill = "Phylum", x = "Treatment", facet_grid = "Time")
p = p + scale_fill_manual(values = Phylum_color)
p

#Family level_Unplanted
DIMITRI_U_data_Family <- tax_glom(DIMITRI_rar_U, taxrank = "Family")
Top_N_Fams <- names(sort(taxa_sums(DIMITRI_U_data_Family), TRUE)[1:20])
DIMITRI_U_data_Family_Top <- prune_taxa(Top_N_Fams, DIMITRI_U_data_Family)
DIMITRI_U_data_Family_prune <- filter_taxa(DIMITRI_U_data_Family_Top, function(x) mean(x) >0.1, TRUE)
Fam_color <- c("#1D91C0", "#67001F", "#F7FCFD", "#CB181D", "#78C679", "#F46D43", "#A6CEE3", "#FD8D3C", "#A6D854", "#6A51A3", "#7F0000", "#FFF7BC", "#000000", "#F0F0F0", "#C7EAE5", "#003C30", "#F16913", "#FFF7FB", "#8C6BB1", "#C7E9B4", "#762A83", "#FC9272", "#AE017E", "#F7F7F7", "#DF65B0", "#EF3B2C", "#74C476")
p <- plot_bar(DIMITRI_U_data_Family_prune, fill = "Family", x = "Treatment", facet_grid = "Time")
p = p + scale_fill_manual(values = Fam_color)
p