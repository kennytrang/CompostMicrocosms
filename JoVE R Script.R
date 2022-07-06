library("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(ggforce)
theme_set(theme_bw())
library(RColorBrewer)
library(dplyr)
library(vegan)
library(indicspecies)
library(DECIPHER)
library(phangorn)
library(decontam)
library(microshades)
library(hrbrthemes)
library(scales)
library(picante)


path <- "C:\\Users\\kenne\\Documents\\Academics\\Research\\Shapira\\Data\\Illumina\\5.24.22"
list.files(path)


fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#plotQualityProfile(fnFs[3:6])
#plotQualityProfile(fnRs[1:5])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)


errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = 6)

seqtab <- makeSequenceTable(mergers)
#dim(seqtab)
#table(nchar(getSequences(seqtab)))


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, 
                       "C:\\Users\\kenne\\Documents\\Academics\\Research\\Shapira\\Data\\Illumina\\silva_nr_v132_train_set.fa.gz", 
                       multithread=TRUE)
#taxa.sp <- addSpecies(taxa, "C:\\Users\\kenne\\Documents\\Academics\\Research\\Shapira\\Data\\Illumina\\silva_species_assignment_v132.fa.gz")


#phyloseq
row1 <- rownames(seqtab.nochim)
row2 <- c("BP Worm #1", "BP Worm #2", "BP Worm #3", 
          "AP Worm #1", "AP Worm #2", "AP Worm #3", 
          "POT Worm #1", "POT Worm #2", "POT Worm #3", 
          "BAN Worm #1", "BAN Worm #2", "BAN Worm #3", 
          "OR Worm #1", "OR Worm #2", "OR Worm #3",
          "C", "C",
          "BP Soil #1", "BP Soil #2", "BP Soil #3", 
          "AP Soil #1", "AP Soil #2", "AP Soil #3", 
          "POT Soil #1", "POT Soil #2", "POT Soil #3",
          "BAN Soil #1", "BAN Soil #2", "BAN Soil #3", 
          "OR Soil #1", "OR Soil #2", "OR Soil #3")
row3 <- c("Worm #1", "Worm #2", "Worm #3", 
          "Worm #1", "Worm #2", "Worm #3", 
          "Worm #1", "Worm #2", "Worm #3", 
          "Worm #1", "Worm #2", "Worm #3",  
          "Worm #1", "Worm #2", "Worm #3", 
          "C", "C",
          "Soil #1", "Soil #2", "Soil #3", 
          "Soil #1", "Soil #2", "Soil #3",  
          "Soil #1", "Soil #2", "Soil #3", 
          "Soil #1", "Soil #2", "Soil #3",  
          "Soil #1", "Soil #2", "Soil #3")
row4 <- c("Worm", "Worm", "Worm", "Worm", "Worm", "Worm", "Worm", "Worm", "Worm", 
          "Worm", "Worm", "Worm", "Worm", "Worm", "Worm", 
          "C", "C",
          "Soil", "Soil", "Soil","Soil", "Soil", "Soil","Soil", "Soil", "Soil",
          "Soil", "Soil", "Soil","Soil", "Soil", "Soil")
row5 <- c("Bell Pepper", "Bell Pepper", "Bell Pepper", 
          "Apple", "Apple", "Apple", 
          "Potato", "Potato", "Potato", 
          "Banana", "Banana", "Banana", 
          "Orange", "Orange", "Orange", 
          "C", "C",
          "Bell Pepper", "Bell Pepper", "Bell Pepper", 
          "Apple", "Apple", "Apple", 
          "Potato", "Potato", "Potato", 
          "Banana", "Banana", "Banana", 
          "Orange", "Orange", "Orange")

row6 <- c("True Sample", "True Sample", "True Sample", 
          "True Sample", "True Sample", "True Sample", 
          "True Sample", "True Sample", "True Sample", 
          "True Sample", "True Sample", "True Sample",  
          "True Sample", "True Sample", "True Sample", 
          "Control Sample", "Control Sample", 
          "True Sample", "True Sample", "True Sample", 
          "True Sample", "True Sample", "True Sample", 
          "True Sample", "True Sample", "True Sample", 
          "True Sample", "True Sample", "True Sample",  
          "True Sample", "True Sample", "True Sample")

row7 <- c("BP Worm", "BP Worm", "BP Worm", 
          "AP Worm", "AP Worm", "AP Worm", 
          "POT Worm", "POT Worm", "POT Worm", 
          "BAN Worm", "BAN Worm", "BAN Worm #3", 
          "OR Worm ", "OR Worm", "OR Worm #3",
          "C", "C",
          "BP Soil", "BP Soil", "BP Soil", 
          "AP Soil", "AP Soil", "AP Soil", 
          "POT Soil", "POT Soil", "POT Soil",
          "BAN Soil", "BAN Soil", "BAN Soil", 
          "OR Soil", "OR Soil", "OR Soil ")

samdf <- data.frame(Subject=row1, SampleName = row2, SampleCode = row3, Sample = row4, Compost = row5,
                    Sample_or_Control = row6, SampleType = row7)
rownames(samdf) <- row1

n <- 200
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colvec = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colvec <- colvec[c(1:7,9:200)][c(1:14,16:199)]


#97% ASVs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa),
               sample_data(samdf))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#decontam
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.3)
table(contamdf.prev$contaminant)
ps <- prune_taxa(!contamdf.prev$contaminant, ps)
ps <- subset_samples(ps, Sample_or_Control != "Control Sample")
#ps <- subset_samples(ps, Compost != "Banana")

#plot_richness(ps, x="SampleName", measures=c("Shannon"), color="Compost", shape = "Sample") + 
  #geom_point(size = 5)

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color = "CompostType", 
                shape = "SampleType", title="NMDS") + geom_point(size = 8) +
  geom_mark_ellipse(aes(color = CompostType, fill = CompostType))


top <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:500]
ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps <- prune_taxa(top, ps)

plot_bar(ps, x="SampleCode", fill="Family") + facet_wrap(~CompostType, scales="free_x")+
  scale_fill_manual(values=colvec)


#relabun
taxa.df <- as.data.frame(ps@tax_table)
seqtab.nochim.df <- t(as.data.frame(ps@otu_table))

TaxaCounts <- cbind(taxa.df, seqtab.nochim.df)
colnames(TaxaCounts)[7:36] <- row2[c(1:15,18:32)]
TaxaCounts <- TaxaCounts[!is.na(TaxaCounts$Kingdom),]

#test example for Firmicutes
total_reads <- numeric(30)
for(i in 7:36) {
  total_reads[i-6] = sum(TaxaCounts[,i])
}

Firmicutes <- TaxaCounts[!is.na(TaxaCounts[,2]),][TaxaCounts[!is.na(TaxaCounts[,2]),][,2] == "Firmicutes",]

test <- numeric(30)
for(i in 1:30) {
  test[i] = sum(Firmicutes[,i+6])/total_reads[i]
}

#Function outputs vector for relative abundances
relabun <- function(Level, Taxon){
  #Translate Level to Number
  for(i in 1:7){
    taxonomic <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    if(Level == taxonomic[i]){Level <- i} 
    Level
  }
  #Define Total Reads
  total_reads <- numeric(30)
  for(i in 7:36) {
    total_reads[i-6] = sum(TaxaCounts[,i])
  }
  
  
  Tax <- TaxaCounts[!is.na(TaxaCounts[,Level]),][TaxaCounts[!is.na(TaxaCounts[,Level]),][,Level] == Taxon,]
  
  output <- numeric(30)
  for(i in 1:30) {
    output[i] = sum(Tax[,i+6])/total_reads[i]
  }
  output
}

#enrichment

total <- numeric(30)
for (i in 1:30){
  total[i] <- sum(TaxaCounts[,i+6])
}
enrichment <- function(Family){
  results <- numeric(15)
  table <- na.omit(TaxaCounts[TaxaCounts$Family == Family, 7:36])
  for (i in 1:15){
    results[i] <- (sum(table[i,])/total[i]) / (sum(table[i+15,])/total[i+15])
  }
  return(results)
}
fam <- as.vector(levels(as.factor(TaxaCounts$Family)))
results <- matrix(0, nrow=30, ncol = length(fam))
for (i in seq_along(fam)){
  results[,i] <- enrichment(as.character(fam[i])) 
}
rownames(results) <- row2[c(1:15, 18:32)]
colnames(results) <- as.character(fam)

write.csv(results, 
          "C:\\Users\\kenne\\Documents\\Academics\\Research\\Shapira\\Data\\Microbiome Composition\\5.24.22\\enrich.csv")

table <- read.csv("C:\\Users\\kenne\\Documents\\Academics\\Research\\Shapira\\Data\\Microbiome Composition\\5.24.22\\Enrichment.csv")
#collapse max value to 3
for (i in 1:nrow(table)){
  if(table[i,3] > 3){
    table[i, 3] <- 3 
  }
}

table$y <- as.character(table$y)
table$y <- factor(table$y, levels=rev(as.character(levels(as.factor(table$y)))))
my_breaks <- c(0, 1, 2, 3, 4)
ggplot(table, aes(x, y, fill = z)) + 
  geom_tile() + scale_fill_gradient2(name = "log2(enrichment)", breaks = my_breaks, labels = my_breaks,
                                     low = "pink", high = "red")

#microshades
prep <- prep_mdf(ps, subgroup_level = "Family")
color <- create_color_dfs(prep, group_level = "Phylum", subgroup_level = "Family", cvd = TRUE,
                          selected_groups = c("Firmicutes", "Actinobacteria",
                                              "Bacteroidetes", "Planctomycetes", "Proteobacteria"))

mdf <- color$mdf
cdf <- color$cdf
cdf <- color_reassign(cdf,  group_assignment = c("Firmicutes", "Actinobacteria",
                                                 "Bacteroidetes", "Planctomycetes", "Proteobacteria"),
                      color_assignment = c("micro_cvd_purple", "micro_cvd_blue", 
                                           "micro_purple", "micro_cvd_turquoise", "micro_cvd_orange"))

plot <- plot_microshades(mdf, cdf, group_label = "Phylum-Family", x= "SampleCode")
plot + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(1.5, "cm"), text=element_text(size=10)) +
  theme(axis.text.x = element_text(size= 6)) + facet_wrap(~CompostType, scales="free_x")




#Proteobacteria Only
#ps <- subset_samples(ps, CompostType != "Banana")
ps_prot <- subset_taxa(ps, Phylum == "Proteobacteria")

prep <- prep_mdf(ps_prot, subgroup_level = "Family")
color <- create_color_dfs(prep, group_level = "Class", subgroup_level = "Family", cvd = TRUE,
                          selected_groups = c("Alphaproteobacteria", "Gammaproteobacteria", 
                                              "Deltaproteobacteria"))
mdf <- color$mdf
cdf <- color$cdf
cdf <- color_reassign(cdf, group_assignment = c("Alphaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria"),
                      color_assignment = c("micro_cvd_turquoise", "micro_cvd_orange","micro_cvd_blue"), group_level = "Class")

plot_1 <- plot_microshades(mdf, cdf, group_label = "Class-Family", x= "SampleCode")
plot_1 + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(1.5, "cm"), text=element_text(size=10)) +
  theme(axis.text.x = element_text(size= 6)) + facet_wrap(~CompostType, scales="free_x")


#Unifraq
sequences<-getSequences(seqtab.nochim)[1:250]
names(sequences)<-sequences
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

ps <- phyloseq(otu_table(seqtab.nochim[,1:250], taxa_are_rows=FALSE), 
                  tax_table(taxa), phy_tree(fitGTR$tree), sample_data(samdf))

phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.3)
table(contamdf.prev$contaminant)
ps <- prune_taxa(!contamdf.prev$contaminant, ps)
ps <- subset_samples(ps, Sample_or_Control != "Control Sample")
#ps <- subset_samples(ps, CompostType != "Banana")

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method = "PCoA", distance = "wunifrac")

plot_ordination(ps.prop, ord.nmds.bray, color = "CompostType", 
                shape = "SampleType", title="Unifrac") + geom_point(size = 8) +
  geom_mark_ellipse(aes(color = CompostType, fill = CompostType))


#Faith's Phylogenetic Diversity
ps <- phyloseq(otu_table(seqtab.nochim[,1:250], taxa_are_rows=FALSE), 
               tax_table(taxa), phy_tree(fitGTR$tree), sample_data(samdf))

phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.3)
table(contamdf.prev$contaminant)
ps <- prune_taxa(!contamdf.prev$contaminant, ps)
ps <- subset_samples(ps, Sample_or_Control != "Control Sample")
ps <- subset_samples(ps, Compost != "Banana")

psOTU <- as.data.frame(ps@otu_table)
pstree <- ps@phy_tree

phylo_div <- pd(psOTU, pstree,include.root=TRUE)
phylo_div <- cbind(phylo_div$PD, samdf[(samdf$Compost != "Banana" & samdf$Compost != "C"), ])
colnames(phylo_div)[1] <- "PD"
phylo_div <- cbind(phylo_div[,1:2], lapply(phylo_div[,3:8], factor))


ggplot(phylo_div, aes(x=SampleName, y=PD)) + geom_point(aes(color= Compost, shape = Sample), size = 5) + 
  theme(axis.text.x = element_text(angle = 90)) 
