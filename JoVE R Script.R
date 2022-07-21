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

plotQualityProfile(fnFs[3:6])
plotQualityProfile(fnRs[1:5])

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

plotErrors(errF, nominalQ=TRUE)

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

colnames(track) <- c("reads", "filtered", "dadaFs", "dadaRs", "merged", "nonchim")
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


#97% ASVs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa),
               sample_data(samdf))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#plot_richness(ps, x="SampleName", measures=c("Shannon"), color="Compost", shape = "Sample") + 
  #geom_point(size = 5)

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color = "CompostType", 
                shape = "SampleType", title="NMDS") + geom_point(size = 8) +
  geom_mark_ellipse(aes(color = CompostType, fill = CompostType))


top <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:800]
ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps <- prune_taxa(top, ps)

plot_bar(ps, x="SampleCode", fill="Family") + facet_wrap(~Compost, scales="free_x") +
  scale_fill_manual(values=colvec) + theme(legend.position="none")


#Unifraq
sequences<-getSequences(seqtab.nochim)[1:250]
names(sequences) <-sequences
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

plot_ordination(ps.prop, ord.nmds.bray, 
		color = "CompostType", shape = "SampleType", 
		title="Unifrac") + geom_point(size = 8) +
  geom_mark_ellipse(aes(color = CompostType, fill = CompostType))