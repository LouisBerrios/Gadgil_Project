# Gadgil 16S R Script# Louis Berrios

# In terminal data sorting
# Extract all files from Miseq subfolders and move them to new folders (made in terminal)
mkdir /Documents/Gadgil_16S
mkdir /Documents/Gadgil_ITS
find /Documents/Gadgil_Project/MS_2x301_LB-1-PRFB-Gadgil2_LouisBerrios022323-382493243/FASTQ_Generation_2023-02-26_08_35_28Z-656895240/ -type f -name '*16s*' -exec mv {} /Documents/Gadgil_16S/ \;
find /Documents/Gadgil_Project/MS_2x301_LB-1-PRFB-Gadgil2_LouisBerrios022323-382493243/FASTQ_Generation_2023-02-26_08_35_28Z-656895240/ -type f -name '*ITS*' -exec mv {} /Documents/Gadgil_ITS/ \;

# Now that all the files have been sorted, begin processing the fastq files using DADA2

# DADA2 workflow for processing 16S raw sequences
# Follows https://benjjneb.github.io/dada2/tutorial.html

library(dada2)
library(ShortRead)
library(Biostrings)


#############################

path <- "~/Documents/Gadgil_16S"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq.gz and SAMPLENAME_R2_001.fastq.gz
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq.gz", full.names = TRUE)) # If set to be FALSE, then working directory must contain the files
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq.gz", full.names = TRUE))

# Remove any forward files that don't have reverse counterparts, and vise versa
# (filterAndTrim will throw an error if fnFs and fnRs have any mismatches)
basefilenames_Fs <- sub("_L001_R1_001.fastq.gz","",basename(fnFs))
basefilenames_Rs <- sub("_L001_R2_001.fastq.gz","",basename(fnRs))
rm_from_fnFs <- basefilenames_Fs[which(!(basefilenames_Fs %in% basefilenames_Rs))]
rm_from_fnRs <- basefilenames_Rs[which(!(basefilenames_Rs %in% basefilenames_Fs))]

for(name in rm_from_fnFs) {
  print(paste(name, "does not have a reverse-reads counterpart. Omitting from this analysis."))
  fnFs <- fnFs[-which(fnFs == paste0(path, "/", name, "_L001_R1_001.fastq.gz"))]
}
for(name in rm_from_fnRs) {
  print(paste(name, "does not have a forward-reads counterpart. Omitting from this analysis."))
  fnRs <- fnRs[-which(fnRs == paste0(path, "/", name, "_L001_R2_001.fastq.gz"))]
}


# Identify primers - used 515F & 806R from Hiro's spreadsheet
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer seq
REV <- "GGACTACNVGGGTWTCTAAT"  ## CHANGE ME...


# Get all orientations of primers, just to be safe
# Note - changed this for the dimensions project due to different sequencing primers
# Used

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# “pre-filter” the sequences just to remove those with Ns, but perform no other filtering
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 1, multithread = TRUE)

# (From tutorial) We are now ready to count the number of times the primers appear in the 
# forward and reverse read, while considering all possible primer orientations. 
# Identifying and counting the primers on one set of paired end FASTQ files is
# sufficient, assuming all the files were created using the same library preparation,
# so we’ll just process the first sample.

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# If you see the reverse-complement of the forward primer in the reverse reads (cells [2,4] and [3,4]),
# it's because the ITS region is short and it is reading part of the forward primer.

# Remove primers using cutadapt

cutadapt <- "/Users/louisberrios/Documents/Cutadapt/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs.filtN)) {
  # for(i in 1:10) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i], # input files; fnFs.filtN replaced by fnFs.filtN, etc.
                             "--minimum-length", "1")) # min length of cutadapted reads: >0 
}

# Count primers in first post-cutadapt sample (should all be 0):
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[50]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[50]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[50]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[50]]))

# Since they are zero, skip step to remove other orientations of primers

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) {
  paste(strsplit(basename(fname), split="_")[[1]][1:3], collapse="_")
}
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Inspect read quality profiles of forward reads #1-2
plotQualityProfile(cutFs[1:12])

# Inspect read quality profiles of reverse reads #1-2
plotQualityProfile(cutRs[1:12])

# Filter and trim

# Assigning the filenames for the output of the filtered reads 
# to be stored as fastq.gz files.
filtFs <- file.path(path, "filtered", basename(fnFs.filtN))
filtRs <- file.path(path, "filtered", basename(fnRs.filtN))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

filtFs.out <- list.files(paste(path, "filtered", sep="/"), pattern="_L001_R1_001.fastq.gz", full.names=TRUE)
filtRs.out <- list.files(paste(path, "filtered", sep="/"), pattern="_L001_R2_001.fastq.gz", full.names=TRUE)

# Learn the error rates
errF <- learnErrors(filtFs.out, multithread = TRUE)
errR <- learnErrors(filtRs.out, multithread = TRUE)

# Visualize estimated error rates
plotErrors(errF, nominalQ = TRUE)

# Dereplicate identical reads
derepFs <- derepFastq(filtFs.out, verbose = TRUE)
derepRs <- derepFastq(filtRs.out, verbose = TRUE)
# Name the derep-class objects by the sample names
get.sample.name <- function(fname) {
  paste(strsplit(basename(fname), "_")[[1]][1:3], collapse="_")
}
sample.names <- unname(sapply(filtFs.out, get.sample.name))

# DADA2's core sample inference algorithm
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)


# Merge pairs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,trimOverhang = TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

saveRDS(seqtab,"Gadgil_16S_seqtab.bac.rds")
saveRDS(seqtab.nochim,"Gadgil_16S_seqtab.bac.nochim.rds")

# Inspect distribution of sequence lengths
hist(nchar(getSequences(seqtab.nochim)))

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))

#format out to accommodate dropped samples
raw.sample.names <- unname(sapply(row.names(out), get.sample.name))

track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))

track2<-cbind(out,track[match(row.names(out),row.names(track)),])

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                      "nonchim")
write.csv(track2,"Gadgil_16S_summary.csv")

rownames(track) <- sample.names
head(track2)


# Assign taxonomy using the SILVA database
silva.ref<-"silva_nr99_v138.1_train_set.fa.gz"
silva.species<-"silva_species_assignment_v138.1.fa.gz"

taxa <- assignTaxonomy(seqtab.nochim, silva.ref, multithread = TRUE, tryRC = TRUE)

taxa <- addSpecies(taxa, silva.species)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Save OTU table and taxonomic table as RDS files
# to hand off to dada2_to_phyloseq.R
saveRDS(seqtab.nochim, "Gadgil_16S_seqtab_nochim_taxa.rds")
saveRDS(taxa, "Gadgil_16S_taxa.rds")

########Construct phyloseq object with sample data variables#########

#read in environmental sample table

sample <- readRDS("Gadgil_16S_seqtab_nochim_taxa.rds")
SD <- read.csv("Gadgil_T1-T4_SampData16S.csv", row.names = 1)
SAM <- sample_data(SD, errorIfNULL = T)
taxa <- readRDS("Gadgil_16S_taxa.rds")

#create a phyloseq object
bac.ps <- phyloseq(otu_table(sample, taxa_are_rows=FALSE), sample_data(SAM), tax_table(taxa))

#filter out unwanted taxa (e.g., mitochondira and chloroplast sequences)
bac.ps.filt<-subset_taxa(bac.ps,Family!="Mitochondria")
bac.ps.filt<-subset_taxa(bac.ps.filt,Genus!="Chloroplast")

# Removing sequence rownames for display only
taxa.print <- tax_table(bac.ps.filt)
rownames(taxa.print) <- NULL
head(taxa.print)

#save the filtered dataset 
saveRDS(bac.ps.filt,"Gadgil.16S.filtered.rds")

#filter out low abundant sequences
bac.ps.filt2 = prune_taxa(taxa_sums(bac.ps.filt) > 10, bac.ps.filt) 
bac.ps.filt2 = prune_samples(sample_sums(bac.ps.filt2)>1000, bac.ps.filt2)

#save the filtered+pruned dataset
saveRDS(bac.ps.filt2, "Gadgil.16S.filt-prune.rds")
bac.filt <- readRDS("Gadgil.16S.filt-prune.rds")
#rarefy the dataset
bac.ps.rare <- rarefy_even_depth(bac.filt, sample.size = min(sample_sums(bac.filt)),
                                 rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
# Remove Archaea from the dataset
bac.ps.rare2 <- subset_taxa(bac.ps.rare, Kingdom != "Archaea")

#save rarefied dataset
saveRDS(bac.ps.rare2, "Rarefied-Gadgil-Bacteria.rds")
#create ASV-TAX table
taxa_names(bac.ps.rare2) <- paste0("Seq", seq(ntaxa(bac.ps.rare2)))
ASV_BAC <- data.frame(otu_table(bac.ps.rare2))
TAX_ASV <- data.frame(tax_table(bac.ps.rare2))
ASV_BAC_T <- t(ASV_BAC)
MERGER2 <- merge(ASV_BAC_T, TAX_ASV, by = "row.names")
write.csv(MERGER2, "ASV_Bacteria_Gadgil_Rare_Final.csv")

#create the combined ASV table plus taxonomy string metacoder needs

#read it in using phyloseq object
tax_table(bac.ps.rare)<-tax_table(bac.ps.filt2)[,1:7]
x1<-parse_phyloseq(bac.ps.rare2)

#transform data [compositional]

pseq_final <- microbiome::transform(bac.ps.rare2, "compositional")

#save as RDS
saveRDS(pseq_final, "Gadgil-Bacteria-PS_filter+rare+tran.rds")
pseq <- readRDS("Gadgil-Bacteria-PS_filter+rare+tran.rds")
#create dataframe 
pseqDF <- psmelt(pseq)
pseqdf1 <- subset(pseqDF, Condition !="NA") # use this data frame for downstream analyses
library(RColorBrewer)
#get appropriate number of colors for data
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(pseqDF$Phylum))

#plot the relative abundance by Phylum
ggplot(pseqdf1, aes(fill=Phylum, y=Abundance, x=Condition)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + 
  scale_fill_manual(values = c(
  "#5f4b42",
                               "#b8575f",
                               "#d38168",
                               "#d5bf73",
                               "#30949e",
                               "#4e6189",
                               "#8c98b0",
                               "#dedede",
                               "violetred4",
                               "#96b9a9",
                               "#e6dbb2",
                               "gray",
                               "#e6c57f",
                               "steelblue",
                               "#565655",
                               "#2e5e60",
                               "#5b8656")) + 
  ylab("Relative Abundance") + theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_y_continuous(labels = scales::percent) + facet_grid(~Timepoint) + xlab("") 

# save plot 
ggsave(
  filename = "Gadgil_16S_Compositional_Stacked_BarPlot.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 7,
  height = 5.5,
  units = c("in"),
  dpi = 300)

##### GENUS level investigations using Curren Bio indicator bacteria as guides
# Conexibacter, Acidocella, Mycobacterium, Sphingomonas, Bradyrhizobium, Reyranella, Puia, Acidobacter, Burkholderia
# Bryobacter, Candidatus Solibacter

indicators <- subset(pseqdf1, Genus == "Conexibacter" | Genus == "Acidocella" | Genus == "Mycobacterium" | 
                       Genus == "Sphingomonas" | Genus == "Bradyrhizobium" | Genus == "Reyranella" |
                       Genus == "Puia" | Genus == "Acidibacter" | Genus == "Burkholderia-Caballeronia-Paraburkholderia" |
                       Genus == "Bryobacter" | Genus == "Candidatus Solibacter") 
#plot the relative abundance by Phylum
ggplot(indicators, aes(fill=Genus, y=Abundance, x=Condition)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + 
  scale_fill_manual(values = c(
    "#5f4b42",
    "#b8575f",
    "#d38168",
    "#d5bf73",
    "#30949e",
    "#4e6189",
    "#8c98b0",
    "#dedede",
    "violetred4",
    "#96b9a9",
    "#e6dbb2",
    "gray",
    "#e6c57f",
    "steelblue",
    "#565655",
    "#2e5e60",
    "#5b8656")) + 
  ylab("Relative Abundance") + theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_y_continuous(labels = scales::percent) + facet_grid(~Timepoint) + xlab("") 

ggsave(
  filename = "Gadgil_16S-Indicators_Compositional_Stacked_BarPlot.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 7,
  height = 5.5,
  units = c("in"),
  dpi = 300)

####### Alpha Diversity ###########
# rarefied data
rare <- readRDS("Rarefied-Gadgil-Bacteria.rds")
rare2 <- subset_samples(rare, Condition !="NA")
ad <- estimate_richness(rare2)
write.csv(ad, "Bacteria_Alpha_Diversity.FINAL.csv")
library(ggpubr)
rare_alpha_div <- read.csv("Bacteria_Alpha_Diversity.FINAL.csv") 

rare_alpha_div$title <- "Bacterial Richness"
# Models
aov_ad.div <- aov(Observed ~ Condition:Timepoint + Site, rare_alpha_div) # rarefied
aov_ad.div2 <- aov(Observed ~ Condition + Timepoint + Condition:Timepoint, rare_alpha_div)
summary(aov_ad.div)
summary(aov_ad.div2)
# Add stats to alpha diversity plot (Observed Richness)
ggplot(rare_alpha_div, aes(x=Condition, y=Observed, color = Condition)) + geom_boxplot() + theme_bw() + 
  scale_color_manual(values=c("firebrick3", "gray1")) + xlab("") + 
  theme(axis.text.x = element_text(size = 12)) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size = 12)) + 
  facet_grid(~title) + theme(strip.text = element_text(size = 14, face = "bold")) + 
  theme(axis.title.y = element_text(size = 16, vjust = 1.5)) + 
  theme(strip.background = element_rect(fill = "lightblue4")) + 
  annotate(geom = "text", x = 0.94, y = 305, size = 5, label = paste("Condition:F[\"1,154\"] ==", 3.892), parse=TRUE) +
  annotate(geom = "text", x=1.61, y=305, size = 5, label = ", p = 0.0503") +
  annotate(geom = "text", x = 0.94, y = 290, size = 5, label = paste("Timepoint:F[\"3,154\"] ==", 3.436), parse=TRUE) +
  annotate(geom = "text", x=1.6, y=290, size = 5, label = ", p = 0.0185")
  
# Save plot
ggsave(
  filename = "Gadgil_16S_ObservedRichnessBoxPlot.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 5,
  units = c("in"),
  dpi = 300)

# plot graph showing changes across time
ggplot(rare_alpha_div, aes(x=Timepoint, y=Observed, color = Condition)) + geom_boxplot() + theme_bw() + 
  scale_color_manual(values=c("firebrick3", "gray1")) + xlab("") + 
  theme(axis.text.x = element_text(size = 12)) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size = 12)) + 
  facet_grid(~title) + theme(strip.text = element_text(size = 14, face = "bold")) + 
  theme(axis.title.y = element_text(size = 16, vjust = 1.5)) + 
  theme(strip.background = element_rect(fill = "lightblue4"))
# Save plot
ggsave(
  filename = "ObservedRichnessBoxPlot_TIMEPOINT.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4,
  height = 5,
  units = c("in"),
  dpi = 300)

####### NMDS ordination ###########
# remove mock
pseq2 <- subset_samples(pseq, Condition !="NA")
#perform an NMDS ordination
ord.bac2 <- ordinate(rare2, "NMDS", "bray", k=5)
#plot NMDS ordination [axes 2:3]
#add facet label to phyloseq object
sample_data(rare2)$BETA <- "Bacterial Beta Diversity"
plot_ordination(rare2, ord.bac2, type="sample", color="Timepoint", shape="Condition", axes=2:3) + 
  theme_bw() + scale_color_manual(values = c("#9a031e", "#0f4c5c", "#F8AFA8", "#f48c06")) + 
  geom_point(size = 4, alpha = 1) + stat_ellipse(type ="norm", linetype = 2) +
  annotate("text", x = -0.95, y = 2, label = "Stress = 0.0838") +
  annotate("text", x = -0.85, y = 1.8, label = "Condition: p = 0.019") +
  annotate("text", x = -0.84, y = 1.6, label = "Timepoint: p = 0.001") + facet_grid(~BETA) +
  theme(strip.text = element_text(size = 12, face = "bold")) + 
  theme(strip.background = element_rect(fill = "#c18c5d"))
# Save
ggsave(
  filename = "Gadgil_16S_BETADIVERSITY.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 5,
  units = c("in"),
  dpi = 300)

########PERMANOVA###########

#calculate bray-curtis distance matrix
library(vegan)
rareT1 <- subset_samples(rare2, Timepoint == "T1")
rareT2 <- subset_samples(rare2, Timepoint == "T2")
ord.bac.bray1 <- phyloseq::distance(rareT1, method = "bray", k=2)
ord.bac.bray <- phyloseq::distance(pseq2, method = "bray", k=4)
ord.bac.bray2 <- phyloseq::distance(rare2, method = "bray", k=4) # Keep this one. This uses rarefied data AND removes NA condition.
ord.bac.bray3 <- phyloseq::distance(rareT2, method = "bray", k=2)

#make a dataframe from the sample_data
ord.bac.bray.sampleDF <- data.frame(sample_data(pseq2))
ord.bac.bray.sampleDF2 <- data.frame(sample_data(rare2))
ord.bac.bray.sampleDF3 <- data.frame(sample_data(rareT1))
ord.bac.bray.sampleDF4 <- data.frame(sample_data(rareT2))

#adonis test
PERMANOVA <- adonis2(ord.bac.bray ~ Condition, data = ord.bac.bray.sampleDF)
PERMANOVA2 <- adonis2(ord.bac.bray ~ Condition / Timepoint, data = ord.bac.bray.sampleDF, strata = ord.bac.bray.sampleDF$Timepoint)
PERMANOVA3 <- adonis2(ord.bac.bray ~ Condition + Timepoint, data = ord.bac.bray.sampleDF)
PERMANOVA4 <- adonis2(ord.bac.bray2 ~ Condition + Timepoint + Timepoint:Condition, data = ord.bac.bray.sampleDF2) # Keep this one.

#write PERMANOVA results to table CSV
PERMANOVA4.df <- data.frame(PERMANOVA4)
write.csv(PERMANOVA4.df, "PHC_Bac_Permanova.csv")
#run betadisper function on distance matrix
beta_adonis <- betadisper(ord.bac.bray2, ord.bac.bray.sampleDF2$Condition, bias.adjust = TRUE)
beta_adonis2 <- betadisper(ord.bac.bray2, ord.bac.bray.sampleDF2$Timepoint, bias.adjust = TRUE)
#stats
stat_disp_anova <- anova(beta_adonis)
stat_disp_anova2 <- anova(beta_adonis2)
plot(beta_adonis, hull=FALSE, ellipse=TRUE)
plot(beta_adonis2, hull=FALSE, ellipse=TRUE)
#pairwise permutation test for homogeneity of multivariate dispersions
permutest(beta_adonis2, pairwise = TRUE, permutations = 1000)
permutest(beta_adonis, pairwise = TRUE, permutations = 1000)

########### Genera-specific changes across treatment and time ##############
# need to do more with these

burk <- subset(pseqdf1, Genus == "Burkholderia-Caballeronia-Paraburkholderia")
burk$title <- "Burkholderia ASVs"
burk2 <- subset(burk, Abundance > 0) # remove zeros from data
burk2.agg <- aggregate(Abundance ~ Sample + Condition + Timepoint, FUN=sum, data = burk2)
burk2.agg$title <- "Burkholderia ASVs"
library(ggpubr)
my_comps <- list(c("Trenched", "Untrenched"))

library(ggplot2)
library(ggpubr)
ggplot(burk2.agg, aes(x=Condition, y=log(Abundance), color=Condition)) + 
  scale_color_manual(values=c("red3", "black")) + geom_boxplot() + theme_bw() + 
  stat_compare_means(comparisons = my_comps) + xlab("") + ylab("Log[Burkholderia Relative Sequence Abundance]") + 
  facet_grid(~Timepoint) + theme(strip.text = element_text(size=12, face="bold")) + 
  theme(strip.background = element_rect(fill="burlywood4")) + theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = "Gadgil_Burk_AbundanceBoxPlot.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 6,
  units = c("in"),
  dpi = 300)

brady <- subset(pseqdf1, Genus == "Bradyrhizobium")
brady1 <- subset(brady, Abundance > 0)
ggplot(brady1, aes(x=Condition, y=log(Abundance), color=Condition)) + 
  scale_color_manual(values=c("red3", "black")) + geom_boxplot() + theme_bw() + 
  stat_compare_means(comparisons = my_comps) + xlab("") + ylab("log[Bradyrhizobium Sequence Abundance]") + 
  facet_grid(~Timepoint) + theme(strip.text = element_text(size=12, face="bold")) + 
  theme(strip.background = element_rect(fill="burlywood4")) + theme(legend.position = "none")


#### Aggregate Burkholderia sequences ####

burk <- subset(pseqdf1, Genus == "Burkholderia-Caballeronia-Paraburkholderia")
burk.agg <- aggregate(Abundance ~ Sample, FUN = sum, data = burk)
write.csv(burk.agg, "Burk.RelAbundance.AGG.csv")

### Correlations between bacteria and fungi datasets ###

master <- read.csv("~/Documents/Gadgil_Project/Bacteria-Fungi/Bacteria_Fungi_MASTER_DF2.csv")
library(ggplot2)
library(ggpubr)
#remove NAs
master2 <- subset(master, EcMF !="na")
master2 <- subset(master2, SoilSapros !="na")
master2$EcMF <- as.numeric(as.character(master2$EcMF))
master2$SoilSapros <- as.numeric(as.character(master2$SoilSapros))

master2$title <- "EcMF Rel. Abundance ~ Burkholderia Rel. Abundance"
master2$title2 <- "Ectomycorrhizal Fungi"
master2$title3 <- "Soil Saprotroph"

#plot bacterial richness ~ EcMF relative abundance
ggplot(master2, aes(x = EcMF, y = log(Observed))) + geom_point(size = 5, alpha = 0.7, col = "black") + 
  theme_bw() + geom_smooth(method = "lm", col = "#C42126", se = TRUE, size = 2, linetype = "dashed")  + 
  theme(axis.text.x = element_text(size = 12, face = "bold")) + 
  theme(axis.text.y = element_text(size = 12, face = "bold")) + 
  theme(axis.title.x = element_text(size = 16, face = "bold")) + 
  theme(axis.title.y = element_text(size = 16, face = "bold")) +
  ylab("Log[Observed Bacterial Richness]") + xlab("EcMF Relative Abundance") +
  theme(strip.text = element_text(size=12, face="bold")) + facet_grid(~title2) +
  annotate("text", x = 0.265, y = 6.5, label = "EcMF: p = 0.030") +
  annotate("text", x = 0.286, y = 6.3, label = "Condition: p = 0.007") +
  annotate("text", x = 0.241, y = 6.1, label = expression("Adj."~r^2 == 0.28), parse = TRUE, size = 4) +
  scale_x_continuous(labels = scales::percent)
  
#save
ggsave(
  filename = "BacterialRichnessXEcMF.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4.5,
  height = 5,
  units = c("in"),
  dpi = 300)

#plot bacterial richness ~ Soil Sapros relative abundance
ggplot(master2, aes(x = SoilSapros, y = log(Observed))) + geom_point(size = 5, alpha = 0.7, col = "dodgerblue3") + 
  theme_bw() + geom_smooth(method = "lm", col = "#C42126", se = TRUE, size = 2, linetype = "dashed")  + 
  theme(axis.text.x = element_text(size = 12, face = "bold")) + 
  theme(axis.text.y = element_text(size = 12, face = "bold")) + 
  theme(axis.title.x = element_text(size = 16, face = "bold")) + 
  theme(axis.title.y = element_text(size = 16, face = "bold")) + ylab("Log[Observed Bacterial Richness]") + xlab("Soil Saprotroph Relative Abundance") +
  theme(strip.text = element_text(size=12, face="bold")) + facet_grid(~title3) +
  annotate("text", x = 0.218, y = 6.5, label = "Soil Saprotroph: p = 0.045") +
  annotate("text", x = 0.176, y = 6.3, label = "Condition: p = 0.003") +
  annotate("text", x = 0.126, y = 6.1, label = expression("Adj."~r^2 == 0.31), parse = TRUE, size = 4) +
  scale_x_continuous(labels = scales::percent)

#save
ggsave(
  filename = "BacterialRichnessXSoilSapros.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4.5,
  height = 5,
  units = c("in"),
  dpi = 300)

# linear regression model
# bacterial richness ~ EcMF abundance 
BR_ecto <- lm(Observed ~ EcMF + Condition + Timepoint + Site + EcMF:Condition + EcMF:Timepoint + EcMF:Site, data = master2)
BR_Ecto_SS <- lm(Observed ~ EcMF + SoilSapros + Condition + Timepoint + Site, data = master2)
BR_soilsapro <- lm(Observed ~ SoilSapros + Condition + Timepoint + Site + SoilSapros:Condition + SoilSapros:Timepoint + SoilSapros:Site, data = master2)
summary(BR_ecto)
anova(BR_ecto)
anova.br.ecto.df <- data.frame(anova(BR_ecto))
write.csv(anova.br.ecto.df, "ANOVA.Table.LinearModel.BacRich-EcMF.csv")
summary(BR_soilsapro)
anova(BR_soilsapro)
anova.br.SS.df <- data.frame(anova(BR_soilsapro))
write.csv(anova.br.SS.df, "ANOVA.Table.LinearModel.BacRich-SoilSapros.csv")
# EcMF abundance ~ Burk abundance
Burk_ecto <- lm(EcMF ~ Burk + Condition + Timepoint + Site, data = master2)
Burk_alone <- lm(Burk ~ Condition + Timepoint + Site, data = master2)

# try the relaimpo package (have to change each of the categorical variables to 'factors')
library(relaimpo)
# set data in 'master2' to factors to enable the calc.relimp to work
names <- c('Condition' ,'Timepoint', 'Site')
master3 <- master2
master3[,names] <- lapply(master3[,names] , factor)
Burk_ecto.RELAMP <- calc.relimp(EcMF ~ Burk + Condition + Timepoint + Site, type = "lmg", rela = TRUE, master3)
Burk_ecto2 <- lm(EcMF ~ Burk + Condition + Timepoint + Site + SoilSapros:Burk + Burk:Timepoint + Burk:Site, data=master2)
burk.plot <- master2
burk.plot$ecto2 <- residuals(Burk_ecto)
burk.plot$ecto3 <- residuals(Burk_alone)
write.csv(burk.plot, "Burk.plot.residuals.csv")
Burk_ecto2.RELAMP <- calc.relimp(EcMF ~ Burk + Condition + Timepoint + Site + SoilSapros:Burk + Burk:Timepoint + Burk:Timepoint + Burk:Site, data=master3)
Burk_ecto2.RELAMP2 <- calc.relimp(EcMF ~ Condition + Timepoint + Site + SoilSapros:Burk + Burk:Timepoint + Burk:Timepoint + Burk:Site, data=master3)

summary(Burk_ecto)
anova(Burk_ecto)
summary(Burk_ecto2)
anova(Burk_ecto2)
anova.burk.ecto2 <- data.frame(anova(Burk_ecto2))
write.csv(anova.burk.ecto2, "ANOVA.Table.Burk.EcMF.SoilSapros.Full.Model.csv")

#plot  Ecmf abundance ~ Burk abundance
ggplot(master2, aes(y = EcMF, x = Burk)) + geom_point(size = 5, alpha = 0.9, col = "#446455") + 
  theme_bw() + geom_smooth(method = "lm", col = "#C42126", se = TRUE, size = 2, linetype = "dashed") +
  theme(axis.text.x = element_text(size = 12, face = "bold")) + 
  theme(axis.text.y = element_text(size = 12, face = "bold")) + 
  theme(axis.title.x = element_text(size = 16, face = "bold")) + 
  theme(axis.title.y = element_text(size = 16, face = "bold")) +
  xlab("Burkholderia Relative Abundance") + ylab("EcMF Relative Abundance") + 
  scale_y_continuous(labels = scales::percent) + scale_x_continuous(labels = scales::percent) +
  facet_grid(~title)+ theme(strip.text = element_text(size=8, face="bold")) +
  annotate("text", x = 0.045, y = 0.98, label = "Burk X Timepoint: p = 0.009,") +
  annotate("text", x = 0.0352, y = 0.93, label = "Condition: p = 0.006,") +
  annotate("text", x = 0.097, y = 0.98, label = expression("Adj."~r^2 == 0.04), parse = TRUE, size = 4) +
  annotate("text", x = 0.079, y = 0.93, label = expression("Adj."~r^2 == 0.05), parse = TRUE, size = 4) 

#save plot
ggsave(
  filename = "EcMFxBurk.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4.5,
  height = 5,
  units = c("in"),
  dpi = 300)

ggplot(burk.plot, aes(y = ecto2, x = ecto3)) + geom_point(size = 5, alpha = 0.9, col = "#446455") + 
  theme_bw() + geom_smooth(method = "lm", col = "#C42126", se = TRUE, size = 2, linetype = "dashed") +
  theme(axis.text.x = element_text(size = 12, face = "bold")) + 
  theme(axis.text.y = element_text(size = 12, face = "bold")) + 
  theme(axis.title.x = element_text(size = 16, face = "bold")) + 
  theme(axis.title.y = element_text(size = 16, face = "bold")) +
  xlab("Burkholderia Residuals") + ylab("EcMF Residuals") +
  facet_grid(~title)+ theme(strip.text = element_text(size=10, face="bold")) +
  annotate("text", x = -0.01, y = 0.5, label = "Burk X Timepoint: p = 0.009,") +
  annotate("text", x = -0.0185, y = 0.45, label = "Condition: p = 0.006,") +
  annotate("text", x = 0.036, y = 0.5, label = expression("Adj."~r^2 == 0.04), parse = TRUE, size = 4) +
  annotate("text", x = 0.019, y = 0.45, label = expression("Adj."~r^2 == 0.05), parse = TRUE, size = 4)


#save
ggsave(
  filename = "EcMFxBurk_Residuals.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 4.5,
  height = 5,
  units = c("in"),
  dpi = 300)


# Check how EcMF:SoilSapro ratios predict bacterial richness
master2$EcMF.SoilSap <- master2$EcMF / master2$SoilSapros
#plot bacterial richness ~ Soil Sapros relative abundance
ggplot(master2, aes(x = log(EcMF.SoilSap), y = log(Observed))) + geom_point(size = 5, alpha = 0.7, col = "dodgerblue3") + 
  theme_bw() + geom_smooth(method = "lm", col = "#C42126", se = TRUE, size = 2, linetype = "dashed")  + 
  theme(axis.text.x = element_text(size = 12, face = "bold")) + 
  theme(axis.text.y = element_text(size = 12, face = "bold")) + 
  theme(axis.title.x = element_text(size = 16, face = "bold")) + 
  theme(axis.title.y = element_text(size = 16, face = "bold")) + ylab("Log[Observed Bacterial Richness]") + xlab("EcMF:Soil Saprotroph Ratio") +
  theme(strip.text = element_text(size=12, face="bold"))
# linear model for the ratios
BR_EcMF.SS <- lm(Observed ~ EcMF.SoilSap + Condition + Timepoint + 
                   Site + EcMF.SoilSap:Condition + EcMF.SoilSap:Timepoint + 
                   EcMF.SoilSap:Site, data = master2)
summary(BR_EcMF.SS)
anova(BR_EcMF.SS)
#save table
BR.EcMF.SS.Model <- data.frame(anova(BR_EcMF.SS))
write.csv(BR.EcMF.SS.Model, "BacRich.EcMF.SS.Model.csv")
# Test the R squared values 
master2[,names] <- lapply(master2[,names] , factor)
BR_EcMF.SS.RELAMP <- calc.relimp(Observed ~ EcMF.SoilSap + Condition + Timepoint + 
                   Site + EcMF.SoilSap:Condition + EcMF.SoilSap:Timepoint + 
                   EcMF.SoilSap:Site, type = "lmg", data = master2)
BR_EcMF.SS.RELAMP
