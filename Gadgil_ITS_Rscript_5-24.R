################## Gadgil ITS R Script #################
##### Louis Berrios #####

library(dada2)
library(Biostrings)
library(ShortRead)
library(microbiome)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(metacoder)
library(phyloseq)

path <- "~/Documents/Gadgil_Project/Gadgil_ITS"
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

# Identify primers - used Hiro's PCR primers
FWD <- "THGGTCATTTAGAGGAASTAA"  ## CHANGE ME to your forward primer seq
REV <- "TTYRCTRCGTTCTTCATC"  ## CHANGE ME...


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

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 8, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on Mac, set multithread = TRUE
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

saveRDS(seqtab,"Gadgil-ITS1.seqtab.R1.RDS")
saveRDS(seqtab.nochim,"Gadgil-ITS1.seqtab.nochim.R1.RDS")

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
write.csv(track2,"Gadgil-ITS1_fungi_summary.csv")

rownames(track2) <- sample.names
head(track2)


####### Assign Taxonomy #########


seqtab.nochim<-readRDS("Gadgil-ITS1.seqtab.nochim.R1.RDS")

# Assign taxonomy using the UNITE database
unite.ref <- "/Users/louisberrios/Documents/Gadgil_Project/Gadgil_ITS/UNITE2022/sh_general_release_dynamic_29.11.2022_dev.fasta" # need to change
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Save OTU table and taxonomic table as RDS files
# to hand off to dada2_to_phyloseq.R
saveRDS(seqtab.nochim, "Gadgil-ITS1-FUNGI_seqtab_nochim_taxa.rds")
saveRDS(taxa, "Gadgil-ITS1-FUNGI-TAXA.rds")



####### Create the phyloseq object ########

#read in the taxonomy & OTU table files that were created from DADA2

seqtab.nochim<-readRDS('Gadgil-ITS1-FUNGI_seqtab_nochim_taxa.rds')
taxa<-readRDS('Gadgil-ITS1-FUNGI-TAXA.rds')

# create sample data 
names <- row.names(seqtab.nochim)
sam <- as.data.frame(names)
sam2 <- sam %>% 
       mutate(Condition = ifelse(grepl("-\\T\\-", names), "Trenched", "Untrenched"))
sam3 <- sam2 %>% 
       mutate(Timepoint = ifelse(grepl("T2", names), "T2", 
                          ifelse(grepl("T3", names), "T3", 
                          ifelse(grepl("T4", names), "T4", "T1"))))
colnames(sam3)[1] <- "sample_names"
write.csv(sam3, "Gadgil_ITS_SampleData.csv") # tweaked in excel
# import sample_data table
sample <- read.csv("Gadgil_ITS_SampleData.csv", row.names = 1)
SAM <- sample_data(sample, errorIfNULL = T)

#add functional guild data
fungal.traits.database<-read.csv('FungalTraits.csv')
tax.table<-data.frame(taxa)
Genus<-gsub("g__","",tax.table$Genus)
traits_table<-fungal.traits.database[match(Genus,fungal.traits.database$GENUS),]
tax.trait.table<-cbind(tax.table,traits_table)

#create a phyloseq object
library(phyloseq)
fun.ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(SAM), tax_table(as.matrix(tax.trait.table)))

#prune taxa + samples 
ps1 = prune_taxa(taxa_sums(fun.ps) > 100, fun.ps) 
ps1 = prune_samples(sample_sums(ps1)>5000, ps1)
ps1 <- subset_samples(ps1, Condition !="NA")

#save phyloseq object
saveRDS(ps1,"Gadgil-ITS1.rds")

#rarefy the dataset
ps.fun.rare <- rarefy_even_depth(ps1, sample.size = min(sample_sums(ps1)),
                                 rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

raredf <- psmelt(ps.fun.rare)

#save rarefied dataset
saveRDS(ps.fun.rare, "Gadgil-Fungi-PS-Rarefied.rds")

# make compositional 
pseq <- microbiome::transform(ps.fun.rare, "compositional")
pseq_df <- psmelt(pseq)

saveRDS(pseq, "Gadgil-Fungi-Rarefied-Transf.rds")

pseq_final <- readRDS("Gadgil-Fungi-Rarefied-Transf.rds")
pseq_final_df <-psmelt(pseq_final)

# subset out ectomycorrhizal ids 
myco_final <- subset_taxa(pseq_final, primary_lifestyle == "ectomycorrhizal")
# make into data frame 
myco_df_final <- psmelt(myco_final)
# remove NAs
myco_df1_final <- subset(myco_df_final, Condition !="NA")

# save rarefied dataframe
write.csv(myco_df1, "Myco_Gadgil_Rarefied_DF.csv")
# abundance by condition and time
library(ggplot2)
myco_df1$title <- "EcMF Relative Sequence Abundance"
# aggregate by sample 
myco_df.agg <- aggregate(Abundance ~ Sample + Timepoint + Condition, FUN=sum, data = myco_df1_final)
myco_df.agg$title <- "EcMF Relative Sequence Abundance"
# ANOVA
myco_aov <- aov(Abundance ~ Condition + Timepoint + Condition:Timepoint, data=myco_df1)
myco_aov2 <- aov(Abundance ~ Condition + Timepoint + Condition:Timepoint, data=myco_df.agg) # aggregated by sample
summary(myco_aov) 
summary(myco_aov2)
ggplot(myco_df.agg, aes(x=Condition, y=log(Abundance), color = Condition)) + 
  geom_boxplot() + theme_bw()  + theme(axis.title.x = element_blank()) +
  scale_color_manual(values=c("firebrick3", "gray1")) + theme(legend.position = "none") +
  facet_grid(~title) + theme(strip.text = element_text(size=12, face="bold")) + 
  theme(strip.background = element_rect(fill = "burlywood4")) + 
  theme(axis.text.x = element_text(size=12, face="bold")) + ylab("Log[EcMF Sequence Abundance]") +
  annotate(geom = "text", x = 0.99, y = 0.3, size = 5, label = paste("Condition:F[\"1,157\"] ==", 10.30), parse=TRUE) +
  annotate(geom = "text", x=1.75, y=0.3, size = 5, label = ", p = 0.0016") +
  annotate(geom = "text", x = 0.99, y = 0.08, size = 5, label = paste("Timepoint:F[\"3,155\"] ==", 2.91), parse=TRUE) +
  annotate(geom = "text", x=1.75, y=0.08, size = 5, label = ", p = 0.0363")

# save plot
ggsave(
  filename = "Gadgil_ITS_MycoAbundanceBoxPlot.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 6,
  units = c("in"),
  dpi = 300)
# compositional plot of all fungal guilds across time
pseqdf2 <- subset(pseq_df, primary_lifestyle !="NA")
pseqdf2 <- subset(pseqdf2, Condition !="NA")
guild.agg <- aggregate(Abundance ~ primary_lifestyle + Condition + Timepoint, FUN=sum, data = pseqdf2)
pseq.guilds <- pseqdf2
pseq.guilds$primary_lifestyle[pseq.guilds$primary_lifestyle == "animal_parasite" | 
                                pseq.guilds$primary_lifestyle == "epiphyte" |
                                pseq.guilds$primary_lifestyle == "foliar_endophyte" |
                                pseq.guilds$primary_lifestyle == "lichen_parasite" |
                                pseq.guilds$primary_lifestyle == "lichenized" |
                                pseq.guilds$primary_lifestyle == "mycoparasite" | 
                                pseq.guilds$primary_lifestyle == "plant_pathogen" |
                                pseq.guilds$primary_lifestyle == "root_endophyte" |
                                pseq.guilds$primary_lifestyle == "sooty_mold" |
                                pseq.guilds$primary_lifestyle == "nectar/tap_saprotroph" |
                                pseq.guilds$primary_lifestyle == "dung_saprotroph"] <- "Other Guilds"

ggplot(pseq.guilds, aes(fill=primary_lifestyle, y=Abundance, x=Condition)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + 
  scale_fill_manual(values = c("darkgoldenrod3", "gray", "black", "darkslategray4",
                               "#A44A3F", "burlywood4")) + facet_grid(~Timepoint) + 
  theme(axis.text.x = element_text(angle=45)) + scale_y_continuous(labels = scales::percent) + 
  theme(axis.text.x = element_text(hjust = 1)) + labs(fill = "Fungal Guilds") + theme(axis.title.x = element_blank())
# save
ggsave(
  filename = "Gadgil_FungalGuildAbundanceCompo2.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 7,
  height = 5.5,
  units = c("in"),
  dpi = 300)

# compositional plot of myco genera across time
myco.agg <- aggregate(Abundance ~ GENUS, FUN = sum, data=myco_df1) # determine the top genera
myco.agg2 <- aggregate(Abundance ~ Sample, FUN = sum, data = myco_df1) # to be used in correlations and lm models
write.csv(myco.agg2, "Ecto_RelAbundance.Sample.csv")
# plot myco
ggplot(myco_df1, aes(fill=GENUS, y=Abundance, x=Condition)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + 
  scale_fill_manual(values = c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236",
                                 "#EABE94", "#00A08A", "#F2AD94", "#C6CDF7",
                                 "#D4E6B5", "#C6C2B6", "#793A2B", "#CCBA72",
                                 "darkslategray4", "#9B9C5A", "#AB947A", "#816C5B",
                                 "#59351F", "#FAFAFA", "lightpink3", "#D9D9D9",
                                 "#7294D4", "#91B3D7", "#B9975B", "#A89F68",
                                 "#7F9C6C", "#A44A3F", "#DD6031", "#FBCB7B",
                                 "gold", "red3")) + facet_grid(~Timepoint) + 
  theme(axis.text.x = element_text(angle=45)) + scale_y_continuous(labels = scales::percent) + 
  theme(axis.text.x = element_text(hjust = 1)) + labs(fill = "Ectomycorrhizal Genera") + theme(axis.title.x = element_blank())
ggsave(
  filename = "Gadgil_ITS_MycoAbundanceCompo.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 7,
  height = 5.5,
  units = c("in"),
  dpi = 300)

#subset saprotrophs
sapros <- subset_taxa(pseq_final, primary_lifestyle == "dung_saprotroph" | primary_lifestyle == "soil_saprotroph" | primary_lifestyle == "litter_saprotroph" | primary_lifestyle == "unspecified_saprotroph" | primary_lifestyle == "wood_saprotroph")
sapros_df <- psmelt(sapros)

# remove NAs
sapros_df1 <- subset(sapros_df, Condition !="NA")
# aggregate by order to determine the top 10 percent
sapros.agg <- aggregate(Abundance ~ Class.1, FUN = sum, data=sapros_df1)
sapros.agg.genus <- aggregate(Abundance ~ Sample + GENUS + Condition + Timepoint, FUN=sum, data=sapros_df1)
sapros.agg.genus.TABLESX <- aggregate(Abundance ~ GENUS + Condition + Timepoint + primary_lifestyle, FUN=sum, data=sapros_df1)
sapros.agg.genus.TABLESX2 <- subset(sapros.agg.genus.TABLESX, Abundance > 0)
write.csv(sapros.agg.genus.TABLESX2, "SAPROS.Genus.Agg.Table.csv")


sapros.agg.genus2 <- subset(sapros.agg.genus, Abundance > 0)

sapros_df2 <- sapros_df1 # duplicate
# group classes that total less than 5% into 'Other Orders'
sapros_df2$Class.1[sapros_df2$Class.1 == "Calcarisporiellomycetes" | 
                     sapros_df2$Class.1 == "Orbiliomycetes" |
                     sapros_df2$Class.1 == "Tritirachiomycetes" |
                     sapros_df2$Class.1 == "Atractiellomycetes" |
                     sapros_df2$Class.1 == "Pezizomycetes"] <- "Classes < 5% Total Relative Abundance"

                          
# compositional plot of sapros across time
ggplot(sapros_df2, aes(fill=Class.1, y=Abundance, x=Condition)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + 
  scale_fill_manual(values = c(
    "#5f4b42",
    "#b8575f",
    "#d38168",
    "#d5bf73",
    "#30949e",
    "#4e6189",
    "#8c98b0",
    "#2e5e60",
    "violetred4",
    "#96b9a9",
    "#e6dbb2",
    "#565655", "#5b8656")) + facet_grid(primary_lifestyle~Timepoint) +
  theme(axis.text.x = element_text(angle=45)) + scale_y_continuous(labels = scales::percent) + 
  theme(axis.text.x = element_text(hjust = 1)) + labs(fill = "Saprotrophic Guild") + 
  theme(axis.title.x = element_blank())
ggsave(
  filename = "Gadgil_ITS-SAPROS_AbundanceCompo.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 7,
  height = 5.5,
  units = c("in"),
  dpi = 300)

# soil saprotrophs only
soil_sapros <- subset_taxa(pseq_final, primary_lifestyle == "soil_saprotroph")
soil.sapros.df <- psmelt(soil_sapros)
soil.sapros.agg <- aggregate(Abundance ~ GENUS + Condition + Timepoint, FUN = sum, data=soil.sapros.df)
soil.sapros.sample <- aggregate(Abundance ~ Sample, FUN=sum, data = soil.sapros.df)
write.csv(soil.sapros.sample, "SoilSapros.Abundance.csv")
# compositional plot of soil sapros across time
ggplot(soil.sapros.agg, aes(fill=GENUS, y=Abundance, x=Condition)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + 
  scale_fill_manual(values = c(
    "#5f4b42",
    "#b8575f",
    "#d38168",
    "dodgerblue3",
    "#30949e",
    "#4e6189",
    "#8c98b0",
    "#2e5e60",
    "violetred4",
    "#96b9a9", "#00A08A",
    "#565655", "#5b8656", "#EABE94", "#e6dbb2", "#F2AD94", "#C6CDF7",
    "#D4E6B5", "#C6C2B6", "#793A2B", "red3",
    "darkslategray4", "#9B9C5A", "#AB947A", "#816C5B",
    "black", "#FAFAFA")) + 
  theme(axis.text.x = element_text(angle=45)) + scale_y_continuous(labels = scales::percent) + 
  theme(axis.text.x = element_text(hjust = 1)) + labs(fill = "Soil Saprotroph Genus") + 
  theme(axis.title.x = element_blank()) + facet_grid(~Timepoint)
ggsave(
  filename = "Gadgil_SOIL-SAPROS_AbundanceCompo.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 7,
  height = 5.5,
  units = c("in"),
  dpi = 300)


sapros_df1$title <- "Saprotropic Fungi"
# aggregate by sample 
sapros_df1_agg <- aggregate(Abundance ~ Sample + Condition + Timepoint, FUN=sum, data = sapros_df1)
sapros_df1_agg$title <- "Saprotropic Fungi"
# plot
ggplot(sapros_df1_agg, aes(x=Condition, y=log(Abundance), color = Condition)) + 
  geom_boxplot() + theme_bw() + facet_grid(~title) + 
  theme(axis.title.x = element_blank()) + scale_color_manual(values=c("dodgerblue4", "gray1")) + 
  theme(legend.position = "none")
sapros_aov <- aov(Abundance ~ Condition + Timepoint + Condition:Timepoint, data = sapros_df1)
summary(sapros_aov)
sapros_aov.agg <- aov(Abundance ~ Condition + Timepoint + Condition:Timepoint, data = sapros_df1_agg)
summary(sapros_aov.agg)
# plot
ggplot(sapros_df1_agg, aes(x=Condition, y=log(Abundance), color = Condition)) + 
  geom_boxplot() + theme_bw()  + theme(axis.title.x = element_blank()) +
  scale_color_manual(values=c("dodgerblue4", "gray1")) + theme(legend.position = "none") +
  facet_grid(~title) + theme(strip.text = element_text(size=12, face="bold", color = "white")) + 
  theme(strip.background = element_rect(fill = "black")) + 
  theme(axis.text.x = element_text(size=12, face="bold")) + ylab("Log[Saprotrophic Fungi Sequence Abundance]") +
  annotate(geom = "text", x = 1.0, y = 1, size = 5, label = paste("Condition:F[\"1,157\"] ==", 7.96), parse=TRUE) +
  annotate(geom = "text", x=1.75, y=1, size = 5, label = ", p = 0.005") +
  annotate(geom = "text", x = 1, y = 0.79, size = 5, label = paste("Timepoint:F[\"3,155\"] ==", 4.29), parse=TRUE) +
  annotate(geom = "text", x=1.75, y=0.79, size = 5, label = ", p = 0.006")
#save
ggsave(
  filename = "Gadgil_ITS_SaprosAbundanceBoxPlot.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5,
  height = 6,
  units = c("in"),
  dpi = 300)

# NMDS ordination
#perform an ordination
# remove controls
pseq2 <- subset_samples(pseq_final, Condition != "NA") # ALL FUNGI
pseq3 <- subset_taxa(pseq2, primary_lifestyle == "ectomycorrhizal") # ECMF
pseq4 <- subset_taxa(pseq2, primary_lifestyle == "dung_saprotroph" |
                       primary_lifestyle == "litter_saprotroph" |
                       primary_lifestyle == "soil_saprotroph" |
                       primary_lifestyle == "unspecified_saprotroph" |
                       primary_lifestyle == "wood_saprotroph")

ordination <- ordinate(pseq2, "NMDS", "bray", k=5) # stress = 0.1060
ordination2 <- ordinate(pseq3, "NMDS", "bray", k=5) # stress = 0.1074
ordination3 <- ordinate(pseq4, "NMDS", "bray", k=3) # stress = 0.1114
#plot ordination
# all fungi
sample_data(pseq2)$BETA <- "Fungal Beta Diversity"
plot_ordination(pseq2, ordination, type="sample",color="Timepoint", shape = "Condition", axes = 2:3) + 
  theme_bw() + scale_color_manual(values = c("#9a031e", "#0f4c5c", "#F8AFA8", "#f48c06")) + 
  geom_point(size = 4, alpha = 1) + stat_ellipse(type ="norm", linetype = 2) +
  labs(title = "Whole Fungal Community") +
  annotate("text", x = -1.33, y = 3.1, label = "Stress = 0.1060") +
  annotate("text", x = -1.23, y = 2.8, label = "Condition: p = 0.026") +
  annotate("text", x = -1.2, y = 2.5, label = "Timepoint: p = 0.005") + facet_grid(~BETA) +
  theme(strip.text = element_text(size = 12, face = "bold")) + 
  theme(strip.background = element_rect(fill = "#c18c5d"))
# save all fungi plot
ggsave(
  filename = "Gadgil_AllFungi_BETADIVERSITY.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 5,
  units = c("in"),
  dpi = 300)

# saprotrophs
sample_data(pseq4)$BETA <- "Fungal Saprotroph Beta Diversity"
plot_ordination(pseq4, ordination, type="sample",color="Timepoint", shape = "Condition", axes = 2:3) + 
  theme_bw() + scale_color_manual(values = c("#9a031e", "#0f4c5c", "#F8AFA8", "#f48c06")) + 
  geom_point(size = 4, alpha = 1) + stat_ellipse(type ="norm", linetype = 2) +
  labs(title = "Saprotrophic Fungal Community") +
  annotate("text", x = -1.33, y = 3.1, label = "Stress = 0.1114") +
  annotate("text", x = -1.23, y = 2.8, label = "Condition: p = 0.007") +
  annotate("text", x = -1.2, y = 2.5, label = "Timepoint: p = 0.001") + facet_grid(~BETA) +
  theme(strip.text = element_text(size = 12, face = "bold")) + 
  theme(strip.background = element_rect(fill = "#c18c5d"))
# save sapros fungi plot
ggsave(
  filename = "Gadgil_SAPROS_BETADIVERSITY.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 5,
  units = c("in"),
  dpi = 300)

# EcM fungi
sample_data(pseq3)$BETA <- "EcMF Beta Diversity"
plot_ordination(pseq3, ordination2, type="sample",color="Timepoint", shape = "Condition", axes = 2:3) + 
  theme_bw() + scale_color_manual(values = c("#9a031e", "#0f4c5c", "#F8AFA8", "#f48c06")) + 
  geom_point(size = 4, alpha = 1) + stat_ellipse(type ="norm", linetype = 2) +
  labs(title = "Ectomycorrhizal Community") +
  annotate("text", x = -1.33, y = 3.1, label = "Stress = 0.1074") +
  annotate("text", x = -1.23, y = 2.8, label = "Condition: p = 0.064") +
  annotate("text", x = -1.2, y = 2.5, label = "Timepoint: p = 0.711") + facet_grid(~BETA) +
  theme(strip.text = element_text(size = 12, face = "bold")) + 
  theme(strip.background = element_rect(fill = "#c18c5d"))

# Save ECM plot
ggsave(
  filename = "Gadgil_ECM_BETADIVERSITY.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 5,
  units = c("in"),
  dpi = 300)

####### Perform PERMANOVA ######
#calculate bray-curtis distance matrix
library(vegan)
myco.ord.bac.bray <- phyloseq::distance(pseq3, method = "bray", k=5) # ECMF
FUNGI.ord.bac.bray <- phyloseq::distance(pseq2, method = "bray", k=5) # ALL FUNGI
sapros.ord.bac.bray <- phyloseq::distance(pseq4, method = "bray", k=3) # Saprotrophic FUNGI

#make a dataframe from the sample_data
myco.ord.bac.bray.sampleDF <- data.frame(sample_data(pseq3)) # ECMF
FUNGI.ord.bac.bray.sampleDF <- data.frame(sample_data(pseq2)) # ALL FUNGI
sapros.ord.bac.bray.sampleDF <- data.frame(sample_data(pseq4)) # Sapros

#adonis test
PERMANOVA <- adonis2(myco.ord.bac.bray ~ Condition + Timepoint + Condition:Timepoint, data = myco.ord.bac.bray.sampleDF) # ECTOS
PERMANOVA2 <- adonis2(FUNGI.ord.bac.bray ~ Condition + Timepoint + Condition:Timepoint, data = FUNGI.ord.bac.bray.sampleDF) # ALL FUNGI
PERMANOVA3 <- adonis2(sapros.ord.bac.bray ~ Condition + Timepoint + Condition:Timepoint, data = sapros.ord.bac.bray.sampleDF) # sapros

#write PERMANOVA results to table CSV
PERMANOVA.df <- data.frame(PERMANOVA) # ectos
PERMANOVA2.df <- data.frame(PERMANOVA2) # all fungi
PERMANOVA3.df <- data.frame(PERMANOVA3) # sapros
write.csv(PERMANOVA.df, "Gad_Myco_Permanova.csv")
write.csv(PERMANOVA2.df, "Gad_ALL-FUNGI_Permanova.csv")
write.csv(PERMANOVA3.df, "Gad_SAPROS_Permanova.csv")

#run betadisper function on distance matrix
# ECTOS
beta_adonis_myco_condition <- betadisper(myco.ord.bac.bray, myco.ord.bac.bray.sampleDF$Condition, bias.adjust = TRUE) # myco_condition
beta_adonis_myco_timepoint <- betadisper(myco.ord.bac.bray, myco.ord.bac.bray.sampleDF$Timepoint, bias.adjust = TRUE) # myco_timepoint
#stats
stat_disp_anova_myco <- anova(beta_adonis_myco_condition)
stat_disp_anova_myco2 <- anova(beta_adonis_myco_timepoint)
plot(beta_adonis_myco_condition, hull=FALSE, ellipse=TRUE)
plot(beta_adonis_myco_timepoint, hull=FALSE, ellipse=TRUE)
#pairwise permutation test for homogeneity of multivariate dispersions
permutest(beta_adonis_myco_condition, pairwise = TRUE, permutations = 1000)
permutest(beta_adonis_myco_timepoint, pairwise = TRUE, permutations = 1000)
# ALL FUNGI
beta_adonis_allfungi_condition <- betadisper(FUNGI.ord.bac.bray, FUNGI.ord.bac.bray.sampleDF$Condition, bias.adjust = TRUE) # all fungi_condition
beta_adonis_allfungi_timepoint <- betadisper(FUNGI.ord.bac.bray, FUNGI.ord.bac.bray.sampleDF$Timepoint, bias.adjust = TRUE) # all fungi_timepoint
#stats
stat_disp_anova_AF <- anova(beta_adonis_allfungi_condition)
stat_disp_anova2_AF <- anova(beta_adonis_allfungi_timepoint)
plot(beta_adonis_allfungi_condition, hull=FALSE, ellipse=TRUE)
plot(beta_adonis_allfungi_timepoint, hull=FALSE, ellipse=TRUE)
#pairwise permutation test for homogeneity of multivariate dispersions
permutest(beta_adonis_allfungi_condition, pairwise = TRUE, permutations = 1000)
permutest(beta_adonis_allfungi_timepoint, pairwise = TRUE, permutations = 1000)
# SAPROS
beta_adonis_sapros_condition <- betadisper(sapros.ord.bac.bray, sapros.ord.bac.bray.sampleDF$Condition, bias.adjust = TRUE) # sapros_condition
beta_adonis_sapros_timepoint <- betadisper(sapros.ord.bac.bray, sapros.ord.bac.bray.sampleDF$Timepoint, bias.adjust = TRUE) # sapros_timepoint
#stats
stat_disp_anova_SAPROS <- anova(beta_adonis_sapros_condition)
stat_disp_anova2_SAPROS <- anova(beta_adonis_sapros_timepoint)
plot(beta_adonis_sapros_condition, hull=FALSE, ellipse=TRUE)
plot(beta_adonis_sapros_timepoint, hull=FALSE, ellipse=TRUE)
#pairwise permutation test for homogeneity of multivariate dispersions
permutest(beta_adonis_sapros_condition, pairwise = TRUE, permutations = 1000)
permutest(beta_adonis_sapros_timepoint, pairwise = TRUE, permutations = 1000)

###### Alpha Diversity / Observed Richness ####### 
#use non-transformed data (ps.fun.rare)
# try for just ectos and sapros
ad <- estimate_richness(ps.fun.rare) # all fungi
write.csv(ad, "Fungi_Alpha_Diversity.csv")
AD_MOD <- read.csv("Fungi_Alpha_Diversity.csv")
AD_MOD$title <- "Fungal Richness"
aov_ad.div <- aov(Observed ~ Condition + Timepoint + Condition:Timepoint, AD_MOD)
summary(aov_ad.div) 
ggplot(AD_MOD, aes(x=Condition, y=Observed, color=Condition)) + geom_boxplot() + theme_bw() + 
  scale_color_manual(values=c("firebrick3", "gray1")) + xlab("") + 
  theme(axis.text.x = element_text(size = 12)) + 
  theme(legend.position = "none") + theme(axis.text.y = element_text(size = 12)) + 
  facet_grid(~title) + theme(strip.text = element_text(size = 14, face = "bold")) + 
  theme(axis.title.y = element_text(size = 16, vjust = 1.5)) + facet_grid(~Timepoint)

  theme(strip.background = element_rect(fill = "lightblue4")) + annotate(geom = "text", x = 1.04, y = 210, size = 5, label = paste("Condition:F[\"1,151\"] ==", 3.244), parse=TRUE) +
  annotate(geom = "text", x=1.71, y=210, size = 5, label = ", p = 0.0737") +
  annotate(geom = "text", x = 1.02, y = 200, size = 5, label = paste("Timepoint:F[\"3,151\"] ==", 1.777), parse=TRUE) +
  annotate(geom = "text", x=1.71, y=200, size = 5, label = ", p = 0.1539")

# Save plot
ggsave(
  filename = "Gadgil_ITS_ObservedRichnessBoxPlot.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 6,
  height = 5,
  units = c("in"),
  dpi = 300)

