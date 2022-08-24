debug <- FALSE
if (debug){
  setwd("Y:/")
  DB <- "NCBIGenomes06"
  SAMPLE <- "L36_M6"
  DIR <- "/Infections/Runs/V1/Samples/L36_M6"
  REF <- "/Infections/Data/GenomeDB/NCBIGenomes06/LUTGrammy/taxids_names_lengths_tax.tab"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  DIR <- args[1]
  SAMPLE <- args[2]
  REF <- args[3]
  STATS <- args[4]
  OUTPUT = args[5]
}


grammy.file <- paste(DIR, "/", SAMPLE, ".tab", sep = "")
grammy.tab <- read.table(grammy.file, header = FALSE, fill = TRUE)
colnames(grammy.tab) <- c("SAMPLE", "Taxid", "GrAb", "GrEr")
#
blast.file <- paste(DIR, "/", SAMPLE,".tblat.1", sep = "")
blast <- read.table(blast.file, header = FALSE, fill = TRUE)
total.blast <- nrow(blast)
#
align.stats.file <- paste(STATS, SAMPLE, "_stats.align.tab",  sep = "")
align.stats <- read.table(align.stats.file, header = FALSE, fill = TRUE)
hg.coverage <- align.stats[align.stats$V2 == "chr21_coverage",]$V3
#
grammy.LUT <- read.table(REF, header = TRUE, fill = TRUE)

grammy.tab.info <- merge(grammy.tab, grammy.LUT, by = "Taxid")
# weighted genome size
grammy.tab.info$hgcoverage <- hg.coverage
grammy.tab.info$WeightedGenome <- sum(grammy.tab.info$Length * grammy.tab.info$GrAb)
grammy.tab.info$AdjustedBlast <- total.blast*(grammy.tab.info$Length*grammy.tab.info$GrAb/grammy.tab.info$WeightedGenome)
grammy.tab.info$Coverage <- 75*grammy.tab.info$AdjustedBlast/grammy.tab.info$Length
grammy.tab.info$RelCoverage <- 2*75*grammy.tab.info$AdjustedBlast/grammy.tab.info$Length/hg.coverage
#
output.file <- paste(OUTPUT, sep = "")
write.table(grammy.tab.info, output.file,sep ="\t", row.names = FALSE)
