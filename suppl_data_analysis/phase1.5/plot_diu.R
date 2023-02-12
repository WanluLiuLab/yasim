# Import the required packages
library(IsoformSwitchAnalyzeR)
library(tidyverse)

sg1 <- read.table("ce11_pbsim2_P4C2_diu1.fq.bam.fc.txt", sep = "\t", header = T)
sg1 <- sg1[, c(1, 7)]
sg2 <- read.table("ce11_pbsim2_P4C2_diu2.fq.bam.fc.txt", sep = "\t", header = T)
sg2 <- sg2[, c(1, 7)]
#sg3<-read.table("na_rep3Aligned.out.sam.bam.txt",sep = "\t",header = T)
#sg3<-sg3[,c(1,7)]
w1 <- read.table("ce11_pbsim3_ONT_diu1.fq.bam.fc.txt", sep = "\t", header = T)
w1 <- w1[, c(1, 7)]
w2 <- read.table("ce11_pbsim3_ONT_diu2.fq.bam.fc.txt", sep = "\t", header = T)
w2 <- w2[, c(1, 7)]
#w3<-read.table("pr_rep3Aligned.out.sam.bam.txt",sep = "\t",header = T)
#w3<-w3[,c(1,7)]
merge <- dplyr::full_join(sg1, sg2, by = "Geneid")
#merge<-dplyr::full_join(merge,sg3,by="Geneid")
merge <- dplyr::full_join(merge, w1, by = "Geneid")
merge <- dplyr::full_join(merge, w2, by = "Geneid")
#merge<-dplyr::full_join(merge,w3,by="Geneid")
colnames(merge) <- c("isoform_id", "P4C2_1", "P4C2_2", "ONT_1", "ONT_2")
merge <- replace(merge, is.na(merge), 0)
#write.csv(merge,"./Mix_DIU_merge.csv")
sampleID <- c("P4C2_1", "P4C2_2", "ONT_1", "ONT_2")
condition <- c("P4C2", "P4C2", "ONT", "ONT")
designMatrix <- cbind(data.frame(sampleID), data.frame(condition))

### Create switchAnalyzeRlist
aSwitchList <- importRdata(
    isoformNtFasta = "stringtie_P4C2_pbsim3ONT.fa",
    showProgress = FALSE,
    isoformCountMatrix = merge,
    designMatrix = designMatrix,
    isoformExonAnnoation = "merge_P4C2_pbsim3ONT.gtf",
    ignoreAfterPeriod = FALSE #if using the reference gtf, then this should set TRUE.
)


#SwitchListFiltered <- preFilter(
# switchAnalyzeRlist = aSwitchList,
#geneExpressionCutoff = 5,
#isoformExpressionCutoff = 5,
#removeSingleIsoformGenes = TRUE)

SwitchListAnalyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = aSwitchList,
    reduceToSwitchingGenes = TRUE,
    reduceFurtherToGenesWithConsequencePotential = FALSE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    onlySigIsoforms = FALSE
)

switchListO <- analyzeORF(SwitchListAnalyzed,
                          showProgress = FALSE)
switchListS <- extractSequence(switchListO,
)
switchListsR <- analyzeAlternativeSplicing(switchListS, quiet = TRUE, onlySwitchingGenes = FALSE)
consequencesOfInterest <- c('intron_retention', 'NMD_status', 'ORF_seq_similarity')

exampleSwitchListAnalyzed <- analyzeSwitchConsequences(
    switchListsR,
    consequencesToAnalyze = consequencesOfInterest,
    dIFcutoff = 0.1,
    alpha = 0.05,
    showProgress = FALSE
)

##Obtain all switching genes
All_significant_DIU <- extractTopSwitches(
    exampleSwitchListAnalyzed,
    filterForConsequences = TRUE,
    n = NA,
    extractGenes = FALSE,
    sortByQvals = TRUE
)
