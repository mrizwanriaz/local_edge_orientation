library(MASS)	
library(class)
library(cluster)	
library(impute)
library(dplyr)
library(tibble)
library(GeneCycle)
library(GeneNet)
library(sem)
source("neo.SMA.txt")
source("sma_package.R") # Code from the no longer supported sma package.
#Load network functions 
source("NetworkFunctions_Jan2010.R"); 
source("CausalityFunctions.R")

dummy_samples <- c("samp1","samp2","samp3","samp4","samp5","samp6","samp7","samp8","samp9","samp10")

#Read in the SNPs data file
datSNP <- readRDS("SNPs.RDS")[c(1:10),]

# Display the column names after the transformations
colnames(datSNP)

# Read module eigengenes (from WGCNA output)
MEs <- readRDS("MEs.RDS")

# Read traits data 
datTraits <- readRDS("traits.RDS")[c(1:10),]


# Calculate number of samples, SNPs, traits, and modules
SampNum <- length(dummy_samples)  # number of samples
SNP_Num <- ncol(datSNP)  # number of SNPs (columns in SNP data)
Traits_Num <- ncol(datTraits)  # number of traits (columns in trait data)
MEs_Num <- ncol(MEs)  # number of modules (columns in MEs data)



####### causality analysis using LEO (single.marker.analysis).

# Now we specify the 5 single anchor models
# specify.model is now a synonym for specifyModel in the package sem
# Single anchor model 1: M1->A->B
CausalModel1=specifyModel(text="
                          A -> B, betaAtoB, NA
                          M1 -> A, gammaM1toA, NA
                          A <-> A, sigmaAA, NA
                          B <-> B, sigmaBB, NA
                          M1 <-> M1, sigmaM1M1,NA                        
                          ")


# Single anchor model 2: M1->B->A
CausalModel2=specifyModel(text="
                          B -> A, betaBtoA, NA
                          M1 -> B, gammaM1toB, NA
                          A <-> A, sigmaAA, NA
                          B <-> B, sigmaBB, NA
                          M1 <-> M1, sigmaM1M1,NA                     
                          ")

# Single anchor model 3: A<-M1->B
CausalModel3=specifyModel(text="
                          M1 -> A, gammaM1toA, NA
                          M1 -> B, gammaM1toB, NA
                          A <-> A, sigmaAA, NA
                          B <-> B, sigmaBB, NA
                          M1 <-> M1, sigmaM1M1,NA")


# Single anchor model 4: M1->A<-B
CausalModel4=specifyModel(text="
                          M1 -> A, gammaM1toA, NA
                          B -> A, gammaBtoA, NA
                          A <-> A, sigmaAA, NA
                          B <-> B, sigmaBB, NA
                          M1 <-> M1, sigmaM1M1,NA")


# Single anchor model 5: M1->B<-A
CausalModel5=specifyModel(text="
                          M1 -> B, gammaM1toB, NA
                          A -> B, gammaAtoB, NA
                          A <-> A, sigmaAA, NA
                          B <-> B, sigmaBB, NA
                          M1 <-> M1, sigmaM1M1,NA
                          ")


# Function definitions for obtaining the model fitting p-value
# from a sem object:
ModelFittingPvalue=function(semObject)
{
  ModelChisq=summary(semObject)$chisq
  df= summary(semObject)$df
  ModelFittingPvalue=1-pchisq(ModelChisq,df) 
  ModelFittingPvalue
}
# this function calculates the single anchor score
LEO.SingleAnchor=function(M1,A,B)
{
  datObsVariables= data.frame(M1,A,B)
  # this is the observed correlation matrix
  S.SingleAnchor=cor(datObsVariables,use="p")
  m=dim(datObsVariables)[[1]]
  semModel1 =sem(CausalModel1,S=S.SingleAnchor,N=m)
  semModel2 =sem(CausalModel2,S=S.SingleAnchor,N=m)
  semModel3 =sem(CausalModel3,S=S.SingleAnchor,N=m)
  semModel4 =sem(CausalModel4,S=S.SingleAnchor,N=m)
  semModel5 =sem(CausalModel5,S=S.SingleAnchor,N=m)
  # Model fitting p values for each model
  p1=ModelFittingPvalue(semModel1)
  p2=ModelFittingPvalue(semModel2)
  p3=ModelFittingPvalue(semModel3)
  p4=ModelFittingPvalue(semModel4)
  p5=ModelFittingPvalue(semModel5)
  MaxP<- max(p2,p3,p4,p5)
  #if(MaxP==0)
  #{
  #  MaxP<- 0.0000000001
  #}
  LEO.SingleAnchor=log10(p1/MaxP)
  data.frame(LEO.SingleAnchor,p1,p2,p3,p4,p5) 
  
} # end of function
#

#############################

# Prepare output dataframe 
output_df <- data.frame(SNP = character(0),
                        ME = character(0),
                        Trait = character(0),
                        LEO.score = numeric(0),
                        Cor_ME_SNP = numeric(0),
                        Cor_ME_Trait = numeric(0),
                        Cor_SNP_Trait = numeric(0),
                        stringsAsFactors = FALSE)

# Loop through SNPs, MEs, and Traits
for (i in 1:SNP_Num) {
  SNP <- as.numeric(datSNP[, i])
  SNPName <- colnames(datSNP)[i]
  
  # Skip SNPs with no variation
  if (all(SNP == 0) | all(SNP == 1)) {
    next
  }
  
  for (j in 1:MEs_Num) {
    ME <- MEs[, j]
    MEName <- colnames(MEs)[j]
    
    for (k in 1:Traits_Num) {
      Trait <- datTraits[, k]
      TraitName <- colnames(datTraits)[k]
      
      # Skip traits with no variation
      if (all(Trait == 0)) {
        next
      }
      
      # Calculate correlations
      CorInteg <- cor(cbind(SNP, ME, Trait), use = "p")
      Cor_ME_SNP <- CorInteg[2, 1]
      Cor_ME_Trait <- CorInteg[2, 3]
      Cor_SNP_Trait <- CorInteg[1, 3]
      
      # Calculate LEO score
      Anchor <- LEO.SingleAnchor(M1 = SNP, A = ME, B = Trait)[[1]]
      
      # If Anchor is positive, store the result
      if (Anchor > 0) {
        output_df <- rbind(output_df, data.frame(
          SNP = SNPName, 
          ME = MEName, 
          Trait = TraitName, 
          LEO.score = Anchor, 
          Cor_ME_SNP = Cor_ME_SNP, 
          Cor_ME_Trait = Cor_ME_Trait, 
          Cor_SNP_Trait = Cor_SNP_Trait
        ))
      }
    }
  }
}

# View the output
head(output_df)

write.csv(output_df,"LEO_SNP_ME_Trait.csv", row.names = FALSE)

