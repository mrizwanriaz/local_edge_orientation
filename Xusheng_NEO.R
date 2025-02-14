
########################################################################
###  12.3   Systems Genetic Analysis with NEO
########################################################################

##  R Code for the LEO.SingleAnchor Score
library(sem)


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
                          M1 <-> M1, sigmaM1M1,NA"
)


# Single anchor model 4: M1->A<-B
CausalModel4=specifyModel(text="
                          M1 -> A, gammaM1toA, NA
                          B -> A, gammaBtoA, NA
                          A <-> A, sigmaAA, NA
                          B <-> B, sigmaBB, NA
                          M1 <-> M1, sigmaM1M1,NA"
)


# Single anchor model 5: M1->B<-A
CausalModel5=specifyModel(text="
                          M1 -> B, gammaM1toB, NA
                          A -> B, gammaAtoB, NA
                          A <-> A, sigmaAA, NA
                          B <-> B, sigmaBB, NA
                          M1 <-> M1, sigmaM1M1,NA
                          ")


# Function for obtaining the model fitting p-value
# from an sem object:
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

#LEO.SingleAnchor(M1=x1,A=y1,B=y2)


# Get SNP markers which encode the mQTLs of the blue module
###  the file BXD.geno -> BXD.geno_v2.txt:
###  1) delete lines which begin with @; 2) values B -> 1, D ->-1, H/U -> 0.

SNP_origin=read.table("./BXDphenotype/BXD.geno_v2.txt", sep="\t",header = TRUE)
SNPdataFemale=SNP_origin[,-c(1,3,4)]
row.names(SNPdataFemale)<- SNPdataFemale$Locus
SNPdataFemale<- SNPdataFemale[,-1]
SNPdataFemale<- as.data.frame(t(SNPdataFemale))
write.csv(SNPdataFemale,"./Integ/SNPdataFemale.csv")
# find matching row numbers to line up the SNP data
SNPdataFemale_integ<- read.csv("./Integ/SNPdataFemale_integ.csv", sep=",", header = TRUE)
row.names(SNPdataFemale_integ)<- SNPdataFemale_integ[,1]
SNPdataFemale_integ<- SNPdataFemale_integ[,-1]

###  samples C57BL.6J and DBA.2J in datExprFemale don't exist in SNPdataFemale, these two samples would be deleted in advance.
datExprFemale_v2<- datExprFemale[-c(1,2,36,37,73,100,106,107),]
row.names(datExprFemale_v2)
snpRows=match(dimnames(datExprFemale_v2)[[1]],row.names(SNPdataFemale_integ))
# define a data frame whose rows correspond to those of datExpr
datSNP = SNPdataFemale_integ[snpRows,]
ChooseRows<- row.names(SNPdataFemale_integ)[snpRows]
rownames(datSNP)=ChooseRows
# show that row names agree
table(rownames(datSNP)==rownames(datExprFemale_v2))



Exam_datSNP<- datSNP
#Exam_datSNP<- datSNP[,1:10]
Exam_MEsFemale<- MEsFemale_delgrey[-c(1,2,36,37,73,100,106,107),]
Exam_datTraits<- datTraits[-c(1,2,36,37,73,100,106,107),]
colnames(Exam_datTraits)<- gsub("X", "", colnames(Exam_datTraits))
ExpSampNum<- length(row.names(datExprFemale))
SampFilt<- apply(Exam_datTraits, 2, function(x) (length(x[!is.na(x)])/ExpSampNum)>0.5)
PhenoFiltered<- Exam_datTraits[,SampFilt]
write.csv(PhenoFiltered, "./Integ/PhenoFiltered.csv")
CorTest<- cor(PhenoFiltered$"12450", PhenoFiltered, use = "p")
#write.csv(CorTest, "./Integ/CorTest.csv")
CorTest<- t(CorTest)
names(CorTest)<- colnames(PhenoFiltered)
CorTest<- CorTest[!is.na(CorTest)]
CorTest<- CorTest[CorTest>0.5]
write.csv(CorTest, "./Integ/PhenoID_12450_cor_PhenoFilterd_05.csv")
Exam_datTraits<- Exam_datTraits[,(colnames(Exam_datTraits) %in% names(CorTest))]
write.csv(Exam_datTraits, "./Integ/Exam_datTraits.csv")

SNP_Num<- dim(Exam_datSNP)[[2]]
Traits_Num<- dim(Exam_datTraits)[[2]]
MEs_Num<- dim(Exam_MEsFemale)[[2]]

out_title<- "SNP,ME,Trait,Anchor"


write.table(out_title, "./Integ/LEO_SingleAnchor.csv", col.names = FALSE, row.names = FALSE)

for(i in 1:SNP_Num)
{
  SNP<- as.numeric(Exam_datSNP[,i])
  SNPName<- names(Exam_datSNP)[i]
  for(j in 1:MEs_Num)
  {
    ME<- Exam_MEsFemale[,j]
    MEName<- names(Exam_MEsFemale)[j]
    for(k in 1:Traits_Num)
    {
      Trait<- Exam_datTraits[,k]
      TraitName<- names(Exam_datTraits)[k]
      if(all(Trait==0))
      {
        next
      }
      Anchor<- LEO.SingleAnchor(M1=SNP,A=ME,B=Trait)[[1]]   ### if Anchor is negative, there is no evidence for a causal relationship between the module and the trait when the SNP is used as a causal anchor.
      Anchor
      if(Anchor>0)
      {
        add<- paste(SNPName, MEName, TraitName, Anchor, sep=",")
        write.table(add, "./Integ/LEO_SingleAnchor.csv", sep="", append = TRUE, col.names = FALSE, row.names = FALSE)
      }
    }
  }
  
}

#write.table(out, "./Integ/LEO_SingleAnchor.csv", sep="", col.names = FALSE, row.names = FALSE)


LEO_SingleAnchor<- read.csv("./Integ/LEO_SingleAnchor.csv", sep=",", header = TRUE)
LEO_Sele<- aggregate(LEO_SingleAnchor$Anchor, by=list(SNP=LEO_SingleAnchor$SNP, ME=LEO_SingleAnchor$ME), FUN=max)
colnames(LEO_Sele)<- c("SNP", "ME", "Anchor")
write.csv(LEO_Sele, "./Integ/LEO_Sele_first.csv", row.names = FALSE)
MatchID<- match(LEO_SingleAnchor$Anchor, LEO_Sele[,3])
LEO_Sele_final<- LEO_SingleAnchor[!is.na(MatchID),]
write.csv(LEO_Sele_final, "./Integ/Integ_LEO_Sele.csv", row.names = FALSE)
#write.csv(base::unique(LEO_Sele_final[!is.na(LEO_Sele_final[,1]),]), "./Integ/Integ_LEO_Sele.csv", row.names = FALSE)

setwd("F:/WXSheng/Issue2_V3")
Integ_LEO_Sele_v2<- read.csv("./Integ/Integ_LEO_Sele.csv", sep=",", header = TRUE)
Integ_LEO_Sele_v3<- Integ_LEO_Sele_v2[(Integ_LEO_Sele_v2$Anchor>=0.8),]
write.csv(Integ_LEO_Sele_v3, "./Integ/Integ_LEO_Sele_08.csv", row.names = FALSE)
Integ_LEO_Sele_v4<- Integ_LEO_Sele_v2[(Integ_LEO_Sele_v2$Anchor>=0.85),]
write.csv(Integ_LEO_Sele_v4, "./Integ/Integ_LEO_Sele_085.csv", row.names = FALSE)
Integ_LEO_Sele_v5<- Integ_LEO_Sele_v2[(Integ_LEO_Sele_v2$Anchor>=0.9),]
write.csv(Integ_LEO_Sele_v5, "./Integ/Integ_LEO_Sele_09.csv", row.names = FALSE)
Integ_LEO_Sele_v6<- Integ_LEO_Sele_v2[(Integ_LEO_Sele_v2$Anchor>=0.95),]
write.csv(Integ_LEO_Sele_v6, "./Integ/Integ_LEO_Sele_095.csv", row.names = FALSE)
Integ_LEO_Sele_v7<- Integ_LEO_Sele_v2[(Integ_LEO_Sele_v2$Anchor>=1),]
write.csv(Integ_LEO_Sele_v7, "./Integ/Integ_LEO_Sele_1.csv", row.names = FALSE)



### Cytoscape In
#CytoIn<- Integ_LEO_Sele_v3
CytoIn<- Integ_LEO_Sele_v4
#Nodes<- cbind(CytoIn$SNP, base::unique(CytoIn$ME), as.character(base::unique(CytoIn$Trait)))
#Nodes<- paste(CytoIn$SNP, base::unique(CytoIn$ME), as.character(base::unique(CytoIn$Trait)))
SNP_ME<- data.frame(paste(CytoIn[,1], "=",round(CytoIn[,4], 2), sep=""), as.character(rep(1, dim(CytoIn)[1])), CytoIn[,2])
colnames(SNP_ME)<- c("Source", "Type", "Target")
ME_Pheno<- data.frame(CytoIn[,2], as.character(rep(2, dim(CytoIn)[1])), as.character(CytoIn[,3]))
colnames(ME_Pheno)<- c("Source", "Type","Target")
SNP_ME_Pheno<- rbind(SNP_ME, ME_Pheno)
write.table("Nodes\tType", "./Integ/CytoIn_Nodes.txt", sep="\t", row.names = FALSE, col.names = FALSE)
write.table(paste(CytoIn$SNP, as.character(rep(1, length(CytoIn$ME))), sep="\t"), "./Integ/CytoIn_Nodes.txt", sep="\t", row.names = FALSE, col.names = FALSE, append = TRUE)
write.table(paste(base::unique(CytoIn$ME),as.character(rep(2, length(base::unique(CytoIn$ME)))), sep="\t"), "./Integ/CytoIn_Nodes.txt", sep="\t", row.names = FALSE, col.names = FALSE, append = TRUE)
write.table(paste(as.character(base::unique(CytoIn$Trait)), as.character(rep(3, length(base::unique(CytoIn$Trait)))), sep="\t"), "./Integ/CytoIn_Nodes.txt", sep="\t", row.names = FALSE, col.names = FALSE, append = TRUE)
write.table(SNP_ME_Pheno, "./Integ/CytoIn_Net.txt", sep="\t", row.names = FALSE)



SNP=as.numeric(datSNP$rs3677817)
weight=as.numeric(datTraits$X10005)
weight_v2<- weight[-c(1,2,36,37,73,100,106,107)]
MEblue=MEsFemale$MEblue
MEblue_v2<- MEblue[-c(1,2,36,37,73,100,106,107)]

# evaluate the relative causal fit
# of SNP -> MEblue -> weight

LEO.SingleAnchor(M1=SNP,A=MEblue_v2,B=weight_v2)[[1]]
# output -2.823433
# The negative LEO score indicates that there is no evidence for a causal relationship MEblue_v2 -> weight_v2, if SNP (as.numeric(datSNP$rs3677817)) is used as causal anchor.



###  identify genes inside the blue module that have a causal effect on weight


SNP_Net<- datSNP[,which(colnames(datSNP) %in% CytoIn$SNP)]  ### dim(SNP_Net)==[35,25]
ME_Names<- base::unique(CytoIn$ME)
ME_Net<- MEsFemale[-c(1,2),which(colnames(MEsFemale) %in% ME_Names)]   ### dim(ME_Net)==[35,2]
Trait_Names<- base::unique(CytoIn$Trait)
Trait_Net<- datTraits[-c(1,2),which(gsub("X", "", colnames(datTraits)) %in% Trait_Names)]   ### dim(Trait_Net)==[35,3]
Module_Colors<-gsub("ME", "", ME_Names) 
write.table("SNP,Gene,Module,Trait,Cor_Gene_SNP,Cor_Gene_Trait,Stress,Alcohol", "./Integ/Cor_Gene_SNP_Trait.csv", col.names = FALSE, row.names = FALSE)

Gene_on_alcohol_stress <- read.table("./Integ/Gene_on_alcohol_stress.txt", sep = "\t", header = TRUE)
Gene_on_stress<- Gene_on_alcohol_stress[which(Gene_on_alcohol_stress[,3] %in% "Y"),2]
Gene_on_alcohol<- Gene_on_alcohol_stress[which(Gene_on_alcohol_stress[,4] %in% "Y"),2]


for(i in 1:dim(SNP_Net)[2])
{
  SNP_Name<- colnames(SNP_Net)[i]
  for(k in 1:dim(Trait_Net)[2])
  {
    Trait_Name<- gsub("X", "", colnames(Trait_Net))[k]
    for(m in 1:length(Module_Colors))
    {
      whichmodule<- Module_Colors[m]
      Module_Name<- whichmodule
      restGenes=moduleColorsFemale==whichmodule
      datExprModule=datExprFemale_v2[,restGenes]
      for(n in 1:dim(datExprModule)[2])
      {
        Gene_Name<- colnames(datExprModule)[n]
        CorInteg<- signif(cor(data.frame(SNP_Net[,i],datExprModule[,n],Trait_Net[,k]),use="p"),2)
        Cor_Gene_SNP<- CorInteg[2,1]
        Cor_Gene_Trait<- CorInteg[2,3]
        if((abs(Cor_Gene_SNP)>=0.5) && (abs(Cor_Gene_Trait)>=0.5))
        {
          if(Gene_Name %in% Gene_on_stress)
          {
            Judge_stress<- "Y"
          }
          else
          {
            Judge_stress<- ""
          }
          if(Gene_Name %in% Gene_on_alcohol)
          {
            Judge_alcohol<- "Y"
          }
          else
          {
            Judge_alcohol<- ""
          }            
          write.table(paste(SNP_Name, Gene_Name, Module_Name, Trait_Name, Cor_Gene_SNP, Cor_Gene_Trait,Judge_stress, Judge_alcohol, sep=","), "./Integ/Cor_Gene_SNP_Trait.csv", sep="", col.names = FALSE, row.names = FALSE, append = TRUE)   
        }
        
      }
    }
  }
}



#whichmodule=Module_Colors
whichmodule<- "blue"
restGenes=moduleColorsFemale==whichmodule   ###  moduleColorsFemale equals moduleColorsAutomatic, and has one column which is comprised of colors.The restGenes is comprised of logical values, the genes in the blue module have been set TRUE.
datExprModule=datExprFemale_v2[,restGenes]   ### obtain the exprssion values of genes in the blue module
#write.csv(datExprModule,"./datExprModule.csv")
attach(datExprModule)###   the function attach is to call objects in the datExprModule directly, not through $.e.g. datExprModule$MMT00000549 also can be called using MMT00000549
CorInteg<- signif(cor(data.frame(SNP,MEblue_v2,datExprModule$`101023`,weight_v2),use="p"),2)
#CorInteg<- signif(cor(data.frame(SNP,datExprModule$`101023`,weight_v2),use="p"),2)
Cor_Gene_SNP<- CorInteg[3,1]
Cor_Gene_Trait<- CorInteg[3,4]


LEO.SingleAnchorWeight=rep(NA, dim(datExprModule)[[2]] )    ### In datExprModule, rows are sample names and columns are gene names.
for (i in c(1:dim(datExprModule)[[2]]) )
{
  #printFlush(i)    ### Returns the value of the print function.

  LEO.SingleAnchorWeight[i]=LEO.SingleAnchor(M1=SNP,A=datExprModule[,i],B=weight_v2)[[1]]
}

# Find gens with a very significant LEO score
restCausal=LEO.SingleAnchorWeight>0.1 # could be lowered to 1
names(datExprModule)[restCausal]

# this provides us the corresponding gene symbols
data.frame(datOutput[restGenes,c(1,6)])[restCausal,]
#output

