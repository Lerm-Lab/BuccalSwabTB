## Author details
## Script name: Script_MouthSwabStudy_pipeline_v1.R
## Purpose of the script: To analyze the samples using minfi and ChAMP
##
## Author(s): Shumaila Sayyab, Ph.D.
## Date Created: 2022-03-15
## Date Modified: 2022-03-30
## Copyright statement: MIT open license
## Contact information: shumaila.sayyab@liu.se
## Please cite: @article{,
##  title={MouthSwab_Study},
##  author={},
##  journal={},
##  year={2022},
##  publisher={Cold Spring Harbor Laboratory Press}}
## Notes: The user needs to have the IDAT files.
##
## github link:


####################################################################

## Pipeline for DNA methylation analysis

# 1. Processing of IDAT files with minfi
# 2. Preparing the data and Pheno files and Annotations
# 3. Minfi QC step: Normalization, quality check, filtering (poor Quality probes/samples, chr XY, probes over SNPs, cross reactive probes)
# 4. SVD analysis using ChAMP
# 5. Batch affect correction using Combat
# 6. Density plots, MDS plots, PCA plots to visualize the data (before and after QC/correction)
# 7. Co-variates / confounder identification
# 8. Cell type adjustments and analysis
# 9. DMC/DMGs identification
#10. Volcano, Heatmaps and Pathway Enrichment (ORA)/Network analysis
#11.Estimated genomic inflation and bais using BACON package in R.


##
getwd()
setwd("/Users/lovka39/Documents/Linkoping2021/Results/DNAmethylation/")
pack_R <- c("minfi", "gmqn", "ChAMP", "bnstruct","missMethyl","limma","RColorBrewer","stringr","ggpubr","ggrepel","EnhancedVolcano","calibrate","FlowSorted.Blood.EPIC")
for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
  
  col1=c("Control"="#648FFF","Exposed"="#15958a","Patient"="#DC267F")
}

#####################################################################
# 1. Processing of IDAT files with minfi                           #
#####################################################################

#--------------------------- loading the rgSet object or read in the raw dataset
# set the data directory
dataDir <- setwd(getwd())
list.files(dataDir)
baseDir<-"*_Maria_Lerm"
targets<-read.metharray.sheet(baseDir,pattern="Sample_sheet.csv")
#read in the raw data from the IDAT files; warnings can be ignored.
rgSet <- read.metharray.exp(targets = targets)
sampleNames(rgSet) <- targets$Sample_Name
getManifest(rgSet)
colnames(rgSet)

############################################################
## 2. Preparing the data and Pheno files, Annotations Hg38
#############################################################

## subset dataset
dim(rgSet)
sampleNames(rgSet)
pData(rgSet)

MET_pd <- read.csv("../DNAmethylation/sample_Sheet_validation.csv",sep=";")
MET_pd$SampleName
sampleNames(rgSet)<-MET_pd$SampleName
MET_pd$SampleName<-sampleNames(rgSet)
rgSet$Sample_Name<-sampleNames(rgSet)


## ################################

#-------------------------- Prepare the Pheno file

match(colnames(rgSet),MET_pd$SampleName)
rgSet$Age<-MET_pd$Age
rgSet$Sex<-MET_pd$Sex
rgSet$BMI<-MET_pd$BMI
rgSet$IGRA<-MET_pd$IGRA.status

samples_info<-MET_pd
samples_info$Slide<-rgSet$Slide
samples_info$Array<-rgSet$Array
rgSet$Sample_Group<-samples_info$Group
samples_info$Sample_Group<-samples_info$Group

save(rgSet, file = "rgSet_MS_validation.RData")
save(samples_info,file="samples_info_MS_validation.RData")

 
write.table(samples_info,"samples_info_validation.csv",sep=",",quote = F,row.names = F)
head(pData(rgSet))
names(samples_info)
samples_info
samples<-samples_info
## ##############


library(minfi)
library(devtools)

devtools::install_github("achilleasNP/IlluminaHumanMethylationEPICmanifest") 
devtools::install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")
annotate_hg38<-function(rgSet)
   library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
  rgSet@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")
  # get the 850k annotation data
  ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
  
  
############################################################
## 3. Pipeline Minfi QC steps
#############################################################
minfi_pipeline<-function(rgSet,samples)
{
{ #-------------- Quality Control of the dataset

  detP <- detectionP(rgSet)
  qcReport(rgSet, sampNames=samples$Sample_Name, sampGroups=samples$Sample_Group,
           pdf="qcReport.pdf")
  keep <- colMeans(detP) < 0.05

  rgSet <- rgSet[,keep]
  table(keep)
  
  ## remove poor quality samples from targets data
  samples <- samples[keep,]
  dim(samples)
  samples[,1:5]

  ## remove poor quality samples from detection p-value table
  detP <- detP[,keep]

  # Normalization
  mSetSq <- preprocessSWAN(rgSet) 
  dim(mSetSq)
  
  # Visualise what the data looks like before and after normalisation

  ##----------------------Plot density

  NormPlot <- function(normPlot){
    par(mfrow=c(1,2))
    
    densityPlot(rgSet, sampGroups=samples$Sample_Group,main="Raw", legend=FALSE,pal = col1)
    legend("top", legend = levels(factor(samples$Sample_Group)),
           text.col=col1)
    densityPlot(getBeta(mSetSq), sampGroups=samples$Sample_Group,
                main="SWAN.Normalized", legend=FALSE,pal=col1)
    legend("top", legend = levels(factor(samples$Sample_Group)),
           text.col=col1)
  }
  name<-paste0(mycell,"_Raw_Density.png")
  png(name,height = 10, width = 10, units = "in", res = 300)
    NormPlot()
  dev.off()

  ##---------------------Plot MDS customized function

  png("MDS_SWAN.png",height = 10, width = 10, units = "in", res = 300)
  mdsPlot(dat =getBeta(mSetSq) , numPositions = 1000, sampNames = samples$SampleName, sampGroups = samples$Sample_Group,
          pch = 19, pal = col1,
          legendPos ="topright",
          main = "Normalized beta matrix b/f QC")
  dev.off()
  
  }
  
  ##-------------------- Filtering of the data
  detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

  ## remove any probes that have failed in one or more samples
  keep <- rowSums(detP < 0.01) == ncol(mSetSq)
  table(keep)

   mSetSqFlt <- mSetSq[keep,]
  dim(mSetSqFlt)

  ## Remove probes on the sex chromosomes
  keep <- !(featureNames(mSetSqFlt) %in% ann850k$Name[ann850k$chr %in%
                                                        c("chrX","chrY")])

  mSetSqFlt <- mSetSqFlt[keep,]
  dim(mSetSqFlt)
  
  ## remove probes with SNPs at CpG site
  gset <- mapToGenome(mSetSqFlt)
  mSetSwFlt <- dropLociWithSnps(gset)
  dim(mSetSwFlt)

  ## exclude cross reactive probes
 
  library(RCurl)
  nsp <- read.csv("NSP_48639-non-specific-probes-Illumina450k 2.csv", sep = ";")
  colnames(nsp) 
  targets <- nsp$TargetID
  keep <- !(featureNames(mSetSwFlt) %in% nsp$TargetID)
  table(keep)
  mSetSwFlt <- mSetSwFlt[keep,]
  dim(mSetSwFlt)
 

  # Calculation of M values and beta values
  ## calculate M-values for statistical analysis
 
   mVals <- getM(mSetSqFlt)
  bVals <- getBeta(mSetSqFlt)
  
  colnames(bVals)<-samples$SampleName
  match(colnames(bVals),samples$SampleName)

  ## MDS cleaned 
  png("MDS_bVals.png",height = 10, width = 10, units = "in", res = 300)
  mdsPlot(dat =bVals , numPositions = 1000, sampNames = NULL, sampGroups = samples$Sample_Group,
             pch = 19, pal = col1,
             legendPos ="topright",
             main = "Normalized beta matrix After QC")
  
##---------- save the R objects
  save(bVals,file="MS_bVals.RData")
  save(mVals,file="MS_mVals.RData")
  

} ##end of minfi

############################################################
## 4. SVD analysis using ChAMP
## 5. Batch affect correction using Combat
## 6. Confounders and cov. identification
#############################################################

SVDplot.champ<-function(bVals,samples)

  library(ChAMP)
  
  ##-------------------- SVD and batch effects

  tmp_samples<-samples[,c(7,10,3:6,8,9)]
  colnames(tmp_samples)
  match(samples$SampleName,colnames(bVals))

  bVals_f <- as.data.frame(bVals)
  champ.SVD(beta=bVals_f,pd=tmp_samples)
  
  
  bVals.corrected<-champ.runCombat(beta=bVals,pd=tmp_samples)
  bVals_corrected_f <- as.data.frame(bVals.corrected)
  champ.SVD(beta=bVals_corrected_f,pd=tmp_samples)

  save(bVals.corrected,file="MS_bVals.corrected.RData")

## plot MDS of the batch corrected beta values with shape indicator for the group
  shapes= c(19, 17, 15)
  shapes <- shapes[samples$Group]
  
  bVals.corrected <- as.matrix(bVals.corrected)
  png("MDS_bVals.corrected.png",height = 10, width = 10, units = "in", res = 300)
  mdsPlot(dat =bVals.corrected, numPositions = 1000, sampNames = NULL, sampGroups = samples$Sample_Group,
          pch = shapes, pal = col1, xlim = c(-10,10), ylim= c(-10,10),
          legendPos ="topright",
          main = "MDS of Norm. corrected beta matrix")
  dev.off()
 
 ## end of batch correction and SVD

#####################################################################################################
## 6. Density plots, MDS plots, PCA plots to visualize the data (before and after QC/correction) ##
######################################################################################################

##----------------------  PCAPlot.merge
PCAPlot.merge <-function(bVals.corrected,tmp_samples,p1,p2)
install.packages("FactoMineR")
install.packages("pander")
  
  
  library(FactoMineR)
  library(factoextra)
  library(gridExtra)
  library(devtools)
  library(jcolors)
  library(ggpubr)
  library(ggplot2)
  library(pander)
  library(gridExtra)
library(dplyr)

#################################################################################
##PCA 
################################################################################
sample<- samples_info

data1_top <- as.matrix(bVals.corrected)
dim(bVals.corrected)
betaPCAt<-as.data.frame(t(data1_top))
dim(betaPCAt)
dim(betaPCAt)
betaPCAfig <- betaPCAt[,1:841687]

PCA <- PCA(betaPCAfig, graph = F,ncp=10) # beta PCA

png("PCA_bVals.correctd.png", height = 9, width = 10, units = "in", res = 300)

p <- fviz_pca_ind(PCA, repel = T, 
                  geom.ind =  "point", #c("point","text"),
                  pch = shapes,
                  pointsize = 5, 
                  col.ind = betaPCAt$sample, 
                  mean.point = FALSE,
                  palette = c("Control"="#648FFF","Exposed"="#15958a","Patient"="#DC267F"),
                  title = "PCA", legend.title = list(col = "Group"),
                  addEllipses = T, ellipse.level = 0.85) + xlim(-1500,2000)+ ylim(-1200, 1000)


p + theme(text = element_text(family = "Arial", size = 34), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(size = 1.5))

dev.off()

###############################################################################

# 8. Cell type adjustments and analysis
############### 2. EPIDISH : est. cell prop. for Epi, Fib, Fat and IC.
## used the robust partial correction method

###############################################################################

Epidish_analysis<-function()
{
  library(EpiDISH)
  library(gridExtra)
  data(centEpiFibIC.m)
  
  frac.m <- hepidish(beta.m = bVals, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m[,c(1,2,3,4,5,6,7)], h.CT.idx = 3, method = 'RPC')
  boxplot(frac.m)
  
  png("Epidish_MS.png",height = 10, width = 10, units = "in", res = 300 )
  boxplot(frac.m)
  dev.off()
  
  write.table(frac.m, "Epidish_MS_Frac.m_incl_7_blood_celltypes.csv",
              sep = ",",
              quote = F,
              row.names = T)
  
}

######################################################################

## 9.DMC idenfication using linear model (Limma package in R)

######################################################################

DMC_identifcation<-function()
{
  library(limma)
  
  library(minfi)
  
  targets<-samples
  match(colnames(bVals.corrected),samples$SampleName)
  group <- factor(targets$Sample_Group)

  design_2<-model.matrix(~0+group,data=targets)
  colnames(design_2) <- c(levels(group))
  contMatrix <- makeContrasts(Patient-Control, Exposed-Control, Patient-Exposed,levels=design_2)
  dim(design_2)
  fit<-lmFit(bVals.corrected,design_2)
  colnames(fit$coefficients)
  rownames(contMatrix)
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)
  summary(decideTests(fit2))

  ##--- Annotate the DMCs

  ann850kSub_combat <- ann850k[match(rownames(bVals.corrected),ann850k$Name),
                               c(1:ncol(ann850k))]

  Patient_Control<-topTable(fit2,coef = 1,number = Inf,genelist = ann850kSub_combat)
  Exposed_Control<-topTable(fit2,coef = 2,number = Inf,genelist = ann850kSub_combat)
  Patient_Exposed<-topTable(fit2,coef = 3,number = Inf,genelist = ann850kSub_combat)

  getwd()
  save(Patient_Control,file="Patient_Control.MouthSwab.validation.RData")
  save(Exposed_Control,file="Exposed_Control.MouthSwab.validation.RData")
  save(Patient_Exposed,file="Patient_Exposed.MouthSwab.validation.RData")
  #

  #################

  topset1<-subset(Patient_Control, Patient_Control$adj.P.Val<0.05 & abs(Patient_Control$logFC)>0.2)
  topset2<-subset(Exposed_Control, Exposed_Control$adj.P.Val<0.05 & abs(Exposed_Control$logFC)>0.2)
  topset3<-subset(Patient_Exposed, Patient_Exposed$adj.P.Val<0.05 & abs(Patient_Exposed$logFC)>0.2)
  
  top1_top2<-merge(topset1, topset2, by="Name",all=TRUE, sort = TRUE,no.dups = TRUE)
  top1_top2_top3<-merge(top1_top2, topset3, by="Name",all=TRUE, sort = TRUE,no.dups = TRUE)
  dim(top1_top2_top3)
  save(top1_top2_top3,file="top1_top2_top3.MouthSwab.RData")
  head(top1)
  #


  write.table(topset1,"DMC_Patient_Control_padj0.05_0.2.csv",sep = ",",quote=F, row.names = F)
  write.table(topset2,"DMC_Exposed_Control_padj0.05_0.2.csv",sep = ",",quote=F, row.names = F)
  write.table(topset3,"DMC_Patient_Exposed_padj0.05_0.2.csv",sep = ",",quote=F, row.names = F)
  write.table(top1_top2_top3,"DMC_MERGED_Top1_Top2_Top3_padj0.05_0.2.csv",sep = ",",quote=F, row.names = F)
  #
  betaPCA <- bVals.corrected
  tp1 <- top1_top2_top3
  samples <- samples_info


}

save.image(file = "MS.RData")


######################################################################
######################################################################


##--------------- end of pipeline


