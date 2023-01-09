source("3-IPM_modeling/Scripts/IPM_functions.R") # load IPM functions for this script!
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(dplyr)
#library(scales)
#library(tidyr)
#library(reshape2)
library(stringr)
#library(ggfortify)
library(lubridate)

################################################################################
############# Functions #####################
################################################################################

RunIPMfromParamMatrix=function(ModPrm,RecValCol="MN.RecVal_Sec",IPMstructure){
  modlist= vector(mode = "list", length = nrow(ModPrm))
  for(i in which(!is.na(ModPrm[,RecValCol]))){
    MP=list(n = IPMstructure$n,
            ds = IPMstructure$ds,
            y = IPMstructure$y,
            rec.size = IPMstructure$rec.size,
            
            rec = ModPrm[,RecValCol][i],
            g.int = ModPrm$g.int[i], 
            g.slp = ModPrm$g.slp[i], 
            g.var = ModPrm$g.var[i], 
            s.int = ModPrm$s.int[i], 
            s.slp = ModPrm$s.slp[i],
            Interval_Years=ModPrm$Interval_Years[i])
    modlist[[i]]=bigmatrix(MP)
  }
  return(modlist)
}

Region_RunIPMfromParamMatrix=function(ModPrm,RecValCol="SiteStock_mean",IPMstructure){
  modlist= vector(mode = "list", length = nrow(ModPrm))
  for(i in which(!is.na(ModPrm[,RecValCol]))){
    MP=list(n = IPMstructure$n,
            ds = IPMstructure$ds,
            y = IPMstructure$y,
            rec.size = IPMstructure$rec.size,
            
            rec = ModPrm[,RecValCol][i],
            g.int = ModPrm$g.int[i], 
            g.slp = ModPrm$g.slp[i], 
            g.var = ModPrm$g.var[i], 
            s.int = ModPrm$s.int[i], 
            s.slp = ModPrm$s.slp[i])
    modlist[[i]]=bigmatrix_regional(MP)
  }
  return(modlist)
}

lamFromModList=function(ModList){
  N=length(ModList)
  lamvec=rep(NA,N)
  i_s=which(!unlist(lapply(ModList,is.null)))
  for(i in i_s){
    lamvec[i]=ModList[[i]]$lam
  }
  return(lamvec)
}


XFromModList=function(ModList,outparam){
  N=length(ModList)
  outvec=rep(NA,N)
  i_s=which(!unlist(lapply(ModList,is.null)))
  for(i in i_s){
    outvec[i]=(ModList[[i]][outparam])
  }
  return(outvec)
}

PlotKnorm=function(Kmat,Nr=4,ImageTitle="Elasticity",NORM_REC=TRUE,colmax=1,sizeclasses=y){ #colmax = 1 is orig
  M=dim(Kmat)[1]
  Kn=Kmat
  if(NORM_REC){Kn[1:Nr,(Nr+1):M]=Kmat[1:Nr,(Nr+1):M]/max(Kmat[1:Nr,(Nr+1):M])}
  #Kn[(Nr+1):M,]=Kmat[(Nr+1):M,]/max(Kmat[(Nr+1):M,])
  image(x=sizeclasses,y=sizeclasses,z=t(Kn),
        xlab="Log Size (time = T)",ylab="Log Size (time = T+1)",
        cex.lab=1.5,
        main=ImageTitle,zlim=c(0,colmax))
  abline(a=0,b=1)
  abline(h=(Nr-.5)/M)
  abline(v=(Nr-.5)/M)
}

FullModelReport=function(ModList,ModDF,ModNum,lambda_str="lambda_SiteStock",lA=y){
  par(mfrow=c(4,3))
  A=10^lA
  D=2*sqrt(A/pi)
  plot(lA,XFromModList(ModList,outparam = "w")[[ModNum]],xlab="Log Size (time=t)",ylab="Probability Density",
       main=paste0("Stable Size Distribution (w):\n",ModDF$SIG[ModNum],"\nL= ",
                   round(ModDF[,lambda_str][ModNum],3)),type="b",cex=.1)
  plot(lA,XFromModList(ModList,outparam = "v")[[ModNum]],xlab="Log Size (time=t)",ylab="Probability Density",
       main=paste0("Reproductive Value (v):\n",ModDF$SIG[ModNum],"\nL= ",
                   round(ModDF[,lambda_str][ModNum],3)),type="b",cex=.1)
  plot.new()
  text(0.5,0.5,paste0("Lambda: ",round(XFromModList(ModList,outparam = "lam")[[ModNum]],3),
                      "\nv.dot.w: ",round(XFromModList(ModList,outparam = "v.dot.w")[[ModNum]],3)),cex=1.5)
  PlotKnorm(XFromModList(ModList,outparam = "Gk")[[ModNum]],
            ImageTitle=paste0("G: ",ModDF$SIG[ModNum],"\nL= ",
                              round(ModDF[,lambda_str][ModNum],3)),NORM_REC = F)
  plot(lA,XFromModList(ModList,outparam = "Sk")[[ModNum]],xlab="Log Size (time=t)",ylab="Survival Probability",
       main=paste0("S: ",ModDF$SIG[ModNum],"\nL= ",
                   round(ModDF[,lambda_str][ModNum],3)),type="b",cex=.1)
  PlotKnorm(XFromModList(ModList,outparam = "Pk")[[ModNum]],
            ImageTitle=paste0("P: ",ModDF$SIG[ModNum],"\nL= ",
                              round(ModDF[,lambda_str][ModNum],3)),NORM_REC = F)
  PlotKnorm(XFromModList(ModList,outparam = "K")[[ModNum]],
            ImageTitle=paste0("K: ",ModDF$SIG[ModNum],"\nL= ",
                              round(ModDF[,lambda_str][ModNum],3)),NORM_REC = F)
  PlotKnorm(XFromModList(ModList,outparam = "sens")[[ModNum]],
            ImageTitle=paste0("sens: ",ModDF$SIG[ModNum],"\nL= ",
                              round(ModDF[,lambda_str][ModNum],3)),NORM_REC = F)
  Ek=XFromModList(ModList,outparam = "K.elas")[[ModNum]]
  Ep=XFromModList(ModList,outparam = "P.elas")[[ModNum]]
  Er=XFromModList(ModList,outparam = "R.elas")[[ModNum]]
  PlotKnorm(Ek,NORM_REC = F,colmax=max(max(Ek)),
            ImageTitle=paste0("K.elas: ",ModDF$SIG[ModNum],"\nL= ",
                              round(ModDF[,lambda_str][ModNum],3)))
  PlotKnorm(Ep,NORM_REC = F,colmax=max(max(Ep)),
            ImageTitle=paste0("P.elas: ",ModDF$SIG[ModNum],"\nL= ",
                              round(ModDF[,lambda_str][ModNum],3)))
  PlotKnorm(Er,NORM_REC = F,colmax=max(max(Er)),
            ImageTitle=paste0("R.elas: ",ModDF$SIG[ModNum],"\nL= ",
                              round(ModDF[,lambda_str][ModNum],3)))
  Elim=c(min(c(Ep,Er)),max(c(Ep,Er)))
  plot(colSums(Ep),colSums(Er),
       xlim=Elim,ylim=Elim,
       xlab="Elascity to Pi",ylab="Elasticity to Ri",
       main="Elasticity to P vs R")
  abline(0,1);
}


################################################################################
############# Loading in Data and defining mesh parameters #####################
################################################################################

# Dataframe named: ColonyLevel
# load to update ColonyLevel
load("2-VitalRateFunctions/Output/Colony_Data_edited.rdata") 
#View(ColonyLevel)


# Global mesh variables 

max.size <- 1.1*max(c(log10(ColonyLevel$StartingSize), 
                      log10(ColonyLevel$EndingSize)), na.rm = T)
# max ~ 3.75 m colony diameter

n = 50 #mesh size / number of cells in the discretized kernel
rec.size <- -0.1  
min.size <- -0.5  

# boundary points (the edges of the cells defining the kernel 
bin_size <- min.size + c(0:n) * (max.size - min.size)/n 

# mesh points (midpoints of cells)
y <- 0.5 * (bin_size[1:n]+bin_size[2:(n+1)])

I <- y >= rec.size
delta_size <- y[2] - y[1] #width of cells (h)


## Managing and Organizing Data ------------------------------------------------

Scenario= "Regional" #Change based on analysis: "Regional" "SIG"  

if(Scenario=="SIG"){
  load("2-VitalRateFunctions/Output/VR_Models__allmodfits_noNAs.rdata") #SIG VR functions used for IPM: ModelParams_SIG
}else if (Scenario=="Regional"){
  load("2-VitalRateFunctions/Output/VR_Models__allModFits_noNAs_regional.rdata") # regional VR functions used for IPM: ModelParams_Regional
}
#View(ModelParams_SIG)
View(ModelParams_Regional)


#FOR REGIONAL MODEL SKIP TO LINE 251

# # Check that we only want to keep Interval Years with more than N>MinTrans colonies
# # Find low representation intervals to drop

#minTrans <- 20
minTotTrans <- 42
# drop any extra rows with poor representation...
#ModelParams_SIG_Enough=subset(ModelParams_SIG, g.N>=minTrans&s.N>=minTrans)
ModelParams_SIG_Enough=subset(ModelParams_SIG, g.N+s.N >= minTotTrans)


################################################################################
###### Calculating Lambdas/SSDs and sensitivities ##############################
################################################################################


###### RUN SIG MODELS ##############################


BigMatParamList=list(n,delta_size,y,rec.size)
names(BigMatParamList)=c("n","ds","y","rec.size")

#run IPM (outputs lam,w,v,v.dot.w,sens,elas,K,Gk,Sk,Pk,K.elas,P.elas,R.elas,eK,eP,eR)


#calculate lambda values using MEAN site and sector rec to SIG models
MnSiteMods=RunIPMfromParamMatrix(ModPrm = ModelParams_SIG_Enough,RecValCol = "MN.RecVal_Site",IPMstructure = BigMatParamList)
ModelParams_SIG_Enough$lambda_SiteStock=unlist(XFromModList(MnSiteMods,outparam = "lam"))

MnSecMods=RunIPMfromParamMatrix(ModPrm = ModelParams_SIG_Enough,RecValCol = "MN.RecVal_Sec",IPMstructure = BigMatParamList)
ModelParams_SIG_Enough$lambda_SecStock=lamFromModList(MnSecMods)

LoSecMods=RunIPMfromParamMatrix(ModPrm = ModelParams_SIG_Enough,RecValCol = "LOW_CI95.RecVal_Sec",IPMstructure = BigMatParamList)
ModelParams_SIG_Enough$lambda_LowCI95_SecStock=lamFromModList(LoSecMods)

HiSecMods=RunIPMfromParamMatrix(ModPrm = ModelParams_SIG_Enough,RecValCol = "HIGH_CI95.RecVal_Sec",IPMstructure = BigMatParamList)
ModelParams_SIG_Enough$lambda_HighCI95_SecStock=lamFromModList(HiSecMods)
#combine site and sector rec calculated lambas into 1 dataframe
ModelParams_SIG_Enough$lambda_SiteSecMean=colMeans(rbind(ModelParams_SIG_Enough$lambda_SiteStock,
                                                         ModelParams_SIG_Enough$lambda_SecStock),na.rm=T)


#calculate lambda values using MEDIAN site and sector rec to SIG models
MnSecAllMods=RunIPMfromParamMatrix(ModPrm = ModelParams_SIG_Enough,RecValCol = "MD.RecVal_Sec_All",IPMstructure = BigMatParamList)
ModelParams_SIG_Enough$lambda_SecStockAll=lamFromModList(MnSecAllMods)

LoSecAllMods=RunIPMfromParamMatrix(ModPrm = ModelParams_SIG_Enough,RecValCol = "MD.CI95_LO.RecVal_Sec_All",IPMstructure = BigMatParamList)
ModelParams_SIG_Enough$lambda_LowCI95_SecStockAll=lamFromModList(LoSecAllMods)

HiSecAllMods=RunIPMfromParamMatrix(ModPrm = ModelParams_SIG_Enough,RecValCol = "MD.CI95_HI.RecVal_Sec_All",IPMstructure = BigMatParamList)
ModelParams_SIG_Enough$lambda_HighCI95_SecStockAll=lamFromModList(HiSecAllMods)

#calculate mid point date between sampling periods
ModelParams_SIG_Enough$MidPtDate=ModelParams_SIG_Enough$StartingDate+(ModelParams_SIG_Enough$EndingDate-ModelParams_SIG_Enough$StartingDate )/2
View(ModelParams_SIG_Enough)


#Save SIG lambda values to use in Temperature Stress Modeling
Date <- 01012023
save(ModelParams_SIG_Enough, file = sprintf("3-IPM_modeling/Output/all_SIGlambdas_%s.rdata",Date))
write.csv(ModelParams_SIG_Enough, "3-IPM_modeling/Output/SIGlambdas_allrectypes.csv", row.names = FALSE)

#load("./data/all_Lambdas_120621.rdata") #SIG lambda values: FINAL DATASET




###### RUN REGIONAL MODELS ##############################

BigMatParamList=list(n,delta_size,y,rec.size)
names(BigMatParamList)=c("n","ds","y","rec.size")

MnSiteMods=Region_RunIPMfromParamMatrix(ModPrm = ModelParams_Regional,RecValCol = "SiteStock_mean",IPMstructure = BigMatParamList)
ModelParams_Regional$lambda_SiteStock=unlist(XFromModList(MnSiteMods,outparam = "lam"))

LoSiteMods=Region_RunIPMfromParamMatrix(ModPrm = ModelParams_Regional,RecValCol = "SiteStock_q05",IPMstructure = BigMatParamList)
ModelParams_Regional$lambda_LowCI95_SiteStock=lamFromModList(LoSiteMods)

HiSiteMods=Region_RunIPMfromParamMatrix(ModPrm = ModelParams_Regional,RecValCol = "SiteStock_q95",IPMstructure = BigMatParamList)
ModelParams_Regional$lambda_HighCI95_SiteStock=lamFromModList(HiSiteMods)

MnSecMods=Region_RunIPMfromParamMatrix(ModPrm = ModelParams_Regional,RecValCol = "SectorStock_mean",IPMstructure = BigMatParamList)
ModelParams_Regional$lambda_SecStock=lamFromModList(MnSecMods)

LoSecMods=Region_RunIPMfromParamMatrix(ModPrm = ModelParams_Regional,RecValCol = "SectorStock_MN_CI95_LO",IPMstructure = BigMatParamList)
ModelParams_Regional$lambda_LowCI95_SecStock=lamFromModList(LoSecMods)

HiSecMods=Region_RunIPMfromParamMatrix(ModPrm = ModelParams_Regional,RecValCol = "SectorStock_MN_CI95_HI",IPMstructure = BigMatParamList)
ModelParams_Regional$lambda_HighCI95_SecStock=lamFromModList(HiSecMods)


#SecStockAll
MnSecAll_Mods=Region_RunIPMfromParamMatrix(ModPrm = ModelParams_Regional,RecValCol = "MN.RecVal_Sec_All",IPMstructure = BigMatParamList)
ModelParams_Regional$lambda_SecStockAll=lamFromModList(MnSecAll_Mods)



#subkernel Regional elas values
#MnSiteMods[[2]] #prints ALL of the IPM outputs
XFromModList(MnSecMods,outparam = "eK") #prints summed elas. should be 1
XFromModList(MnSecMods,outparam = "eP") #print the growth/surv elas
XFromModList(MnSecMods,outparam = "eR") #print the repro contribution (elas)


#Save Regional lambda values
Date <- 01012023
save(ModelParams_Regional, file = sprintf("3-IPM_modeling/Output/all_RegionalLambdas_%s.rdata",Date))

#load("3-IPM_modeling/Output/all_RegionalLambdas_120321.rdata") #Regional lambda values: FINAL DATASET


################################################################################
###### Plotting kernels, elasticities and sensitivities ##############################
################################################################################

#Kernel plots for 3 high and 3 low lambda values
par(mfrow=c(2,3))
PlotN=2
PlotKnorm(XFromModList(MnSecAllMods,outparam = "K")[[PlotN]],
          ImageTitle=paste0("K: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SecStockAll[PlotN],3)),NORM_REC = T)
PlotN=29
PlotKnorm(XFromModList(MnSecAllMods,outparam = "K")[[PlotN]],
          ImageTitle=paste0("K: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SecStockAll[PlotN],3)),NORM_REC = T)
PlotN=31
PlotKnorm(XFromModList(MnSecAllMods,outparam = "K")[[PlotN]],
          ImageTitle=paste0("K: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SecStockAll[PlotN],3)),NORM_REC = T)
PlotN=10
PlotKnorm(XFromModList(MnSecAllMods,outparam = "K")[[PlotN]],
          ImageTitle=paste0("K: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SecStockAll[PlotN],3)),NORM_REC = T)
PlotN=15
PlotKnorm(XFromModList(MnSecAllMods,outparam = "K")[[PlotN]],
          ImageTitle=paste0("K: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SecStockAll[PlotN],3)),NORM_REC = T)
PlotN=8
PlotKnorm(XFromModList(MnSecAllMods,outparam = "K")[[PlotN]],
          ImageTitle=paste0("K: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SecStockAll[PlotN],3)),NORM_REC = T)



#REGIONAL KERNELS
par(mfrow=c(2,3))
PlotN=1
PlotKnorm(XFromModList(MnSecMods,outparam = "K")[[PlotN]],
          ImageTitle=paste0(ModelParams_Regional$Genus_Code[PlotN],", ",ModelParams_Regional$REGION[PlotN],
                            "\nLambda= ",round(ModelParams_Regional$lambda_SecStock[PlotN],3)))
PlotN=2
PlotKnorm(XFromModList(MnSecMods,outparam = "K")[[PlotN]],
          ImageTitle=paste0(ModelParams_Regional$Genus_Code[PlotN],", ",ModelParams_Regional$REGION[PlotN],
                            "\nLambda= ",round(ModelParams_Regional$lambda_SecStock[PlotN],3)))
PlotN=3
PlotKnorm(XFromModList(MnSecMods,outparam = "K")[[PlotN]],
          ImageTitle=paste0(ModelParams_Regional$Genus_Code[PlotN],", ",ModelParams_Regional$REGION[PlotN],
                            "\nLambda= ",round(ModelParams_Regional$lambda_SecStock[PlotN],3)))
PlotN=5
PlotKnorm(XFromModList(MnSecMods,outparam = "K")[[PlotN]],
          ImageTitle=paste0(ModelParams_Regional$Genus_Code[PlotN],", ",ModelParams_Regional$REGION[PlotN],
                            "\nLambda= ",round(ModelParams_Regional$lambda_SecStock[PlotN],3)))
PlotN=4
PlotKnorm(XFromModList(MnSecMods,outparam = "K")[[PlotN]],
          ImageTitle=paste0(ModelParams_Regional$Genus_Code[PlotN],", ",ModelParams_Regional$REGION[PlotN],
                            "\nLambda= ",round(ModelParams_Regional$lambda_SecStock[PlotN],3)))
PlotN=6
PlotKnorm(XFromModList(MnSecMods,outparam = "K")[[PlotN]],
          ImageTitle=paste0(ModelParams_Regional$Genus_Code[PlotN],", ",ModelParams_Regional$REGION[PlotN],
                            "\nLambda= ",round(ModelParams_Regional$lambda_SecStock[PlotN],3)))


#elasticity plots
par(mfrow=c(2,3))
PlotN=29
PlotKnorm(XFromModList(MnSecAllMods,outparam = "K.elas")[[PlotN]],
          ImageTitle=paste0("K.elas: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SecStock[PlotN],3)),)
PlotN=31
PlotKnorm(XFromModList(MnSecAllMods,outparam = "K.elas")[[PlotN]],
          ImageTitle=paste0("K.elas: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SecStock[PlotN],3)))
PlotN=10
PlotKnorm(XFromModList(MnSecAllMods,outparam = "K.elas")[[PlotN]],
          ImageTitle=paste0("K.elas: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SecStockAll[PlotN],3)))
PlotN=15
PlotKnorm(XFromModList(MnSecAllMods,outparam = "K.elas")[[PlotN]],
          ImageTitle=paste0("K.elas: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SecStockAll[PlotN],3)))
PlotN=8
PlotKnorm(XFromModList(MnSecAllMods,outparam = "K.elas")[[PlotN]],
          ImageTitle=paste0("K.elas: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SecStockAll[PlotN],3)))

#REGIONAL elasticities
#Change
#P.elas   #R.elas   #K.elas

#png(filename = "Output/RegionalElasticities_SectorStock.png")
par(mfrow=c(2,3))
PlotN=1
PlotKnorm(XFromModList(MnSecMods,outparam = "K.elas")[[PlotN]],
          ImageTitle=paste0(ModelParams_Regional$Genus_Code[PlotN],", ",ModelParams_Regional$REGION[PlotN],
                            "\nLambda= ",round(ModelParams_Regional$lambda_SecStock[PlotN],3)))

PlotN=2
PlotKnorm(XFromModList(MnSecMods,outparam = "K.elas")[[PlotN]],
          ImageTitle=paste0(ModelParams_Regional$Genus_Code[PlotN],", ",ModelParams_Regional$REGION[PlotN],
                            "\nLambda= ",round(ModelParams_Regional$lambda_SecStock[PlotN],3)))

PlotN=3
PlotKnorm(XFromModList(MnSecMods,outparam = "K.elas")[[PlotN]],
          ImageTitle=paste0(ModelParams_Regional$Genus_Code[PlotN],", ",ModelParams_Regional$REGION[PlotN],
                            "\nLambda= ",round(ModelParams_Regional$lambda_SecStock[PlotN],3)))
PlotN=5
PlotKnorm(XFromModList(MnSecMods,outparam = "R.elas")[[PlotN]],
          ImageTitle=paste0(ModelParams_Regional$Genus_Code[PlotN],", ",ModelParams_Regional$REGION[PlotN],
                            "\nLambda= ",round(ModelParams_Regional$lambda_SecStock[PlotN],3)))
PlotN=4
PlotKnorm(XFromModList(MnSecMods,outparam = "R.elas")[[PlotN]],
          ImageTitle=paste0(ModelParams_Regional$Genus_Code[PlotN],", ",ModelParams_Regional$REGION[PlotN],
                            "\nLambda= ",round(ModelParams_Regional$lambda_SecStock[PlotN],3)))
PlotN=6
PlotKnorm(XFromModList(MnSecMods,outparam = "R.elas")[[PlotN]],
          ImageTitle=paste0(ModelParams_Regional$Genus_Code[PlotN],", ",ModelParams_Regional$REGION[PlotN],
                            "\nLambda= ",round(ModelParams_Regional$lambda_SecStock[PlotN],3)))
#dev.off()




#PLOT sensitivities for 3 high and 3 low lambda values
par(mfrow=c(2,3))
PlotN=2
PlotKnorm(XFromModList(MnSiteMods,outparam = "sens")[[PlotN]],
          ImageTitle=paste0("sens: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SiteStock[PlotN],3)))
PlotN=29
PlotKnorm(XFromModList(MnSiteMods,outparam = "sens")[[PlotN]],
          ImageTitle=paste0("sens: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SiteStock[PlotN],3)))
PlotN=31
PlotKnorm(XFromModList(MnSiteMods,outparam = "sens")[[PlotN]],
          ImageTitle=paste0("sens: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SiteStock[PlotN],3)))
PlotN=10
PlotKnorm(XFromModList(MnSiteMods,outparam = "sens")[[PlotN]],
          ImageTitle=paste0("sens: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SiteStock[PlotN],3)))
PlotN=15
PlotKnorm(XFromModList(MnSiteMods,outparam = "sens")[[PlotN]],
          ImageTitle=paste0("sens: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SiteStock[PlotN],3)))
PlotN=8
PlotKnorm(XFromModList(MnSiteMods,outparam = "sens")[[PlotN]],
          ImageTitle=paste0("sens: ",ModelParams_SIG_Enough$SIG[PlotN],"\nL= ",
                            round(ModelParams_SIG_Enough$lambda_SiteStock[PlotN],3)))



################################################################################
###### Plotting full model reports ##############################
################################################################################


outpath= "3-IPM_modeling/Output/"
imgsc=1.25
dev.off()

#print full model reports for all SIG models
#Using Site sector rec
ModDF=ModelParams_SIG_Enough
ModList=MnSiteMods
lam="lambda_SiteStock"
i_ord=order(ModDF[,lam])
i_ord=i_ord[which(!is.na(ModDF[,lam][i_ord]))]
for(ModNum_i in i_ord){
  jpeg(filename = paste0(outpath,"L",sprintf("%1.3f",ModDF[,lam][ModNum_i]),
                         ModDF$SIG[ModNum_i],"_RecSite.jpg"),
       width=imgsc*800,height=imgsc*800)
  FullModelReport(ModList = ModList,ModDF = ModDF,ModNum=ModNum_i,lambda_str = lam)
  dev.off()
  print(ModNum_i)
}


ModList = MnSiteMods

names(ModList[[15]])
Ek=XFromModList(ModList,outparam = "K.elas")[[15]]
Ep=XFromModList(ModList,outparam = "P.elas")[[15]]
Er=XFromModList(ModList,outparam = "R.elas")[[15]]
#Eg=XFromModList(ModList,outparam = "G.elas")[[15]]
sum(Ek)
sum(Ep)
sum(Er)
#sum(Eg)
sum(Ep)+sum(Er)


FullModelReport(ModList = MnSiteMods,ModDF = ModelParams_SIG_Enough,ModNum=28)

w=XFromModList(MnSiteMods,outparam = "w")[[28]]
v=XFromModList(MnSiteMods,outparam = "v")[[28]]
v.dot.w=XFromModList(MnSiteMods,outparam = "v.dot.w")[[28]]
sum(w*v)
v.dot.w



#Regional model full model reports
#outpath="C:/Users/Thomas.Oliver/WORK/Projects/FY22/Vital Rates - SFM/SIG_FullModelReports/"
imgsc=1.25
dev.off()

RegMods=Region_RunIPMfromParamMatrix(ModPrm = ModelParams_Regional,RecValCol = "SiteStock_mean",IPMstructure = BigMatParamList)

ModDF=ModelParams_Regional
ModList=RegMods 
lam="lambda_SiteStock"
i_ord=order(ModDF[,lam])
i_ord=i_ord[which(!is.na(ModDF[,lam][i_ord]))]
for(ModNum_i in i_ord){
  jpeg(filename = paste0(outpath,"L",sprintf("%1.3f",ModDF[,lam][ModNum_i]),
                         ModDF$Genus_Code[ModNum_i],"_Regional.jpg"),
       width=imgsc*800,height=imgsc*800)
  FullModelReport(ModList = ModList,ModDF = ModDF,ModNum=ModNum_i,lambda_str = lam)
  dev.off()
  print(ModNum_i)
}


######################################################################################




##########################
#PLOTTING Rec values: Site-stock rec, sector-stock rec, and modeled sector-stock rec values

SecSite_Rec=ggplot(ModelParams_SIG_Enough)+
  geom_jitter(aes(x=0.001,y=MN.RecVal_Site),height=0,width=0.015,color="gray",shape=15)+
  geom_jitter(aes(x=MN.RecVal_Sec,y=0.001),height=0.015,width=0,color="gray",shape=15)+
  geom_point(aes(x=MN.RecVal_Sec,y=MN.RecVal_Site),color="blue")+
  #geom_errorbarh(aes(y=MN.RecVal_Site,xmin=LOW_CI95.RecVal_Sec,xmax=HIGH_CI95.RecVal_Sec),color="blue")+
  geom_abline()+ theme_bw()+
  xlab("Sector-Level")+ylab("Site-Level")+
  ggtitle("Site vs Sector Recruitment")+ coord_equal()+
  xlim(c(-.01,.08))+ylim(c(-.01,.02))+
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        strip.text.y= element_text(size = 12),
        plot.title = element_text(size = 15) )
  #scale_x_continuous(trans = 'log10')+scale_y_continuous(trans = 'log10')


SecAll_Site_Rec=ggplot(ModelParams_SIG_Enough)+
  geom_jitter(aes(x=0.001,y=MN.RecVal_Site),height=0,width=0.015,color="gray",shape=15)+
  geom_jitter(aes(x=MD.RecVal_Sec_All,y=0.001),height=0.015,width=0,color="gray",shape=15)+
  geom_point(aes(x=MD.RecVal_Sec_All,y=MN.RecVal_Site),color="blue")+
  #geom_errorbarh(aes(y=MN.RecVal_Site,xmin=MD.CI95_LO.RecVal_Sec_All,xmax=MD.CI95_HI.RecVal_Sec_All),color="blue",width=0.1)+ #xvalue
  geom_abline()+ theme_bw()+
  xlab("All-Sectors")+ylab("Site-Level")+
  ggtitle("Site vs All-Sectors Recruitment")+ coord_equal()+
  xlim(c(-.01,.08))+ylim(c(-.01,.02))+
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        strip.text.y= element_text(size = 12),
        plot.title = element_text(size = 15) )
  #scale_x_continuous(trans = 'log10')+scale_y_continuous(trans = 'log10')+
  

Sec_SecAll_Rec=ggplot(ModelParams_SIG_Enough)+
  geom_jitter(aes(x=0.001,y=MD.RecVal_Sec_All),height=0,width=0.015,color="gray",shape=15)+
  geom_jitter(aes(x=MN.RecVal_Sec,y=0.001),height=0.015,width=0,color="gray",shape=15)+
  geom_point(aes(x=MN.RecVal_Sec,y=MD.RecVal_Sec_All),color="blue")+
  #geom_errorbarh(aes(y=MD.RecVal_Sec_All,xmin=LOW_CI95.RecVal_Sec,xmax=HIGH_CI95.RecVal_Sec),color="blue")+ #xvalue error
  #geom_errorbar(aes(x=MN.RecVal_Sec,ymin=MD.CI95_LO.RecVal_Sec_All,ymax=MD.CI95_HI.RecVal_Sec_All),color="blue")+
  geom_abline()+ theme_bw()+
  xlab("Sector-Level")+ylab("All-Sectors")+
  ggtitle("Sector vs All-Sectors Recruitment")+ coord_equal()+
  xlim(c(-.01,.08))+ylim(c(-.01,.02))+
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        strip.text.y= element_text(size = 12),
        plot.title = element_text(size = 15) )
  #scale_x_continuous(trans = 'log10')+scale_y_continuous(trans = 'log10')

SecSite_Rec+SecAll_Site_Rec+Sec_SecAll_Rec

hist(ModelParams_SIG_Enough$MN.RecVal_Site)
hist(ModelParams_SIG_Enough$MN.RecVal_Sec)
hist(ModelParams_SIG_Enough$MD.RecVal_Sec_All)


###########################


#Comparing lambdas calculated using Site-stock rec, sector-stock rec, and modeled sector-stock rec values

Sec_Site=ggplot(ModelParams_SIG_Enough)+
  geom_jitter(aes(x=0.6,y=lambda_SiteStock),height=0,width=0.025,color="gray")+
  geom_jitter(aes(x=lambda_SecStock,y=0.6),height=0.025,width=0,color="gray")+
  geom_point(aes(x=lambda_SecStock,y=lambda_SiteStock),color="blue")+
  geom_errorbarh(aes(y=lambda_SiteStock,xmin=lambda_LowCI95_SecStock,xmax=lambda_HighCI95_SecStock),color="blue")+
  geom_vline(xintercept = 1)+
  geom_hline(yintercept = 1)+
  geom_abline()+
  xlab("Lambda Using Sector-Level Observed Recruitment")+
  ylab("Lambda Using Site-Level Recruitment")+
  ggtitle("Population Growth Rate: Site vs Sector (obs.) Recruitment")+xlim(c(.5,1.5))+ylim(c(.5,1.5))+coord_equal()

SecAll_Site=ggplot(ModelParams_SIG_Enough)+
  geom_jitter(aes(x=0.6,y=lambda_SiteStock),height=0,width=0.025,color="gray")+
  geom_jitter(aes(x=lambda_SecStockAll,y=0.6),height=0.025,width=0,color="gray")+
  geom_point(aes(x=lambda_SecStockAll,y=lambda_SiteStock),color="blue")+
  geom_errorbarh(aes(y=lambda_SiteStock,xmin=lambda_LowCI95_SecStockAll,xmax=lambda_HighCI95_SecStockAll),color="blue")+
  geom_vline(xintercept = 1)+
  geom_hline(yintercept = 1)+
  geom_abline()+
  xlab("Lambda Using Sector-Level Modeled Recruitment")+
  ylab("Lambda Using Site-Level Recruitment")+
  ggtitle("Population Growth Rate: Site vs Sector (mod.) Recruitment")+xlim(c(.5,1.5))+ylim(c(.5,1.5))+coord_equal()

Sec_SecAll=ggplot(ModelParams_SIG_Enough)+
  geom_jitter(aes(x=0.6,y=lambda_SecStockAll),height=0,width=0.025,color="gray")+
  geom_jitter(aes(x=lambda_SecStock,y=0.6),height=0.025,width=0,color="gray")+
  geom_point(aes(x=lambda_SecStock,y=lambda_SecStockAll),color="blue")+
  geom_errorbarh(aes(y=lambda_SecStockAll,xmin=lambda_LowCI95_SecStock,xmax=lambda_HighCI95_SecStock),color="blue")+
  geom_errorbar(aes(x=lambda_SecStock,ymin=lambda_LowCI95_SecStockAll,ymax=lambda_HighCI95_SecStockAll),color="blue")+
  geom_vline(xintercept = 1)+
  geom_hline(yintercept = 1)+
  geom_abline()+
  xlab("Lambda Using Sector-Level Observed Recruitment")+
  ylab("Lambda Using Sector-Level Modeled Recruitment")+
  ggtitle("Population Growth Rate: Sector (obs.) vs Sector (mod.) Recruitment")+xlim(c(.5,1.5))+ylim(c(.5,1.5))

Sec_Site + SecAll_Site + Sec_SecAll



setcols=c(brewer.pal(n = 3,name = "Dark2")[1:2],"gray50")

LamPlot=ggplot(ModelParams_SIG_Enough)+
  geom_rect(aes(xmin=ymd("2014-09-15"),xmax=ymd("2014-11-15"),ymin=-Inf,ymax=Inf),fill="gray85")+
  geom_rect(aes(xmin=ymd("2015-09-15"),xmax=ymd("2015-11-15"),ymin=-Inf,ymax=Inf),fill="gray85")+
  geom_rect(aes(xmin=ymd("2019-09-15"),xmax=ymd("2019-11-15"),ymin=-Inf,ymax=Inf),fill="gray85")+
  geom_errorbar(aes(x=MidPtDate,ymin=lambda_LowCI95_SecStockAll,ymax=lambda_HighCI95_SecStockAll,color="gray75"))+
  geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SecStockAll,yend=lambda_SecStockAll,color="gray75"),alpha=.25)+
  geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SiteStock,yend=lambda_SiteStock,color="orange"),alpha=.25)+
  geom_errorbar(aes(x=MidPtDate,ymin=lambda_LowCI95_SecStock,ymax=lambda_HighCI95_SecStock,color="darkgreen"))+
  geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SecStock,yend=lambda_SecStock,color="darkgreen"),alpha=.25)+
  geom_point(aes(x=MidPtDate,y=lambda_SecStockAll,color="gray75"))+ #,shape=Genus_Code
  geom_point(aes(x=MidPtDate,y=lambda_SiteStock,color="orange"))+
  geom_point(aes(x=MidPtDate,y=lambda_SecStock,color="darkgreen"))+
  geom_hline(yintercept = 1)+
  scale_shape_discrete(name="Genus",labels=c("Montipora","Pocillopora","Porites"))+
  scale_color_identity(guide="legend",name="Recruitment Parameterization:",
                       breaks=c("orange","darkgreen","gray75"),
                       labels=c("Site-Scale (observed)","Sector-Scale (observed)","Sector-Scale (regional genus median)"))+
  xlab("Date")+ ylab("Population Growth Rate (Lambda)")+
  scale_y_continuous(breaks=c(0.25,0.50,0.75,1.0,1.25,1.5))+
  facet_grid(Genus_Code~.)+
  #facet_grid(REGION~Genus_Code)+
  ggtitle("Reef Coral Population Growth Rates: Site, Interval, Genus")+
  theme_bw()+
  stat_smooth(aes(x=MidPtDate,y=lambda_SecStockAll),alpha=.25,method = "loess",span=.7)+
  theme(legend.position = "bottom",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"))


#SIG lambdas temporally: lambdas from 3 rec values by region
mhi= ModelParams_SIG_Enough %>%
  filter(REGION=="MHI") %>%
  ggplot()+
  geom_rect(aes(xmin=ymd("2014-09-15"),xmax=ymd("2014-11-15"),ymin=-Inf,ymax=Inf),fill="gray85")+
  geom_rect(aes(xmin=ymd("2015-09-15"),xmax=ymd("2015-11-15"),ymin=-Inf,ymax=Inf),fill="gray85")+
  geom_rect(aes(xmin=ymd("2019-09-15"),xmax=ymd("2019-11-15"),ymin=-Inf,ymax=Inf),fill="gray85")+
  geom_errorbar(aes(x=MidPtDate,ymin=lambda_LowCI95_SecStockAll,ymax=lambda_HighCI95_SecStockAll,color="blue"))+
  geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SecStockAll,yend=lambda_SecStockAll,color="blue"),alpha=.25)+
  geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SiteStock,yend=lambda_SiteStock,color="orange"),alpha=.25)+
  geom_errorbar(aes(x=MidPtDate,ymin=lambda_LowCI95_SecStock,ymax=lambda_HighCI95_SecStock,color="darkgreen"))+
  geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SecStock,yend=lambda_SecStock,color="darkgreen"),alpha=.25)+
  geom_point(aes(x=MidPtDate,y=lambda_SecStockAll,color="blue"))+ #,shape=Genus_Code
  geom_point(aes(x=MidPtDate,y=lambda_SiteStock,color="orange"))+
  geom_point(aes(x=MidPtDate,y=lambda_SecStock,color="darkgreen"))+
  geom_hline(yintercept = 1)+
  coord_cartesian(ylim = c(0.25,1.5))+
  scale_y_continuous(breaks=c(0.25,0.50,0.75,1.0,1.25,1.5))+
  scale_color_identity(guide="legend",name="Recruitment Parameterization:",
                       breaks=c("orange","darkgreen","blue"),
                       labels=c("Site-Level","Sector-Level","All-Sectors"))+
  xlab("Date")+ ylab("Population Growth Rate (Lambda)")+
  facet_grid(Genus_Code~.)+
  ggtitle("MHI Site-Interval-Genus Population Growth Rates")+
  theme_bw()+
  stat_smooth(aes(x=MidPtDate,y=lambda_SecStockAll),alpha=.25,method = "loess",span=.7,color="blue", se=FALSE)+
  stat_smooth(aes(x=MidPtDate,y=lambda_SiteStock),alpha=.25,method = "loess",span=.7,color="orange",se=FALSE)+
  #stat_smooth(aes(x=MidPtDate,y=lambda_SecStock),alpha=.25,method = "loess",span=.7,color="darkgreen", se=FALSE)+
  theme(legend.position = "bottom",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text= element_text(size = 10))
nwhi= ModelParams_SIG_Enough %>%
  filter(REGION=="NWHI") %>%
  ggplot()+
  geom_rect(aes(xmin=ymd("2014-09-15"),xmax=ymd("2014-11-15"),ymin=-Inf,ymax=Inf),fill="gray85")+
  geom_rect(aes(xmin=ymd("2015-09-15"),xmax=ymd("2015-11-15"),ymin=-Inf,ymax=Inf),fill="gray85")+
  geom_rect(aes(xmin=ymd("2019-09-15"),xmax=ymd("2019-11-15"),ymin=-Inf,ymax=Inf),fill="gray85")+
  geom_errorbar(aes(x=MidPtDate,ymin=lambda_LowCI95_SecStockAll,ymax=lambda_HighCI95_SecStockAll,color="blue"))+ #vertical error bar
  geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SecStockAll,yend=lambda_SecStockAll,color="blue"),alpha=.25)+ #horiz time interval
  geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SiteStock,yend=lambda_SiteStock,color="orange"),alpha=.25)+
  geom_errorbar(aes(x=MidPtDate,ymin=lambda_LowCI95_SecStock,ymax=lambda_HighCI95_SecStock,color="darkgreen"))+
  geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SecStock,yend=lambda_SecStock,color="darkgreen"),alpha=.25)+
  geom_point(aes(x=MidPtDate,y=lambda_SecStockAll,color="blue"))+ #,shape=Genus_Code
  geom_point(aes(x=MidPtDate,y=lambda_SiteStock,color="orange"))+
  geom_point(aes(x=MidPtDate,y=lambda_SecStock,color="darkgreen"))+
  geom_hline(yintercept = 1)+
  coord_cartesian(ylim = c(0.25,1.5))+
  scale_y_continuous(breaks=c(0.25,0.50,0.75,1.0,1.25,1.5))+
  #ylim(0.25,1.5)+
  scale_color_identity(guide="legend",name="Recruitment Parameterization:",
                       breaks=c("orange","darkgreen","blue"),
                       labels=c("Site-Level","Sector-Level","All-Sectors"))+
  xlab("Date")+ ylab("Population Growth Rate (Lambda)")+
  facet_grid(Genus_Code~.)+
  ggtitle("NWHI Site-Interval-Genus Population Growth Rates")+
  #stat_smooth(aes(x=MidPtDate,y=lambda_SecStockAll),alpha=.25,method = "loess",span=.7,color="blue")+
  stat_smooth(aes(x=MidPtDate,y=lambda_SiteStock),alpha=.25,method = "loess",span=.7,color="orange")+
  #stat_smooth(aes(x=MidPtDate,y=lambda_SecStock),alpha=.25,method = "loess",span=.7,color="darkgreen")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.text= element_text(size = 10))

mhi+nwhi
#ggsave(plot = TESTING, filename = "Figures/SIG_lambas_SectorStockAll.png",width = 6, height = 4)


ggplot(ModelParams_SIG_Enough)+
  # geom_rect(aes(xmin=StartingDate,xmax=EndingDate,ymin=lambda_LowCI95_SecStock,ymax=lambda_HighCI95_SecStock),
  #           fill=setcols[2],alpha=0.1)+
  geom_rect(aes(xmin=StartingDate,xmax=EndingDate,ymin=lambda_LowCI95_SecStockAll,ymax=lambda_HighCI95_SecStockAll),
            fill=setcols[3],alpha=0.1)+
  geom_rect(aes(xmin=ymd("2014-09-15"),xmax=ymd("2014-11-15"),ymin=-Inf,ymax=Inf),fill="pink")+
  geom_rect(aes(xmin=ymd("2015-09-15"),xmax=ymd("2015-11-15"),ymin=-Inf,ymax=Inf),
            fill="pink")+
  geom_rect(aes(xmin=ymd("2019-09-15"),xmax=ymd("2019-11-15"),ymin=-Inf,ymax=Inf),
            fill="pink")+
   geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SiteStock,yend=lambda_SiteStock),color=setcols[1])+
  # geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SecStock,yend=lambda_SecStock),color=setcols[2])+
  geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SecStockAll,yend=lambda_SecStockAll),color=setcols[3])+
  geom_hline(yintercept = 1)+
  scale_color_manual(breaks=c("Site","Sec","SecAll"),
                     values = c("Site"=setcols[1],"Sec"=setcols[2],"SecAll"=setcols[3]),
                     aesthetics = c('color','fill'))+
  xlab("Date")+
  ylab("Lambda")+
  facet_grid(REGION~Genus_Code)+
  ggtitle("Lambda")+theme_bw()#+ theme(legend.position = c(0.8, 0.2))

#site rec lambdas for MOSP and POSP
ModelParams_SIG_Enough %>%
  filter(Genus_Code%in%c("MOSP","POSP")&REGION=="MHI") %>%
  ggplot()+
  # geom_rect(aes(xmin=StartingDate,xmax=EndingDate,ymin=lambda_LowCI95_SecStock,ymax=lambda_HighCI95_SecStock),
  #           fill=setcols[2],alpha=0.1)+
  # geom_rect(aes(xmin=StartingDate,xmax=EndingDate,ymin=lambda_LowCI95_SecStockAll,ymax=lambda_HighCI95_SecStockAll),
  #           fill=setcols[3],alpha=0.1)+
  geom_rect(aes(xmin=ymd("2014-09-15"),xmax=ymd("2014-11-15"),ymin=-Inf,ymax=Inf),fill="pink")+
  geom_rect(aes(xmin=ymd("2015-09-15"),xmax=ymd("2015-11-15"),ymin=-Inf,ymax=Inf),fill="pink")+
  geom_rect(aes(xmin=ymd("2019-09-15"),xmax=ymd("2019-11-15"),ymin=-Inf,ymax=Inf),fill="pink")+
  # geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SiteStock,yend=lambda_SiteStock),color=setcols[1])+
  # geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SecStock,yend=lambda_SecStock),color=setcols[2])+
  # geom_segment(aes(x=StartingDate,xend=EndingDate,y=lambda_SecStockAll,yend=lambda_SecStockAll),color=setcols[3])+
  # geom_hline(yintercept = 1)+
  geom_point(aes(x=StartingDate,y=lambda_SiteStock))+
  scale_color_manual(guide="legend",breaks=c("Site","Sec","SecAll"),
                     values = c("Site"=setcols[1],"Sec"=setcols[2],"SecAll"=setcols[3]),
                     aesthetics = c('color','fill'))+
  xlab("Date")+
  ylab("Lambda")+
  facet_grid(Genus_Code~.)+
  ggtitle("Lambda")+theme_bw()+ theme(legend.position = c(0.5, 0.5))


#site vs sec lambda histogram by genus
siteH=ggplot(ModelParams_SIG_Enough)+geom_histogram(aes(x=lambda_SiteStock,fill=Genus_Code))+geom_vline(xintercept = 1)+xlim(c(.35,1.6))
secH=ggplot(ModelParams_SIG_Enough)+geom_histogram(aes(x=lambda_SecStock,fill=Genus_Code))+geom_vline(xintercept = 1)+xlim(c(.35,1.6))
siteH/secH













################################################################################
###### Tom's stuff: restoration # nubbins needed ##############################
################################################################################


#SIG models
A=10^y
D=sqrt(A/pi)
Taxon_i=1
SIG_Elas_DF=NULL
ModDF=ModelParams_SIG_Enough
ModList=MnSecAllMods

#restoration # nubbins stuff for SIG models
ModDF$StartingDate=ymd(ModDF$StartingDate)
ModDF$EndingDate=ymd(ModDF$EndingDate)
ModDF$MidPtYear=year(ModDF$StartingDate+difftime(ModDF$EndingDate,ModDF$StartingDate)/2)
for (SIG_i in 1:length(MnSecAllMods)){
  Ek_i=colSums(MnSecAllMods[[SIG_i]]$K.elas)
  Ekn_i=Ek_i/max(Ek_i)
  SIG_Elas_DF=rbind(SIG_Elas_DF,
                data.frame(Island=ModDF$ISLAND[SIG_i],
                           Site=ModDF$Site[SIG_i],
                           Interval=ModDF$Interval[SIG_i],
                           IntMid=ModDF$MidPtYear[SIG_i],
                           Taxon=ModDF$Genus_Code[SIG_i],Sz_i=y,A_i=A,D_i=D,Ek_i,Ekn_i,N_Nubs=1/Ek_i,N_NubsN=1/Ekn_i))
}

#restoration # nubbins stuff for regional models
Reg_Elas_DF=NULL
ModDF=data_regional
for (Taxon_i in 1:3){
  Ek_i=colSums(RegMods[[Taxon_i]]$K.elas)
  Ekn_i=Ek_i/max(Ek_i)
  Reg_Elas_DF=rbind(Reg_Elas_DF,
                data.frame(Taxon=ModDF$Genus_Code[Taxon_i],
                           Sz_i=y,A_i=A,D_i=D,Ek_i,Ekn_i,N_Nubs=1/Ek_i,N_NubsN=1/Ekn_i))
}

AA=SIG_Elas_DF %>% group_by(Site,Interval,Taxon) %>% summarize(Sz.MaxE=D_i[which.max(Ek_i)])
ggplot(AA,aes(x=Sz.MaxE,fill=Taxon))+geom_histogram()+facet_grid("Taxon")+scale_x_log10()

# 
# ggplot(Elas_DF,aes(x=D,y=Ekn_i))+
#   geom_point()+
#   geom_path()+
#   scale_x_log10()+
#   xlab("Coral Max Diameter (cm)")+
#   ylab("Normalized Elasticity of Lambda")+
#   ggtitle(paste0("Elasticity of Lambda By Size Class for ",ModDF$Genus_Code[Taxon_i]))+
#   theme_bw()
#   

Elas

#restoration stuff: elasticities for regional models
ggplot(Reg_Elas_DF,aes(x=D_i,y=Ekn_i,color=Taxon))+
  geom_point()+
  geom_path()+
  scale_x_log10(breaks=c(.5,2,5,20,50,200))+
#  scale_y_log10(breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000))+
  #geom_vline(xintercept = D[which.max(Ek_i)])+
  #facet_grid(Island~Interval)+
  xlab("Coral Max Diameter (cm)")+
  ylab("Ek_i")+
  ggtitle(paste0("Regional Model Elasticity of Lambda to Colony Size (i.e. Column Sum from Ek, Normalized)"))+
  theme_bw()

#restoration stuff: # nubbins required for effective restoration
ggplot(Reg_Elas_DF,aes(x=D_i,y=N_Nubs,color=Taxon))+
  geom_point(alpha=1)+
  geom_path()+
  #stat_smooth()+
  #scale_x_log10(breaks=c(.5,1,2,5,10,20,50,100,200),limits=c(.3,200))+  
  scale_x_continuous(limits=c(.5,20))+
  scale_y_log10(breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000),limits=c(1,750))+
  #geom_vline(xintercept = D[which.max(Ek_i)])+
  #facet_wrap("IntMid")+
  xlab("Coral Max Diameter (cm)")+
  ylab("N. Nubbins")+
  ggtitle(paste0("Number of Nubbins Required to Match Optimal Size for Increasing Population Growth",
                 "\n'N outplants of size i would be as effective as 1 nubbin of max elasticity'"))+
  theme_bw()

plot(D,Ek_i/max(Ek_i),type="l")

