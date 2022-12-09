library(plyr)
library(lubridate)
library(ggplot2)

load(file = "1-Annotations_ColonyTransitions/Output/PatchTransitions.rdata")
VitalRate_Growth$TransitionType=factor(VitalRate_Growth$TransitionType,levels=c("RECR","MORT","GROWTH","SHRINK","FISSION","FUSION","FUSION_FISSION"))
#Aggregate Colonies
CID=substr(VitalRate_Growth$ColonyID,1,7)
CID2=paste0(substr(VitalRate_Growth$T0_PatchName,1,17),CID,substr(VitalRate_Growth$T0_PatchName,25,26))
RM_i=which(substr(CID2,1,4)%in%c("MORT","RECR"))
CID2[RM_i]=paste0(substr(VitalRate_Growth$T0_PatchName[RM_i],6,22),CID[RM_i],substr(VitalRate_Growth$T0_PatchName[RM_i],30,31))
VitalRate_Growth$ColonyID=CID2
VitalRate_Patch=subset(VitalRate_Growth,DataOrError=="DATA")  

head(VitalRate_Patch,5)

PatchSize=data.frame(Name=c(as.vector(VitalRate_Patch$T0_PatchName),as.vector(VitalRate_Patch$T1_PatchName)),
                     Size=c(VitalRate_Patch$StartingSize,VitalRate_Patch$EndingSize),
                     Perim=c(VitalRate_Patch$StartingPerimeter,VitalRate_Patch$EndingPerimeter),
                     Diam=c(VitalRate_Patch$StartingMaxDiam,VitalRate_Patch$EndingMaxDiam))

dim(PatchSize)
PatchSize=unique(PatchSize)
dim(PatchSize)
PatchSizeLU=PatchSize$Size
names(PatchSizeLU)=PatchSize$Name
PatchPerimLU=PatchSize$Perim
names(PatchPerimLU)=PatchSize$Name
PatchDiamLU=PatchSize$Diam
names(PatchDiamLU)=PatchSize$Name
VitalRate_Colony=ddply(VitalRate_Patch,.(Site,DataOrError,ColonyID,Spec_Code,Genus_Code,StartingDate,EndingDate,Interval_Years),summarize,
                       N_t0=length(unique(T0_PatchName)),
                       N_t1=length(unique(T1_PatchName)),
                       StartingSize=sum(PatchSizeLU[unique(as.vector(T0_PatchName))],na.rm=T),
                       EndingSize=sum(PatchSizeLU[unique(as.vector(T1_PatchName))],na.rm=T),
                       StartingPerim=sum(PatchPerimLU[unique(as.vector(T0_PatchName))],na.rm=T),
                       EndingPerim=sum(PatchPerimLU[unique(as.vector(T1_PatchName))],na.rm=T),
                       StartingSummedPatchDiam=sum(PatchDiamLU[unique(as.vector(T0_PatchName))],na.rm=T),
                       EndingSummedPatchDiam=sum(PatchDiamLU[unique(as.vector(T1_PatchName))],na.rm=T),
                       TransitionMagnitude=EndingSize-StartingSize,
                       TransitionType=paste0(unique(TransitionType),collapse="_"),
                       PercentChange=(EndingSize-StartingSize)/StartingSize,
                       Log2Ratio_Change=log2(EndingSize/StartingSize))
VitalRate_Colony$TransitionTypeSimple="GROWTH"
VitalRate_Colony$TransitionTypeSimple[VitalRate_Colony$StartingSize==0]="RECR"
VitalRate_Colony$TransitionTypeSimple[VitalRate_Colony$EndingSize==0]="MORT"
VitalRate_Colony$EndingPerim[VitalRate_Colony$TransitionTypeSimple=="MORT"]=0
VitalRate_Colony$EndingSummedPatchDiam[VitalRate_Colony$TransitionTypeSimple=="MORT"]=0
VitalRate_Colony$StartingPerim[VitalRate_Colony$TransitionTypeSimple=="RECR"]=0
VitalRate_Colony$StartingSummedPatchDiam[VitalRate_Colony$TransitionTypeSimple=="RECR"]=0
  
VitalRate_Colony$TransitionTypeSimple[VitalRate_Colony$TransitionTypeSimple=="GROWTH"&VitalRate_Colony$EndingSize<VitalRate_Colony$StartingSize]="SHRINK"
VitalRate_Colony$Fragmented=VitalRate_Colony$N_t0>1|VitalRate_Colony$N_t1>1

table(VitalRate_Colony$Fragmented,VitalRate_Colony$TransitionTypeSimple)
VitalRate_Colony[VitalRate_Colony$Fragmented==T&VitalRate_Colony$TransitionTypeSimple=="RECR",]


PlotGenusLU=c("Montipora sp.","Pocillopora sp.","Porites sp.");names(PlotGenusLU)=c("MOSP","POCS","POSP")
VitalRate_Colony$Genus=PlotGenusLU[as.vector(VitalRate_Colony$Genus_Code)]
VitalRate_Colony$Recruit=as.numeric(VitalRate_Colony$TransitionType=="RECR")
VitalRate_Colony$Mortality=as.numeric(VitalRate_Colony$TransitionType=="MORT")
VitalRate_Colony$Fragmentation=as.numeric(VitalRate_Colony$Fragmented)
VitalRate_Colony$ln_SS=log(VitalRate_Colony$StartingSize)
VitalRate_Colony$ln_ES=log(VitalRate_Colony$EndingSize)
VRCd=subset(VitalRate_Colony,DataOrError=="DATA")
VRCd$Genus=factor(VRCd$Genus)
VRCd$StartingYear=year(VRCd$StartingDate)
VRCd$EndingYear=year(VRCd$EndingDate)
VRCd$Interval=paste0(substr(VRCd$StartingYear,3,4),"-",substr(VRCd$EndingYear,3,4))
#VRCd$SiteInterval=paste0(substr(VRCd$Site,5,11),"\n",VRCd$Interval) #old line
VRCd$SiteInterval=paste0(substr(VRCd$Site,1,11),"\n",VRCd$Interval) #Changed SiteInterval to include island name in the site

VRCd$Island=factor(substr(VRCd$Site,1,3),levels=c("FFS","HAW","MAI","OAH","PHR","KUR"))
VRCd$PropMagnitude=VRCd$EndingSize/VRCd$StartingSize
VRCd$AnnualPropRate_E=(VRCd$PropMagnitude)^(1/VRCd$Interval_Years)
VRCd$TransitionRate_L=VRCd$TransitionMagnitude/VRCd$Interval_Years
VRCd$AnnualEndingSize_E=VRCd$StartingSize*VRCd$AnnualPropRate_E
VRCd$TriennialEndingSize_E=VRCd$StartingSize*(VRCd$AnnualPropRate_E^3)
VRCd$TriennialEndingSize_E[VRCd$TransitionTypeSimple=="MORT"]=0
VRCd$ln_AES=log(VRCd$AnnualEndingSize_E)
VRCd$ln_3ES=log(VRCd$TriennialEndingSize_E)
VRCd_gs=subset(VRCd,TransitionTypeSimple%in%c("GROWTH","SHRINK"))

PatchLevel=VitalRate_Patch
ColonyLevel=VRCd
save(list = c("PatchLevel","ColonyLevel"),file = "1-Annotations_ColonyTransitions/Output/Patch_And_Colony_Data.rdata")

##########################################
