rm(list=ls())
# Loading Libraries -------------------------------------------------------
library(tidyverse)
library(lubridate)
library(stringr) 
library(pROC)
library(ggpubr)
library(sf)
library(rgeos)
library(sp)
library(patchwork)


# Loading Functions -------------------------------------------------------
leadz=function(x,n){return(formatC(x,width=n,flag=0))}
A2D=function(A){return(2*sqrt(A/pi))}
D2A=function(D){return((D/2)^2*pi)}
source("2-VitalRateFunctions/Scripts/gcdist.R")

# Loading / Managing DataFrames: ColonyLevel, surv_dat  ----------------------------------------------  ---------

##Manage ColonyLevel and create surv_dat
load("1-Annotations_ColonyTransitions/Output/Patch_And_Colony_Data.rdata") #ColonyLevel dataset 
#(This includes the PatchLevel dataset. We performed all annotations at the colony scale (the sum of all associated patches in cases of fission/fusion/partial mort))
if(!all(c("Survival","log10_ESc")%in%names(ColonyLevel))){
  #LN size
  ColonyLevel$ln_SS <- log(ColonyLevel$StartingSize)
  ColonyLevel$ln_ES <- log(ColonyLevel$EndingSize)
  #Log10 size
  ColonyLevel$log10_SS <- log10(ColonyLevel$StartingSize)
  ColonyLevel$log10_ES <- log10(ColonyLevel$EndingSize)
  #changing ln_ES/ln_SS from Na/-Inf to an actual value for recruitment/mortality
  ColonyLevel$log10_SS <- as.numeric(ColonyLevel$log10_SS)
  ColonyLevel$log10_SS[ColonyLevel$log10_SS == -Inf] <-NA
  
  #Mortality to Survival
  ColonyLevel <- cbind(ColonyLevel, data.frame(Survival = 1 - ColonyLevel$Mortality))
  
  # Making the size term yearly
  GY <- (ColonyLevel$ln_ES - ColonyLevel$ln_SS)*(1/ColonyLevel$Interval_Year)
  ColonyLevel <- cbind(ColonyLevel, data.frame(ln_ESc = ColonyLevel$ln_SS + GY))
  GY <- (ColonyLevel$log10_ES - ColonyLevel$log10_SS)*(1/ColonyLevel$Interval_Year)
  ColonyLevel <- cbind(ColonyLevel, data.frame(log10_ESc = ColonyLevel$log10_SS + GY))
  ColonyLevel$log10_ESc[ColonyLevel$log10_ESc == -Inf] <-NA
}
#Add island code to Lisianski
ColonyLevel$Island <- paste0(substr(ColonyLevel$Site,1,3))
#add region to ColonyLevel datset
RegionLU=c("NWHI","MHI","NWHI","NWHI","MHI","MHI","NWHI") #create a lookup table
names(RegionLU)=c("FFS","HAW","KUR","LIS","MAI","OAH","PHR")
ColonyLevel$REGION = RegionLU[ColonyLevel$Island]


#Prepping surv_dat output
surv_dat=ColonyLevel[,c("ColonyID","log10_SS","Survival","Genus_Code","Interval","SiteInterval","Site","N_t0","TransitionType","Interval_Years","StartingDate","EndingDate")]
surv_dat=subset(surv_dat,TransitionType!="RECR")
names(surv_dat)=c("ColonyID","size","survival","Genus_Code","Interval","SiteInterval","Site","N_t0","TransitionType","Interval_Years")

# Change these file names when saving data to refect the new changes
 save(ColonyLevel, file="2-VitalRateFunctions/Output/Colony_Data_edited.rdata")
 save(surv_dat, file="2-VitalRateFunctions/Output/Colony_Data_edited_survival.rdata")

# Loading DataFrames: ColonyLevel, surv_dat  ----------------------------------------------  ---------
#load("2-VitalRateFunctions/Output/Colony_Data_edited.rdata") #ColonyLevel
#load("2-VitalRateFunctions/Output/Colony_Data_edited_survival.rdata") #surv_dat

# Where we run the models -------------------------------------------------

#Define Site - Interval - Taxa groupings for which to run models
#SiteIntervalGenus
ColonyLevel$SIG=paste0(substr(ColonyLevel$Site,1,3),substr(ColonyLevel$Site,9,11),"_",
                       substr(year(ColonyLevel$StartingDate),3,4),leadz(month(ColonyLevel$StartingDate),2),"-",
                       substr(year(ColonyLevel$EndingDate),3,4),leadz(month(ColonyLevel$EndingDate),2),"_",
                       ColonyLevel$Genus_Code)
#SiteIntervalSpecies
ColonyLevel$SIS=paste0(substr(ColonyLevel$Site,1,3),substr(ColonyLevel$Site,9,11),"_",
                       substr(year(ColonyLevel$StartingDate),3,4),leadz(month(ColonyLevel$StartingDate),2),"-",
                       substr(year(ColonyLevel$EndingDate),3,4),leadz(month(ColonyLevel$EndingDate),2),"_",
                       ColonyLevel$Spec_Code)

#Report on Grouping
Usig=unique(ColonyLevel$SIG)
Usis=unique(ColonyLevel$SIS)
Nsig=length(Usig);Nsig
Nsis=length(Usis);Nsis

Name="VR_Models_"

# Load Data to Get Proportion of "True" Recruitment ------------------------------------
#Assign SFM Sites to Sectors 
sitemd=read.csv("2-VitalRateFunctions/Data/Sites_Metadata.csv");names(sitemd)[1]="Site";names(sitemd)[1]="Site";sitemd$EndingDate=mdy(sitemd$SampleDate.MM.DD.YYYY.)
sitemd <- select(sitemd, Site, Latitude, Longitude) #get 1 lat/long for each site
sitemd <- distinct(sitemd)
#sitemd=sitemd %>% group_by(Site) %>% summarize(Lat=mean(Latitude),Lon=mean(Longitude))

#Reassign OAH_XX_022 to OAH_OCC_005
ColonyLevel$Site[which(ColonyLevel$Site=="OAH_XXX_022")]="OAH_OCC_005"

uSD=ColonyLevel[,c("Site","StartingDate","EndingDate")] %>% 
  pivot_longer(cols=c("StartingDate","EndingDate"),values_to="Date") %>% 
  select(all_of(c("Site","Date"))) %>% 
  unique()

CL.sf=as.data.frame(left_join(uSD,sitemd[,c("Site","Latitude","Longitude")],by=c("Site"))) #Site,Date,Lat,Long dataframe
CL.sf=st_as_sf(na.omit(CL.sf),coords=c("Longitude","Latitude")) #create geometry column
secshp= st_read(dsn = "2-VitalRateFunctions/Data/SectorSHP",layer = "ALLPacific_Sectors_Islands_5km_buffer")

CL.sf = CL.sf %>% st_set_crs(value = st_crs(secshp)) #retrieve coordinate system
Site2Sec=st_join(CL.sf,secshp[,"SEC_NAME"])

#Add Sector to ColonyLevel dataframe!
Site2Sec_ONLY=unique(st_drop_geometry(Site2Sec)[,c("Site","SEC_NAME")])
Site2Sec_ONLY$SEC_NAME[Site2Sec_ONLY$Site=="HAW_OCC_003"]="HAW_PUNA"
Site2Sec_ONLY$SEC_NAME[Site2Sec_ONLY$Site=="HAW_SIO_K08"]="HAW_KONA"
ColonyLevel=ColonyLevel %>% left_join(Site2Sec_ONLY[,c("Site","SEC_NAME")],by="Site")
#Save final dataset of coral colonies
save(ColonyLevel, file="2-VitalRateFunctions/Output/Colony_Data_edited.rdata")


#Get observed proportion of 'juveniles' as true recruits (SfM data)
#gray=juveniles, blue=true recruits
RECs=subset(ColonyLevel,TransitionTypeSimple=="RECR")
bw=0.5
ggplot()+
  geom_histogram(data=subset(ColonyLevel,EndingSize>0),aes(A2D(EndingSize)),fill="grey",binwidth=bw, labels= "Juveniles")+
  xlim(c(0,10))+
  geom_histogram(data=RECs,aes(A2D(EndingSize)),fill="blue",binwidth=bw, labels= "True Recruits")+
  facet_grid(Genus_Code~.,scales = "free_y")+
  geom_vline(xintercept = 5)+
  theme_bw()+
  #ggtitle("Proportion of juveniles compared to true recruits")+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), strip.text.y= element_text(size = 12))+
  labs(x="Diameter (cm)",color = "Legend")

#Check Proportion of Recruits per Taxon - hist calculation
binwidths=bw
uG=c("SSSS","POCS","POSP","MOSP")
PropRec_g=cbind(expand.grid(Genus_Code=uG,BinWidth=binwidths),PropRec=NA)
for(i_b in 1:length(binwidths)){
  binwidth=binwidths[i_b]
  N=5/binwidth
  for(i_G in 1:length(uG)){
    if(uG[i_G]=="SSSS"){
      hcl=hist(A2D(subset(ColonyLevel,EndingSize>0)$EndingSize),breaks=seq(0,1000,by=binwidth),plot = F)
      hrc=hist(A2D(RECs$EndingSize),breaks=seq(0,1000,by=binwidth),plot = F)
      PropRec_v=(hrc$counts[1:(N+1)]/hcl$counts[1:(N+1)])
      out_i=which(PropRec_g$Genus_Code==uG[i_G]&PropRec_g$BinWidth==binwidth)
      PropRec_g$PropRec[out_i]=mean(PropRec_v[1:(N+1)],na.rm=T)
    }else{
      hcl=hist(A2D(subset(ColonyLevel,EndingSize>0&Genus_Code==uG[i_G])$EndingSize),breaks=seq(0,1000,by=binwidth),plot = F)
      hrc=hist(A2D(subset(RECs,Genus_Code==uG[i_G])$EndingSize),breaks=seq(0,1000,by=binwidth),plot = F)
      PropRec_v=(hrc$counts[1:(N+1)]/hcl$counts[1:(N+1)])
      out_i=which(PropRec_g$Genus_Code==uG[i_G]&PropRec_g$BinWidth==binwidth)
      PropRec_g$PropRec[out_i]=mean(PropRec_v[1:(N+1)],na.rm=T)
      # plot(hcl$breaks[1:(N+1)],PropRec_v,type="b",ylim=c(0,1),main=uG[i_G])
      # abline(h=mean(propREC,na.rm=T))}
    }
  }
}
  PropRec_g
  
  detach(package:plyr)
  PropRecMn=PropRec_g %>% 
    group_by(Genus_Code) %>% 
    summarise(meanPropRec=mean(PropRec))
  #Final 0-5 cm diam. Proportion Recruits in the Juvenile size classes (will use for 'pro-rating')
  PropRecMn
  
  
  
  ################################################################################
  ######### Site + Interval Years + Genus models #################################
  ########## Recruitment ! #######################################################
  ################################################################################
  
  ####Case 1 Assume Site - Level Stock Recruitment Relationship
  ####### Unlikely to hold true, but might be a decent approximation in skewed dispersal kernels
  #recval is N recruits per Area of adult biomass - N/cm^2
  # EQN 1: Then recval = ((Observed recruits-N)/ (Area Surveyed-cm^2)) / (Adult Area-cm^2)/Area Surveyed-cm^2)
  # Adult area could come from site-level percent cover*area surveyed; or measured SFM adult area / Area surveyed
  
  ####Case 2 Assume Sector - Level Stock Recruitment Relationship
  
  #recval is N recruits per Area of adult biomass - N/cm^2
  # EQN 2: Then sector SFM recval = (SFM Observed recruits-N)/ Area Surveyed-cm^2) / (sector-level proportional cover) ###?
  # EQN 3: Or sector REA recval = (REA Observed recruits***-N)/ Area Surveyed-cm^2) / (sector-level proportional cover)
  # *** here we should prorate our REA juv density by SFM calibrated prop. of recruits in juv size classes
  
  
  ################################################################################
  #SFM Observed Recruitment Rates 
  ################################################################################
  
  #Case 1: Assume Site-Level Stock-Recruitment
  
  #Area Surveyed for each SIG
  Asurv=read.csv("2-VitalRateFunctions/Data/AreaSurveyed_N_Circrats.csv"); names(Asurv)[1]="Site"
  Asurv$Site[Asurv$Site=="OAH_XXX_022"]="OAH_OCC_005"
  Asurv$EndingDate=mdy(Asurv$Date)
  Asurv_l=Asurv[,c("Site","EndingDate","POCS","POSP","MOSP")]%>%
    pivot_longer(cols=c("POSP","MOSP","POCS"),names_to=c("Genus_Code"),values_to=c("Ncircrats")) %>% 
    mutate(Ncircrats=as.numeric(Ncircrats))
  Asurv_l$A.Surv.m2 =Asurv_l$Ncircrats*0.5 #numb of circrats * 0.5m2 (size of circrat) to get area surveyed in m2
  Asurv_l$A.Surv.cm2 = Asurv_l$Ncircrats*5000 #area surveyed cm2
  
  #Area of Adult Colonies for each SIG
  Atax=ColonyLevel %>%
    group_by(SIG,Site,Interval,Genus_Code,StartingDate,EndingDate,Interval_Years) %>% 
    summarise(AdultCoralArea_cm2=sum(StartingSize)) #Area of adult cm2 = sum all corals in each genus for each year
  
  #Link to Nrecruits from the StartingSize summed area: Adults at beginning of interval generate 
  RecSFMTib=ColonyLevel %>%
    filter(TransitionTypeSimple %in% c("RECR")) %>% 
    group_by(SEC_NAME,SIG,Site,Interval,Genus_Code,REGION,StartingDate,EndingDate,Interval_Years) %>% 
    summarize(Nrec=length(which(TransitionTypeSimple=="RECR"))) %>% 
    left_join(Atax) %>% 
    left_join(Asurv_l) %>% 
    mutate(
      # numb recruits/area surveyed cm2
      RecSFM_p_Survcm2=Nrec/A.Surv.cm2,
      #Case 1: Assume Site-Level Stock-Recruitment
      #Recruit N per adult area and annual rate
      RecSFM_p_SiteAdcm2=Nrec/AdultCoralArea_cm2,
      RecSFM_p_SiteAdcm2_Yr=RecSFM_p_SiteAdcm2/Interval_Years) # numb recruits/total adult area cm2 annualized
  
  RecSFMDataFrame=as.data.frame(RecSFMTib)
  write.csv(RecSFMDataFrame, "2-VitalRateFunctions/Output/SiteLevelStockRec.csv", row.names = FALSE)
  save(RecSFMDataFrame, file = paste0("2-VitalRateFunctions/Output/",Name,"_RecSFMrates.rdata"))
  
  #plot
  hist(RecSFMDataFrame$RecSFM_p_SiteAdcm2_Yr, breaks = 20)
  
  
  #REGIONAL CALCULATIONS BY GENUS AND REGION
  RecSFMDataFrame$RecSFM_p_SiteAdcm2_Yr[RecSFMDataFrame$RecSFM_p_SiteAdcm2_Yr == Inf] <-NA
  #Calculate median, mean, and lower/upper quantiles for site-level stock recruitment GENUS values
  RecSFMSummary=RecSFMDataFrame %>%
    dplyr::select(SEC_NAME,Site,Interval,Genus_Code,REGION,RecSFM_p_SiteAdcm2_Yr)%>%
    group_by(Genus_Code,REGION)%>%
    summarise(SiteStock_median = median(RecSFM_p_SiteAdcm2_Yr, na.rm = TRUE),
              SiteStock_mean = mean(RecSFM_p_SiteAdcm2_Yr, na.rm = TRUE),
              SiteStock_q05= quantile(RecSFM_p_SiteAdcm2_Yr, c(0.05),na.rm = TRUE),
              SiteStock_q95= quantile(RecSFM_p_SiteAdcm2_Yr, c(0.95),na.rm = TRUE)
    )
  save(RecSFMSummary, file = paste0("2-VitalRateFunctions/Output/",Name,"_RecSFMrates_regional.rdata"))
  write.csv(RecSFMSummary, "2-VitalRateFunctions/Output/SiteLevelStockRec_regional.csv", row.names = FALSE)
  

  
  
  ################################################################################
  #REA Recruitment Rates 
  ################################################################################
  #REA Observed Recruitment Rates per adult cover, at Sector and Site Level, with SIG matching site-level data grouped at sector containing SIG
  
  #Case 2: Assume Sector-Level Stock-Recruitment
  #Site Recruitment, Sector Stock
  
  #Get Site-Level and Sector_level REA juv density
  rea_sec=read.csv("2-VitalRateFunctions/Data/NWHI_MHI_REA_Data_Sector.csv")
  
  #Get Percent Cover by Genus at each Sector
  cov1_sec=read.csv("2-VitalRateFunctions/Data/NWHI_MHI_Cover_T1_Data_Sector.csv")
  cov3_sec=read.csv("2-VitalRateFunctions/Data/NWHI_MHI_Cover_T3_Data_Sector.csv")
  uSec=na.omit(unique(ColonyLevel$SEC_NAME))
  names(cov1_sec)[2:15]=substr(names(cov1_sec)[2:15],6,999)
  meta=names(cov1_sec)[2:6]
  classes1=names(cov1_sec)[7:15]
  SEclasses1=names(cov1_sec)[22:30]
  
  names(cov3_sec)[2:111]=substr(names(cov3_sec)[2:111],6,999)
  classes3=names(cov3_sec)[7:111]
  SEclasses3=names(cov3_sec)[118:(104+118)]
  
  #select first 15 cols of cov1_sec and pivot CORAL,MA,TURF,CCA,EMA,SC,I,SED,HAL to
  #a new column titled CAT and values to column titled COVER
  cov1_lse = cov1_sec[,c(meta,SEclasses1)] %>% 
    pivot_longer(cols = all_of(SEclasses1),names_to="CAT",values_to="SE.COVER")
  cov1_lse$CAT=substr(cov1_lse$CAT,10,999) #drop pooled SE
  cov1_l = cov1_sec[,c(meta,classes1)] %>% #get cover by category & N
    pivot_longer(cols = all_of(classes1),names_to="CAT",values_to="COVER") %>% 
    left_join(cov1_lse,by=c(meta,"CAT"))
  cov1_l=cov1_l[,c("REGION","ISLAND","ANALYSIS_SEC","OBS_YEAR","CAT","COVER","SE.COVER","N")]
  cov1_l$SD.COVER=cov1_l$SE.COVER*sqrt(cov1_l$N)
  cov1_l$VAR.COVER=cov1_l$SD.COVER^2
  
  cov3_lse = cov3_sec[,c(meta,SEclasses3)] %>% 
    pivot_longer(cols = all_of(SEclasses3),names_to="CAT",values_to="SE.COVER")
  cov3_lse$CAT=substr(cov3_lse$CAT,10,999)
  cov3_l = cov3_sec[,c(meta,classes3)] %>% 
    pivot_longer(cols = all_of(classes3),names_to="CAT",values_to="COVER") %>% 
    left_join(cov3_lse,by=c(meta,"CAT"))
  cov3_l=cov3_l[,c("REGION","ISLAND","ANALYSIS_SEC","OBS_YEAR","CAT","COVER","SE.COVER","N")]
  cov3_l$SD.COVER=cov3_l$SE.COVER*sqrt(cov3_l$N)
  cov3_l$VAR.COVER=cov3_l$SD.COVER^2
  
  
  MPPcovers=c("MOBR","MOEN","MONE","MOFO","POCS","POEN","POFO","POMA","PONM","POBR") #cover categories for each genus for Tier3 benthic cat
  names(MPPcovers)=c(rep("MOSP",4),"POCS",rep("POSP",5))
  MMPlu=names(MPPcovers)
  names(MMPlu) = MPPcovers
  
  cov_l=rbind(subset(cov1_l,CAT=="CORAL"),subset(cov3_l,CAT%in%MPPcovers)) 
  cov_l$GENUS_CODE=MMPlu[cov_l$CAT] #add genus_code column
  cov_l$GENUS_CODE[cov_l$CAT=="CORAL"]="SSSS" #if category is coral, assign SSSS genus code
  cov_l$SEC_YEAR=paste0(cov_l$ANALYSIS_SEC,"-",cov_l$OBS_YEAR)
  cov_l$REGION=factor(cov_l$REGION,levels=c("MHI","NWHI","PRIAs","MARIAN","SAMOA"))
  #SUM all within genus variation at each site
  keepcols=c("REGION","ISLAND","ANALYSIS_SEC","OBS_YEAR","GENUS_CODE","SEC_YEAR")
  cov_sum=cov_l %>%  #cover by genus code and sector
    group_by(ANALYSIS_SEC,OBS_YEAR,GENUS_CODE) %>% 
    summarize(COVER=sum(COVER),VAR.COVER=sum(VAR.COVER),N.COVER=sum(N))
  cov_sum$SD.COVER=sqrt(cov_sum$VAR.COVER)
  cov_sum$SE.COVER=cov_sum$SD.COVER/sqrt(cov_sum$N.COVER)
  
  cov=left_join(cov_sum,unique(cov_l[,keepcols]),by=c("ANALYSIS_SEC","OBS_YEAR","GENUS_CODE")) #cov_l by distint benth cats. sum all cover by genus
  cov=cov[,c(keepcols,"COVER","SE.COVER","N.COVER")] #reorganize columns
  #cover is percent cover (0-100%)
  
  
  # prorated mean juvenile col density for each sector and each year divided by adult cover sector year genus
  #get mean, 5 and 95 quantiles for each sector year genus
  #for each of those ^^ 3 numbers, want to prorate depending on genus code. 
  
  table(cov$OBS_YEAR,useNA = "always")
  table(rea_sec$ANALYSIS_YEAR,useNA = "always")
  names(rea_sec)[which(names(rea_sec)=="Sector")]="ANALYSIS_SEC"
  names(rea_sec)[which(names(rea_sec)=="n")]="N.JCD"
  rea_sec$OBS_YEAR=rea_sec$ANALYSIS_YEAR
  
  # Propagating error for sector level stock recruitment
  #mJCDp =  mJCD*MeanPropRec
  #seJCDp = seJCD*MeanPropRec
  #sdJCDp = sqrt(nJCD)*seJCDp
  #nJCDp = nJCD
  #sdAC = sqrt(nAC)*seAC
  #mRec = mJCDp/mAC
  #sdRec = mRec*sqrt((sdJCDp/mJCDp)^2 + (sdAC/mAC)^2))
  #nRec = min(nJCD,nAC)????
  #seRec = sdRec/sqrt(nRec)
  #CI95Rec = 1.96*seRec
  
  cov$P_COVER=cov$COVER/100 #convert percent cover to proportional cover
  cov$SE.P_COVER=cov$SE.COVER/100
  rea_sec$Mean_JuvColDen_cm2=  rea_sec$Mean_JuvColDen/10^4 #convert mean juv col density to cm2
  rea_sec$SE_JuvColDen_cm2=  rea_sec$SE_JuvColDen/10^4 #convert juv col den SE to cm2

  #Calculate sector level stock recruitment and propagate error
  #sector level recruitment for sectors for this study
  Jsec_P_Asec= rea_sec %>%
    filter(ANALYSIS_SEC%in%uSec&GENUS_CODE%in%c("SSSS","MOSP","POSP","POCS")) %>% #filter by target taxa and target sectors (uSec is unique sector list from above)
    select(REGION,ISLAND,ANALYSIS_SEC,OBS_YEAR,GENUS_CODE,Mean_JuvColDen_cm2, SE_JuvColDen_cm2,N.JCD)%>%
    left_join(cov,by=c("REGION","ISLAND","ANALYSIS_SEC","GENUS_CODE","OBS_YEAR")) %>% 
    left_join(PropRecMn,by=c("GENUS_CODE"="Genus_Code")) %>%
    group_by(REGION,ISLAND,ANALYSIS_SEC,OBS_YEAR,GENUS_CODE)%>%
    mutate(Mean_JuvColDenP=Mean_JuvColDen_cm2*meanPropRec,
           SE_JuvColDenP=SE_JuvColDen_cm2*meanPropRec,
           N.JCDp=N.JCD,
           SD_JuvColDenP=sqrt(N.JCDp)*SE_JuvColDenP,
           SD.P_COVER=sqrt(N.COVER)*SE.P_COVER,
           MN.RecVal = Mean_JuvColDenP/P_COVER,
           SE.RecVal = MN.RecVal*sqrt((SE_JuvColDenP/Mean_JuvColDenP)^2+(SE.P_COVER/P_COVER)^2),#propogating SE error
           N.RecVal = min(N.JCDp,N.COVER,na.rm=T),
           #SE.RecVal = SD.RecVal/sqrt(N.RecVal),
           CI95.RecVal = 1.96*SE.RecVal,
           LOW_CI95.RecVal= MN.RecVal - CI95.RecVal,
           HIGH_CI95.RecVal= MN.RecVal + CI95.RecVal
           )
  write.csv(Jsec_P_Asec, "2-VitalRateFunctions/Output/SectorLevelStockRec.csv", row.names = FALSE)
  
  #plot site-level stock recruitment values and sector-level stock recruitment values
  Site=ggplot(RecSFMDataFrame,aes(RecSFM_p_SiteAdcm2_Yr))+
    geom_histogram()+ theme_bw()+
    facet_grid("Genus_Code")+
    #facet_wrap(Genus_Code~REGION)+
    xlim(c(0,.12)) + xlab("Site-level stock recruitment") + ylab("# recruits/cm^2 area adult coral")
  
  Sec= Jsec_P_Asec%>%
    filter(GENUS_CODE%in%c("MOSP","POSP","POCS")) %>%
    ggplot(aes(MN.RecVal))+
    geom_histogram()+theme_bw()+
    facet_grid("GENUS_CODE")+
    #facet_wrap(GENUS_CODE~REGION)+
    xlim(c(0,.12))+ xlab("Proportional sector-level stock recruitment") +ylab("# juveniles/cm^2 area adult cover") 
  
  Site/Sec
  
  
  
  #Sector Stock Rec for REGIONAL MODEL
  Regional_SectorStockRec <- as.data.frame(Jsec_P_Asec) 
  Regional_SectorStockRec$MN.RecVal[Regional_SectorStockRec$MN.RecVal == Inf] <- NA
  Regional_SectorStockRec$SE.RecVal[Regional_SectorStockRec$SE.RecVal == Inf] <- NA
  
  Regional_SectorStockRec <- Regional_SectorStockRec %>%
    select(ANALYSIS_SEC,OBS_YEAR,GENUS_CODE,REGION,MN.RecVal,SE.RecVal)%>%
    group_by(GENUS_CODE, REGION) %>%
    summarise(SectorStock_median = median(MN.RecVal,na.rm = TRUE),
              SectorStock_MD_CI95_LO = median(MN.RecVal-1.96*SE.RecVal,na.rm=T),
              SectorStock_MD_CI95_HI = median(MN.RecVal+1.96*SE.RecVal,na.rm=T),
              SectorStock_mean = mean(MN.RecVal,na.rm = TRUE),
              SectorStock_MN_CI95_LO = mean(MN.RecVal-1.96*SE.RecVal, na.rm = TRUE), 
              SectorStock_MN_CI95_HI = mean(MN.RecVal+1.96*SE.RecVal, na.rm = TRUE) 
    )
  #change negative low CI values to 0 since you can't have negative recruitment
  Regional_SectorStockRec[Regional_SectorStockRec<0] = 0
  
  
  
  ##############################
  #####SECTOR-STOCK-ALL REC#####
  #Case 2 (continued): Assume Sector-Level Stock-Recruitment using ALL sectors, not only study sectors)
  #############################
  
  #sector level recruitment for all sectors in the NWHI and MHI (not just study sectors)
  Jsec_P_ALL= rea_sec %>%
    filter(REGION%in%c("NWHI","MHI")&GENUS_CODE%in%c("SSSS","MOSP","POSP","POCS")) %>% #uSec is unique sector list from above
    select(REGION,ISLAND,ANALYSIS_SEC,OBS_YEAR,GENUS_CODE,Mean_JuvColDen_cm2, SE_JuvColDen_cm2,N.JCD)%>%
    left_join(cov,by=c("REGION","ISLAND","ANALYSIS_SEC","GENUS_CODE","OBS_YEAR")) %>% 
    left_join(PropRecMn,by=c("GENUS_CODE"="Genus_Code")) %>%
    group_by(REGION,ISLAND,ANALYSIS_SEC,OBS_YEAR,GENUS_CODE)%>%
    mutate(Mean_JuvColDenP=Mean_JuvColDen_cm2*meanPropRec,
           SE_JuvColDenP=SE_JuvColDen_cm2*meanPropRec,
           N.JCDp=N.JCD,
           SD_JuvColDenP=sqrt(N.JCDp)*SE_JuvColDenP,
           SD.P_COVER=sqrt(N.COVER)*SE.P_COVER,
           MN.RecVal = Mean_JuvColDenP/P_COVER,
           SE.RecVal = MN.RecVal*sqrt((SE_JuvColDenP/Mean_JuvColDenP)^2+(SE.P_COVER/P_COVER)^2),#propogating SE error
           N.RecVal = min(N.JCDp,N.COVER,na.rm=T),
           #SE.RecVal = SD.RecVal/sqrt(N.RecVal),
           CI95.RecVal = 1.96*SE.RecVal,
           LOW_CI95.RecVal= MN.RecVal - CI95.RecVal,
           HIGH_CI95.RecVal= MN.RecVal + CI95.RecVal
    )
  
  #plot the sector-level-all rec value for all islands in the Hawaiian Archipelago
  Jsec_P_ALL %>% #filter(ISLAND%in%c("Hawaii","Maui","Oahu","French Frigate","Lisianski","Kure")) %>% 
  ggplot(aes(MN.RecVal))+geom_histogram(binwidth=0.001)+facet_grid(ISLAND~GENUS_CODE)#+xlim(c(0,.12))
    

  #Rec values by region and genus code (SECTOR-STOCK-ALL REC) for REGIONAL MODEL
  RecVal_Sec_Dists=Jsec_P_ALL %>% filter(!is.infinite(MN.RecVal)) %>% 
    group_by(REGION,GENUS_CODE) %>% 
    summarize(MD.RecVal_Sec_All=median(MN.RecVal,na.rm=T),
              MD.CI95_LO.RecVal_Sec_All=median(MN.RecVal-1.96*SE.RecVal,na.rm=T),
              MD.CI95_HI.RecVal_Sec_All=median(MN.RecVal+1.96*SE.RecVal,na.rm=T),
              MN.RecVal_Sec_All=mean(MN.RecVal,na.rm=T),
              MN.CI95_LO.RecVal_Sec_All=mean(MN.RecVal-1.96*SE.RecVal,na.rm=T),
              MN.CI95_HI.RecVal_Sec_All=mean(MN.RecVal+1.96*SE.RecVal,na.rm=T)
    )

#plot distributions of sector-stock-all recruitment values for REGIONAL MODEL
  Jsec_P_ALL %>% #filter(ISLAND%in%c("Hawaii","Maui","Oahu","French Frigate","Lisianski","Kure")) %>% 
    ggplot(aes(MN.RecVal))+
    geom_histogram(binwidth=0.01)+
    geom_vline(aes(xintercept=MD.RecVal_Sec_All),data=RecVal_Sec_Dists,color="darkcyan",lty=1)+
    geom_vline(aes(xintercept=MD.CI95_LO.RecVal_Sec_All),data=RecVal_Sec_Dists,color="darkcyan",lty=3)+
    geom_vline(aes(xintercept=MD.CI95_HI.RecVal_Sec_All),data=RecVal_Sec_Dists,color="darkcyan",lty=3)+
    facet_grid(REGION~GENUS_CODE)+
    scale_x_sqrt()  
  
  
  ############################
  #REGIONAL MODEL INTEGRATION: Sector-stock and sector-stock-all recruitment
  
  #Include Sec_All values in Regional sector stock
  Regional_SectorStockRec <- left_join(Regional_SectorStockRec,RecVal_Sec_Dists, by= c("GENUS_CODE","REGION"))
  
  write.csv(Regional_SectorStockRec, "2-VitalRateFunctions/Output/SectorLevelStockRec_regional.csv", row.names = FALSE)
  
  
  
  
  ################################################################################
  ######### Site + Interval Years + Genus models ###############################
  ########## GROWTH ! ############################################################
  ################################################################################
  
  GrowthTib=ColonyLevel %>% 
    filter(TransitionTypeSimple %in% c("GROWTH","SHRINK")) %>% 
    group_by(SIG,Site,Interval,Genus_Code,StartingDate,EndingDate,Interval_Years) %>% 
    nest() %>% 
    mutate(
      GrowthMod=map(data,~lm(log10_ESc ~ log10_SS,data=as.data.frame(data))),
      g.N=map(data,~length(.$log10_SS)),
      g.int=map(GrowthMod,~coef(.)[1]),
      g.slp=map(GrowthMod,~coef(.)[2]),
      g.var=map(GrowthMod,~(summary(.)$sigma)^2),
      g.R2=map(GrowthMod,~summary(.)$adj.r.squared),
      g.AIC=map(GrowthMod,~AIC(.)),
    ) %>% 
    unnest(c(SIG, Site,Interval,Genus_Code,StartingDate,Interval_Years,EndingDate,g.N,g.int,g.slp,g.var,g.R2,g.AIC)) %>% 
    select_if(negate(is.list)) 
  
  GrowthDataFrame = as.data.frame(GrowthTib)
  save(GrowthDataFrame, file = paste0("2-VitalRateFunctions/Output/",Name,"_Gmodfits.rdata"))

  
  #Regional growth model
  GrowthRegionalTib=ColonyLevel %>% 
    filter(TransitionTypeSimple %in% c("GROWTH","SHRINK")) %>% 
    group_by(Genus_Code,REGION) %>% 
    nest() %>% 
    mutate(
      GrowthMod=map(data,~lm(log10_ESc ~ log10_SS,data=as.data.frame(data))),
      g.N=map(data,~length(.$log10_SS)),
      g.int=map(GrowthMod,~coef(.)[1]),
      g.slp=map(GrowthMod,~coef(.)[2]),
      g.var=map(GrowthMod,~(summary(.)$sigma)^2),
      g.R2=map(GrowthMod,~summary(.)$adj.r.squared),
      g.AIC=map(GrowthMod,~AIC(.)),
    ) %>% 
    unnest(c(Genus_Code,g.N,g.int,g.slp,g.var,g.R2,g.AIC)) %>% 
    select_if(negate(is.list)) 
  head(GrowthRegionalTib)
  GrowthRegionalDF = as.data.frame(GrowthRegionalTib)
  save(GrowthRegionalDF, file = paste0("2-VitalRateFunctions/Output/",Name,"_Gmodfits_regional.rdata"))
  
  
  ################################################################################
  ######### Genera + Site + Interval Years  models ###############################
  ## SURVIVAL ! ##################################################################
  ################################################################################
  #length(which(is.na(subset(ColonyLevel,TransitionTypeSimple %in% c("GROWTH","SHRINK","MORT"))$log10_SS)))
  
  #To annualize survival:
  #Fit logistic regression at Site-Interval-Genus scae. Use this model to predict estimated survival probability for each individual.
  #Annualize survival probability (raise to 1/IntervalYears) and add to ColonyLevel dataframe
  #Use annualized survival probability to refit logistic models for aggregate regional models 
  #Run SIG models like normal using logistic regression
  
  
  #Set up for pooling larger than SIG survival
  wo=getOption("warn");options(warn = -1)
  AnnSurvTib=ColonyLevel %>% 
    filter(TransitionTypeSimple %in% c("GROWTH","SHRINK","MORT")&!is.na(log10_SS)) %>% 
    group_by(SIG,Site,Interval,Genus_Code,StartingDate,EndingDate,Interval_Years) %>% 
    nest() %>% 
    mutate(
      s.N=map(data,~length(.$Survival)),
      SurvMod=map(data,~glm(Survival ~ log10_SS, family = "binomial" , data = as.data.frame(data))),
      SurvProbs=map(SurvMod,~predict(.,type = "response",na.action="na.pass"))#,
    ) %>% 
    ungroup() %>% 
    unnest(c(data,s.N,SurvProbs)) %>% 
    select_if(negate(is.list)) 
  AnnSurvTib$AnnSurvProbs=AnnSurvTib$SurvProbs^(1/AnnSurvTib$Interval_Years)
  options(warn = wo)
  
  #checking to make sure predict call works and we're getting actual annualized surv
  Usig=unique(AnnSurvTib$SIG)
  temp = subset(AnnSurvTib,SIG==Usig[12])
  ggplot(temp, aes(x=log10_SS))+ geom_point(aes(y=SurvProbs), color = "blue")+ geom_point(aes(y=AnnSurvProbs), color = "red")
  
  
  ColonyLevel_ap=left_join(ColonyLevel,AnnSurvTib[,c("SIG","ColonyID","s.N","AnnSurvProbs")])
  
  #plot annualized survival probability
  ggplot(ColonyLevel_ap,aes(x=log10_SS,y=AnnSurvProbs,fill=Site,size=s.N))+
    geom_point(shape=21,color="white")+
    facet_grid(Island~Genus_Code)
  
  ### JUST REGIONAL by Genus
  RegionalAnnSurvTib=ColonyLevel_ap %>% 
    filter(TransitionTypeSimple %in% c("GROWTH","SHRINK","MORT")) %>% 
    group_by(Genus_Code,REGION) %>% 
    nest() %>% 
    mutate(
      SurvMod=map(data,~glm(AnnSurvProbs ~ log10_SS, family = "binomial" , data = as.data.frame(data))),
      NullSurvMod=map(data,~glm(AnnSurvProbs ~ 1, family = "binomial" , data = as.data.frame(data))),
      s.N=map(data,~length(.$AnnSurvProbs)),
      s.int=map(SurvMod,~coef(.)[1]), 
      allsurv = map(data,~length(.$AnnSurvProbs)),
      s.int=map(SurvMod,~coef(.)[1]),
      s.slp=map(SurvMod,~coef(.)[2]),
      s.pR2=map(SurvMod,~1-(summary(.)$deviance/summary(.)$null.deviance)),
      s.AIC=map(SurvMod,~AIC(.))
    ) %>% 
    unnest(c(Genus_Code,s.N,s.int,s.slp, s.pR2,s.AIC)) %>% 
    select_if(negate(is.list)) 
  RegionalAnnSurvDataFrame = as.data.frame(RegionalAnnSurvTib)
  head(RegionalAnnSurvDataFrame)
  save(RegionalAnnSurvDataFrame, file = sprintf("2-VitalRateFunctions/Output/%s_Smodfits_regional.rdata", Name))
  
  # SIG Survival
  SurvTib=ColonyLevel %>% 
    filter(TransitionTypeSimple %in% c("GROWTH","SHRINK","MORT")) %>% 
    group_by(SIG,Site,Interval,Genus_Code,StartingDate,EndingDate,Interval_Years) %>% 
    nest() %>% 
    mutate(
      SurvMod=map(data,~glm(Survival ~ log10_SS, family = "binomial" , data = as.data.frame(data))),
      NullSurvMod=map(data,~glm(Survival ~ 1, family = "binomial" , data = as.data.frame(data))),
      s.N=map(data,~length(.$Survival)),
      s.int=map(SurvMod,~coef(.)[1]), 
      allsurv = map(data,~length(.$Survival)),
      s.int=map(SurvMod,~coef(.)[1]),
      s.slp=map(SurvMod,~coef(.)[2]),
      s.pR2=map(SurvMod,~1-(summary(.)$deviance/summary(.)$null.deviance)),
      s.AIC=map(SurvMod,~AIC(.))
    ) %>% 
    unnest(c(SIG,Site,Interval,Genus_Code,StartingDate,Interval_Years,EndingDate,s.N,s.int,s.slp, s.pR2,s.AIC)) %>% 
    select_if(negate(is.list)) 
  SurvDataFrame = as.data.frame(SurvTib)
  head(SurvDataFrame)
  save(SurvDataFrame, file = sprintf("2-VitalRateFunctions/Output/%s_Smodfits.rdata", Name))
  
  #compare regional vs SIG slope/int
  ggplot()+
    geom_point(aes(s.int,s.slp,size=s.N,shape=Genus_Code),color="red",data=RegionalAnnSurvDataFrame)+
    geom_point(aes(s.int,s.slp,size=s.N,shape=Genus_Code),data=subset(SurvDataFrame,s.N>=0))+
    facet_wrap("Genus_Code")
  
  ################################################################################
  ######### Site + Interval Years + Genus models #################################
  ########## INTEGRATION ! #######################################################
  ################################################################################
  #region lookup
  RegLU=c("NWHI","MHI","NWHI","NWHI","MHI","MHI","NWHI")
  names(RegLU)=c("FFS","HAW","KUR","LIS","MAI","OAH","PHR")
  
  #SIG models
  combo <- left_join(GrowthDataFrame,SurvDataFrame)
  #save(combo, file = sprintf("2-VitalRateFunctions/Data/%s_allmodfits.rdata",Name))
  no_NAs <- na.omit(combo)
  GrowthSurv_SIG <- data.frame()
  GrowthSurv_SIG <- no_NAs[order(no_NAs$Site),]
  GrowthSurv_SIG$StartingYear=as.numeric(year(GrowthSurv_SIG$StartingDate))
  GrowthSurv_SIG$SEC_NAME=Site2Sec_ONLY[match(GrowthSurv_SIG$Site,Site2Sec_ONLY$Site),"SEC_NAME"]
  
  #sort(table(GrowthSurv_SIG$Site,GrowthSurv_SIG$Genus_Code,GrowthSurv_SIG$StartingYear))
  
  #add SfM site-level stock recruitment data
  GrowthSurvRec_SIG = left_join(GrowthSurv_SIG,
                                RecSFMDataFrame[,c("Site","Genus_Code","EndingDate","Nrec","RecSFM_p_SiteAdcm2_Yr")],
                                by=c("Site","Genus_Code","StartingDate"="EndingDate"))

  #add REA sector-level stock recruitment data
  ModelParams_SIG= left_join(GrowthSurvRec_SIG,
                             Jsec_P_Asec[,c("ANALYSIS_SEC","GENUS_CODE","OBS_YEAR","N.RecVal","MN.RecVal","LOW_CI95.RecVal","HIGH_CI95.RecVal")],
                             by=c("SEC_NAME"="ANALYSIS_SEC","Genus_Code"="GENUS_CODE","StartingYear"="OBS_YEAR"))
  names(ModelParams_SIG)[which(names(ModelParams_SIG)=="Nrec")]="N.Rec_Site"
  names(ModelParams_SIG)[which(names(ModelParams_SIG)=="RecSFM_p_SiteAdcm2_Yr")]="MN.RecVal_Site"
  names(ModelParams_SIG)[which(names(ModelParams_SIG)=="N.RecVal")]="N.RecVal_Sec"
  names(ModelParams_SIG)[which(names(ModelParams_SIG)=="MN.RecVal")]="MN.RecVal_Sec"
  names(ModelParams_SIG)[which(names(ModelParams_SIG)=="LOW_CI95.RecVal")]="LOW_CI95.RecVal_Sec"
  names(ModelParams_SIG)[which(names(ModelParams_SIG)=="HIGH_CI95.RecVal")]="HIGH_CI95.RecVal_Sec"
  #add SSSS (Scleractinia) sector rec values
  Rec_SSSS=subset(Jsec_P_Asec,GENUS_CODE=="SSSS")
  names(Rec_SSSS)[21:26]=paste0(names(Rec_SSSS)[21:26],"_Sec_SSSS")
  ModelParams_SIG= left_join(ModelParams_SIG,
                             Rec_SSSS[,c("ANALYSIS_SEC","OBS_YEAR",
                                         "N.RecVal_Sec_SSSS","MN.RecVal_Sec_SSSS","LOW_CI95.RecVal_Sec_SSSS","HIGH_CI95.RecVal_Sec_SSSS")],
                             by=c("SEC_NAME"="ANALYSIS_SEC","StartingYear"="OBS_YEAR"))
  ModelParams_SIG$ISLAND=substr(ModelParams_SIG$Site,1,3)
  ModelParams_SIG$REGION=RegLU[ModelParams_SIG$ISLAND]
  ModelParams_SIG= left_join(ModelParams_SIG,
                             RecVal_Sec_Dists[,c("REGION","GENUS_CODE",
                                                 "MD.RecVal_Sec_All","MD.CI95_LO.RecVal_Sec_All","MD.CI95_HI.RecVal_Sec_All")],
                             by=c("REGION"="REGION","Genus_Code"="GENUS_CODE"))
  
  table(!is.na(ModelParams_SIG$MN.RecVal_Site))
  table(!is.na(ModelParams_SIG$MN.RecVal_Sec))
  table(!is.na(ModelParams_SIG$MN.RecVal_Site)&!is.na(ModelParams_SIG$MN.RecVal_Sec))
  table(!is.na(ModelParams_SIG$MD.RecVal_Sec_All))
  
  #Save
  save(ModelParams_SIG, file = sprintf("2-VitalRateFunctions/Output/%s_allmodfits_noNAs.rdata",Name))
  #optional: save final vital rate model parameters as csv
  #write.csv(ModelParams_SIG,file = "2-VitalRateFunctions/Output/ModelParams_SIG.csv")
  
  ModelParams_SIG %>% 
    group_by(Site,StartingYear) %>% 
    summarize(N=sum(g.N)) %>% 
    pivot_wider(names_from = StartingYear, values_from = N)
 
     
  
  ################################################################################
  ######### Regional + Genus models #################################
  ########## INTEGRATION ! #######################################################
  ################################################################################
  
  #Regional models (genus and region)
  combo_regional <- left_join(GrowthRegionalDF,RegionalAnnSurvDataFrame)
  combo_regional <- na.omit(combo_regional)
  
  #add SfM site-level stock recruitment data
  ModelParams_Regional = left_join(combo_regional,RecSFMSummary, by= c("Genus_Code","REGION") )
  
  #add REA sector-level stock recruitment data
  ModelParams_Regional = left_join(ModelParams_Regional,Regional_SectorStockRec,by=c("Genus_Code"="GENUS_CODE","REGION"))
  
  save(ModelParams_Regional, file = sprintf("2-VitalRateFunctions/Output/%s_allModFits_noNAs_regional.rdata",Name))
  

