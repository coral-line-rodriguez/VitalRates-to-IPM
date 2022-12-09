rm(list=ls())
# Libraries ---------------------------------------------------------------
library(igraph)
library(lubridate)
library(ggplot2)
library(ggpubr)
library(plyr)
library(stringr)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(flextable)
library(officer)


# ## Functions ------------------------------------------------------------

#Pass dataset of single site, multiple timepoint data.frame "Patches"
#Must Have first column as "PatchID"
Patches2ColonyGraphs=function(Patches){
  require(igraph)
  require(lubridate)
  # Error Handling -------------------------------------------------
  #Force PatchName to first column
  #patch first
  if(which(names(Patches)=="PatchName")!=1){stop("Patches2ColonyGraphs: ERROR: 'Patches' data frame must have first column with unique Patch Names  named 'PatchName'...")}
  #need cols
  neededcols=c("LinkTP_For","LinkTP_Bac","Date","TimePt","Circrat","PatchID","PatchName","Notes","Shape_Area")#,"FID"
  missingcols=setdiff(neededcols,names(Patches))
  if(length(missingcols)>0){
    missingcols=setdiff(neededcols,names(Patches))
    stop(paste0("Patches2ColonyGraphs: ERROR: 'Patches' data frame must have columns: ",paste0(missingcols,collapse=",")))
  }
  # How Many Digits in PatchID -------------------------------------------------
  NIDdigit=mean(nchar(Patches$PatchID)-2) #assuming "P_" then leading zero padded number for patchID
  if(any(nchar(Patches$PatchID)!=(NIDdigit+2))){
    stop(paste0("Patches2ColonyGraphs: ERROR: This code assumes that PatchIDs start with 'P_' and then have a leading-zero-padded number. You've got variable length IDs, which is going to lead to issues."))
  }
  #CleanDates
  Patches$R_DATE=mdy(Patches$Date)
  #Add Type
  Patches$Type="LIVE"
  
  # Build To/From DataFrame -------------------------------------------------
  Transitions=NULL
  print(paste("Building Transitions from",nrow(Patches),"patches."))
  for (i in 1:nrow(Patches)){
    #PatchID-->PatchName
    P_Home=Patches$PatchName[i]
    PN_sp=unlist(strsplit(x=P_Home,split="_P_"))
    PN_pre=paste0(PN_sp[1],"_")
    PN_post=paste0("_",unlist(strsplit(PN_sp[2],"_"))[2])
    P_To=unlist(strsplit(x = as.character(Patches$LinkTP_For[i]),split = ","))
    P_From=unlist(strsplit(x = as.character(Patches$LinkTP_Bac[i]),split = ","))
    if(!all(is.na(P_To))){#If there's anything aside from NA in "P_To"
      P_To=na.omit(P_To)
      P_To[which(P_To!="-99")]=paste0(PN_pre,"P_",leadz(P_To,n=NIDdigit),PN_post)
      P_To[which(P_To=="-99")]=paste0("MORT_",P_Home)
      T_To=data.frame(From=P_Home,To=P_To,stringsAsFactors = F)
      Transitions=rbind(Transitions,T_To)
    }
    if(!all(is.na(P_From))){#If there's anything aside from NA in "P_From"
      P_From=na.omit(P_From)
      P_From[which(P_From!="-99")]=paste0(PN_pre,"P_",leadz(P_From,n=NIDdigit),PN_post)
      P_From[which(P_From=="-99")]=paste0("RECR_",P_Home)
      T_From=data.frame(From=P_From,To=P_Home,stringsAsFactors = F)
      Transitions=na.omit(rbind(Transitions,T_From))
    }
  }
  #Shoudl expect many repeated to-from links
  Transitions=unique(Transitions)
  print(paste("Generated",nrow(Transitions),"unique transitions from",nrow(Patches),"patches."))
  
  #   #Add Patches for -99s (recruitment/mortality) -------------------------
  #Figure out patches to add (names)
  AllPatches=sort(unique(c(Transitions[,1],Transitions[,2])))
  N99raw=setdiff(AllPatches,Patches$PatchName)
  
  #Here NPP should flag any transition that is either RECR MORT or some breakdown between PatchID and PatchName
  #This last group can include (1) Trans species links, and so other issues. For now, I'm deleting the non MORT/RECR transitions, reporting and moving on 
  Problem_i=setdiff(1:length(N99raw),c(grep(pattern = "MORT",x=N99raw),grep(pattern = "RECR",x=N99raw)))
  
  if(length(Problem_i)>0){
    print(paste("Dropping",length(Problem_i),"Transitions with problematic PatchNames from",nrow(Transitions),"total."))
    print(paste("Examples include:",paste(N99raw[sample(Problem_i,20,replace = T)],collapse=", ")))
    
    #Set N99
    N99=N99raw[-Problem_i]
    #Remove from Transitions
    Problem_Trow=sort(unique(c(which(Transitions[,1]%in%N99raw[Problem_i]),which(Transitions[,2]%in%N99raw[Problem_i]))))
    if(length(Problem_Trow)>0){Transitions=Transitions[-Problem_Trow,]}
  }else{
      N99=N99raw
  }
                
  if(length(N99)>0){
    N99ID=substr(x = N99,start = 6,stop = nchar(N99))
    N99TYPE=substr(x = N99,start = 1,stop = 4)
    
    #Build and modify DF
    N99DF=Patches[match(N99ID,Patches$PatchName),]
    #Get TimePoint LookUp
    TPlu=unique(N99DF[,c("Site","TimePt","Date","R_DATE")])
    TPlu$STP=paste0(TPlu$Site,"_",TPlu$TimePt)
    
    #Modify Patch Data - PatchID
    N99DF$PatchName=N99
    #Modify Patch Data - PatchID
    N99DF$PatchID=substr(N99ID,start = 18,nchar(N99ID)-2)
    #Modify Patch Data - FID
    N99DF$FID=max(Patches$FID,na.rm=T)+1:nrow(N99DF)
    #Modify Patch Data - Site
    N99DF$Site=substr(N99ID,1,11)
    #Modify Patch Data - Circat
    N99DF$Circrat=NA
    #Modify Patch Data - Shape_Area
    N99DF$Shape_Area=0
    #Modify Patch Data - LinkTP_Bac
    N99DF[which(N99TYPE=="MORT"),"LinkTP_Bac"]=N99ID[which(N99TYPE=="MORT")]
    N99DF[which(N99TYPE=="RECR"),"LinkTP_Bac"]=NA
    #Modify Patch Data - LinkTP_For
    N99DF[which(N99TYPE=="RECR"),"LinkTP_For"]=N99ID[which(N99TYPE=="RECR")]
    N99DF[which(N99TYPE=="MORT"),"LinkTP_For"]=NA
    #Modify Patch Data - TimePt
    N99DF[which(N99TYPE=="MORT"),"TimePt"]=N99DF[which(N99TYPE=="MORT"),"TimePt"]+1
    N99DF[which(N99TYPE=="RECR"),"TimePt"]=N99DF[which(N99TYPE=="RECR"),"TimePt"]-1
    #Modify Patch Data - Date
    N99STP=paste0(N99DF$Site,"_",N99DF$TimePt)
    N99DF$Date=TPlu[match(N99STP,TPlu$STP),"Date"]
    N99DF$R_DATE=TPlu[match(N99STP,TPlu$STP),"R_DATE"]
    #Modify Patch Data - Notes
    N99DF$Notes=N99TYPE
    N99DF$Type=N99TYPE
    
    Patches_rm=rbind(Patches,N99DF)
  }else{
    Patches_rm=Patches
  }
  
  g1 <- graph_from_data_frame(Transitions, directed=TRUE, vertices=Patches_rm)
  
  E(g1)$AreaChange_cm2=10^4*(head_of(g1,E(g1))$Shape_Area-tail_of(g1,E(g1))$Shape_Area)
  E(g1)$AreaChangeSign=sign(E(g1)$AreaChange_cm2)
  
  #Set Transition Types
  E(g1)$TRANSITION_TYPE=NA
  E(g1)$TRANSITION_TYPE[E(g1)$AreaChangeSign>=0]="GROWTH"
  E(g1)$TRANSITION_TYPE[E(g1)$AreaChangeSign<0]="SHRINK"
  
  MORT_i=which(head_of(g1,E(g1))$Type=="MORT")
  E(g1)$TRANSITION_TYPE[MORT_i]="MORT"
  RECR_i=which(tail_of(g1,E(g1))$Type=="RECR")
  E(g1)$TRANSITION_TYPE[RECR_i]="RECR"
  
  #FUsion/FIssion are slightly harder because they apply in Edge-index space, but assess in Vertex-index space
  #Mark Every Edge From a vertex with >1 out going connections a "FISSION" edge
  FISSION_Vi=which(degree(g1,v=unique(tail_of(graph = g1,es=E(g1))),
                          mode="out")>1)#index in V() space
  E(g1)[.from(names(FISSION_Vi))]$TRANSITION_TYPE="FISSION"
  #Mark Every Edge From a vertex with >1 in-coming connections a "FUSION" edge
  FUSION_Vi=which(degree(g1,
                         v=unique(head_of(graph = g1,es=E(g1))),
                         mode="in")>1)
  E(g1)[.to(names(FUSION_Vi))]$TRANSITION_TYPE="FUSION"
  
  #Mark Every Edge Listed as both FUSION and FISSION a FUSION_FISSION Edge
  FUSFIS_Ei=intersect(E(g1)[.from(names(FISSION_Vi))], E(g1)[.to(names(FUSION_Vi))])
  E(g1)[FUSFIS_Ei]$TRANSITION_TYPE="FUSION_FISSION"
  
  E(g1)$Interval_Years=(head_of(g1,E(g1))$R_DATE-tail_of(g1,E(g1))$R_DATE)/365.25
  E(g1)$DataOrError=as.vector(cut(E(g1)$Interval_Years,breaks=c(0,7/365.25,90/365.25,999),labels=c("ERROR","NEITHER","DATA"),include.lowest=T))
  E(g1)$CrossAnnotator=head_of(g1,E(g1))$Annotation_Round!=tail_of(g1,E(g1))$Annotation_Round
  E(g1)$GrowthRate_cm2_yr=E(g1)$AreaChange_cm2/E(g1)$Interval_Years
  
  #Assign ColonyID to Vertex
  V(g1)$ColonyID="NA"
  sg1=decompose.graph(g1,min.vertices = 2)
  for(icol in 1:length(sg1)){
    thisgraph=sg1[[icol]]
    V(g1)[which(V(g1)$name%in%V(thisgraph)$name)]$ColonyID=paste0("C_",formatC(icol,width=5,flag=0),"_",V(thisgraph)$Site[1])
  }
  return(g1)
}

#Plotting Support
NodeTypeColor=c("green","red","pink")
names(NodeTypeColor)=c("LIVE","MORT","RECR")

EdgeTypeColor=c("darkblue","lightblue","brown","yellow","red","pink","orange")
names(EdgeTypeColor)=c("GROWTH","SHRINK","FUSION","FISSION","MORT","RECR","FUSION_FISSION")

SizeScaleFunc=function(x,lowest=5,highest=15,trans="none"){
  if(trans=="log"){
    x_l=log(x+1)
    x_sc=(x_l-min(x_l,na.rm=T))/(max(x_l,na.rm=T)-min(x_l,na.rm=T))
  }else{
    x_sc=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}
  sz=lowest+x_sc*(highest-lowest)
  return(sz)
}

subg_Vname=function(graph,vertex_names){
  sg1 <- decompose.graph(graph,mode="weak")
  subfollow=function(s){if(any(V(s)$name %in% vertex_names)) V(s)$name else NULL}  
  neighverts <- unique(unlist(sapply(sg1,FUN=subfollow)))
  gout <- induced.subgraph(graph=graph,vids=neighverts)
  return(gout)
}

subg_Vattr=function(graph,attr,comparison,INFO=T){
  sg1 <- decompose.graph(graph,mode="weak")
  eval(parse(text=paste0("vertex_names=V(graph)$name[which(V(graph)$",attr,comparison,")]")))
  subfollow=function(s){if(any(V(s)$name %in% vertex_names)) V(s)$name else NULL}  
  neighverts <- unique(unlist(sapply(sg1,FUN=subfollow)))
  gout <- induced.subgraph(graph=graph,vids=neighverts)
  if(INFO){print(paste0(length(decompose.graph(gout,mode="weak"))," subgraphs of ",length(sg1)," with feature identified."))}
  return(gout)
}

subg_Eind=function(graph,edge_indices){
  verts=V(graph)[.inc(E(graph)[edge_indices])]$name
  return(subg_Vname(graph,verts))
}

as.nv=function(x){return(as.numeric(as.vector(x)))}
leadz=function(x,n){return(formatC(as.numeric(as.vector(x)),width=n,flag=0))} #formats numbers individually
#Converts value to vector and then converts to numeric vector
#if the value is a charcter/word, returns NAs and gives warning message NAs introduced by coercion

mod1=function(x,k){return(1+mod(x-1,k))}
plot.colgraph=function(graph,NodeColor=NodeTypeColor,EdgeColor=EdgeTypeColor,SizeFunc=SizeScaleFunc,layoutfun=layout_nicely){
  shapes=c("circle", "square","pie","sphere", "rectangle","crectangle", "vrectangle")
  #eval(parse(text=paste0("glay=",layoutfun,"(graph)")))
  glay=layoutfun(graph)
  plot.igraph(graph,#layout=glay,
              edge.width=SizeFunc(E(graph)$AreaChange_cm2,lowest=.05,highest=6),
              edge.arrow.width=SizeFunc(abs(E(graph)$AreaChange_cm2),lowest=.05,highest=4),
              edge.arrow.size=SizeFunc(abs(E(graph)$AreaChange_cm2),lowest=.05,highest=1),
              edge.color=EdgeColor[E(graph)$TRANSITION_TYPE],
              vertex.color=NodeColor[V(graph)$Type],
              vertex.pie.color=NodeColor[V(graph)$Type],
              vertex.size=SizeFunc(V(graph)$Shape_Area,lowest=3,highest=10,trans="log"),
              vertex.shape=shapes[mod1(V(graph)$TimePt,length(shapes))],
              vertex.label.cex=0.75,
              vertex.label.color="black",
              vertex.label.dist=1.5)
}

# Load Datasets -----------------------------------------------------------
basepath="/Users/c-rod/Documents/GitHub/VitalRates-to-IPM/1-Annotations_ColonyTransitions/Data/"  


#Set up to loop through files
fl=list.files(path = basepath,pattern = ".csv",full.names = T)
#Need to distinguish annotator by file
annotation_num=rep(1,length(fl))
annotation_num[grep(pattern = "VitalRatesQC",x =fl)]=2

CanonColNames = c("FID","Site","Date","TimePt","Annotation_Round","Circrat","Empty_Skip_Circrat_Flag","PatchID","PatchName",
                  "Spec_Code","Morph_Code","BLE_EXT","BLE_SEV","Shape_Length","Shape_Area","Shape_Diameter","LinkTP_Bac","LinkTP_For","Notes","FileName")
# extra link columns needed becasue they ran out of room, need to merege

#ColonyID is the Empty flag..
#Site _ PatchID should be unique ##### Dups are likely due to double annotatotion

#load 'em all in storing different columnnames/orders in differentlists
for(i_f in 1:length(fl)){
  thissite=read.csv(fl[i_f],stringsAsFactors = FALSE)
  thissite$Annotation_Round=annotation_num[i_f]
  thissite$FileName=strsplit(fl[i_f],"/")[[1]][6]
  
  names(thissite)
#Catch easy naming errors
  if(any(names(thissite)%in%"ColonyID")){names(thissite)[which(names(thissite)=="ColonyID")]="Empty_Skip_Circrat_Flag"}
  if(any(names(thissite)%in%"OBJECTID")){names(thissite)[which(names(thissite)=="OBJECTID")]="FID"}
  if(any(names(thissite)%in%"SpeciesCod")){names(thissite)[which(names(thissite)=="SpeciesCod")]="Spec_Code"}
  if(any(names(thissite)%in%"Poly_Area")){names(thissite)[which(names(thissite)=="Poly_Area")]="Shape_Area"}
  if(any(names(thissite)%in%"Poly_Perim")){names(thissite)[which(names(thissite)=="Poly_Perim")]="Shape_Length"}
  if(any(names(thissite)%in%"Shape_Leng")){names(thissite)[which(names(thissite)=="Shape_Leng")]="Shape_Length"}
  if(any(names(thissite)%in%"MBG_Diamet")){names(thissite)[which(names(thissite)=="MBG_Diamet")]="Shape_Diameter"}
  if(any(names(thissite)%in%"Max_Diam")){names(thissite)[which(names(thissite)=="Max_Diam")]="Shape_Diameter"}
  if(any(names(thissite)%in%"Max_Diameter")){names(thissite)[which(names(thissite)=="Max_Diameter")]="Shape_Diameter"}
  names(thissite)
  
  #Take care of NA in EMPTY
  thissite$Empty_Skip_Circrat_Flag[is.na(thissite$Empty_Skip_Circrat_Flag)]="NOT_EMPTY"
  thissite$Empty_Skip_Circrat_Flag[str_trim(thissite$Empty_Skip_Circrat_Flag)==""]="NOT_EMPTY"
  
  #Drop any BS extra columns
  dropcols=c("Id","Shape_Area.1","Shape_Leng.1","MBG_Diamet","ORIG_FID","Area","MBG_Diameter")
  thissite=thissite[,setdiff(names(thissite),dropcols)]
  names(thissite)
  
  #If there's no bleaching or MorphCode or Diameter data, add NA cols
  if(!"BLE_EXT"%in%names(thissite)){thissite$BLE_EXT=NA}
  if(!"BLE_SEV"%in%names(thissite)){thissite$BLE_SEV=NA}
  if(!"Morph_Code"%in%names(thissite)){thissite$Morph_Code=NA}
  if(!"Shape_Diameter"%in%names(thissite)){thissite$Shape_Diameter=NA}
  names(thissite)
  
  ###Tidy Link Columns:
  #Get blank and NA cells to null string ""
  LinkCols_i=grep(pattern = "LinkTP_",names(thissite))
  for(lc_i in 1:length(LinkCols_i)){thissite[which(thissite[,LinkCols_i[lc_i]]==" "|is.na(thissite[,LinkCols_i[lc_i]])),LinkCols_i[lc_i]]=""}
  
  #Concatenate multiple link columns
  if(length(grep(pattern = "LinkTP_F",names(thissite)))>1){
    print(paste(i_f,":",fl[i_f],names(thissite)[grep(pattern = "LinkTP_F",names(thissite))]))
    LF=paste(thissite$LinkTP_For,thissite$LinkTP_F_1,sep=",")
    LF[which(LF==",")]=""
    thissite$LinkTP_For=LF
    thissite=thissite[,setdiff(names(thissite),"LinkTP_F_1")]
  }
  if(length(grep(pattern = "LinkTP_B",names(thissite)))>1){
    print(paste(i_f,":",fl[i_f],names(thissite)[grep(pattern = "LinkTP_B",names(thissite))]))
    LB=paste(thissite$LinkTP_Bac,thissite$LinkTP_B_1,sep=",")
    LB[which(LB==",")]=""
    thissite$LinkTP_Bac=LB
    thissite=thissite[,setdiff(names(thissite),"LinkTP_B_1")]
  }
  #Make sure all periods are commas, drop timepoint after -99s
  perstr_iF=grep(thissite$LinkTP_For,pattern="\\.")
  perstr_iB=grep(thissite$LinkTP_Bac,pattern="\\.")
  #for each period in a string, either turn any -99.X into just -99s, or change any periods into commas
  #For For
  if(length(perstr_iF)>0){
    for(i_ps in 1:length(perstr_iF)){
      cellstr=thissite$LinkTP_For[perstr_iF[i_ps]]
      cellstr=gsub(cellstr,pattern="\\..",replacement="\\.")
      thisps=strsplit(cellstr,",")
      tpsnum=lapply(thisps,as.numeric)
      #from list find/replace in strings the -99s
      l1_i=which(lapply(thisps,length)==1)
      if(length(l1_i)>0) {neg99s_i=l1_i[tpsnum[[l1_i]]<0]}else{neg99s_i=NULL}
      if(length(neg99s_i)>0){thisps[[neg99s_i]]="-99"} else {print(paste0("Problem with this csv file: ", fl[i_f], " at line ", perstr_iF[i_ps], " for these problematic linkages: ", cellstr))}
      #moosh it back together, then replace any . with ,
      thisstr=paste(unlist(thisps),collapse=",")
      thissite$LinkTP_For[perstr_iF[i_ps]]=gsub(thisstr,pattern="\\.",replacement="\\,")
    }
  }
  
  #For Bac
  if(length(perstr_iB)>0){
    for(i_ps in 1:length(perstr_iB)){
      cellstr=thissite$LinkTP_Bac[perstr_iB[i_ps]]
      cellstr=gsub(cellstr,pattern="\\..",replacement="\\.")
      thisps=strsplit(cellstr,",")
      tpsnum=lapply(thisps,as.numeric)
      #from list find/replace in strings the -99s
      l1_i=which(lapply(thisps,length)==1)
      if(length(l1_i)>0) {neg99s_i=l1_i[tpsnum[[l1_i]]<0]}else{neg99s_i=NULL}
      if(length(neg99s_i)>0){thisps[[neg99s_i]]="-99"} else {print(paste0("Problem with this csv file: ", fl[i_f], " at line ", perstr_iB[i_ps], " for these problematic linkages: ", cellstr))}
      #moosh it back together, then replace any . with ,
      thisstr=paste(unlist(thisps),collapse=",")
      thissite$LinkTP_Bac[perstr_iB[i_ps]]=gsub(thisstr,pattern="\\.",replacement="\\,")
    }
  }
  
  #Prep for final addition
  thissite=thissite[,sort(names(thissite))]
  
  #Here's where we actually add the new data into our datastructure
  if(i_f==1){ # if this is the first time we're doing this, just set up a list for the data frame and one for the column names
    formatlist=list(names(thissite))
    stacklist=list(thissite)
  }else{ #if this isnt the first time, it gets a bit more complicated...
    formatmatch=NULL
    for(i_for in 1:length(formatlist)){#for each element in the format list
      #check if the current dataframe's columnames match the list-stored sets of columnames
      formatmatch=c(formatmatch,all(formatlist[[i_for]]==names(thissite)))
    }
    formi=which(formatmatch)
    if(length(formi)==0){#There is no match format existing, New Foramt!
      print(paste0("New Format! i_f =",i_f,". Total of ",length(formatlist)+1))
      formatlist[[length(formatlist)+1]]=names(thissite)
      stacklist[[length(stacklist)+1]]=thissite
    }else{
      stacklist[[formi]]=rbind(stacklist[[formi]],thissite)
    }
  }

  print(paste0(i_f," of ",length(fl)))
  #write out THISSITE WITH dATE CODE
  
}
length(stacklist) #want stacklist to be 1 b/c this means all the formats match
lapply(stacklist,nrow)

formatlist
Patches=stacklist[[1]][,CanonColNames]

#exclude the pseudo quadrats that are used to calculate area surveyed
dim(Patches)
Patches <- Patches %>% filter(Empty_Skip_Circrat_Flag == "NOT_EMPTY")
#Remove patches w/ no site, date, patch link, etc
Patches <- Patches %>% filter(PatchID>0)
dim(Patches)

Patches$R_DATE=mdy(Patches$Date)
#write out PaTCHES WITH dATE CODE

#Clean Up And Re-Org PatchID
Pnames=names(Patches) 
Pnames=Pnames[-which(Pnames=="PatchID")]
Patches$PatchIDraw=Patches$PatchID 
Patches$PatchID=paste0("P_",leadz(Patches$PatchIDraw,n=5))#,"_",Patches$Site,"_",Patches$Annotation_Round)
Patches=Patches[,c("PatchID",Pnames)] 

#sites need to be in alphabetical order!
SiteCoherenceLU=c("MAI_OCC_002","OAH_XXX_022") 
names(SiteCoherenceLU)=unique(Patches$Site) 
Patches$Site=SiteCoherenceLU[Patches$Site] # if SiteCoherence isn't in alphabetical order all of the Sites get shuffled to the wrong PatchName b/c they share the same PatchID

table(Patches$Site)

Patches$PatchName=paste(Patches$Site,Patches$Spec_Code,Patches$PatchID,Patches$Annotation_Round,sep="_")
#Clean Up And Re-Org PatchName
Pnames=names(Patches)
Pnames=Pnames[-which(Pnames=="PatchName")]
Patches=Patches[,c("PatchName",Pnames)]
PatchesWithEmpty=Patches
Patches=subset(PatchesWithEmpty,Empty_Skip_Circrat_Flag=="NOT_EMPTY")

GenusCodeLUnames=c("PLIG","PMEA","PGRA","POCS","POSP","PLUT","PLIC","PLOB","MFLA","MPAT","MOSP","MCAP")
GenusCodeLU=     c("POCS","POCS","POCS","POCS","POSP","POSP","POSP","POSP","MOSP","MOSP","MOSP","MOSP")
names(GenusCodeLU)=GenusCodeLUnames

Patches$Spec_Code=str_trim(Patches$Spec_Code)
drop_SpBlank=which(Patches$Spec_Code=="")
if(length(drop_SpBlank)>0){Patches=Patches[-drop_SpBlank,]}
Patches$Genus_Code=GenusCodeLU[Patches$Spec_Code]
table(Patches$Genus_Code,useNA="always")

BadDates_i=c(which(year(Patches$R_DATE)<2000),which(year(Patches$R_DATE)>2020))
table(Patches$Site[BadDates_i])


# Patches2ColonyGraphs Call -----------------------------------------------

#For every Ann_Round==2 Patch, match it to the linked patches from Ann_Round==1
A2_i=which(Patches$Annotation_Round==2)
#Check Annotator's Patch Notes
AndNote=sort(str_trim(unique(Patches$Notes[A2_i])))
GoodNotes=AndNote[c(1,2,5,6,9,22,23,25,26,27)]
BadNotes=setdiff(AndNote,GoodNotes)
BadNote_i=which(Patches$Notes%in%BadNotes)
BadNote_PatchName=Patches$PatchName[BadNote_i]

#MatchUp Patch Up
PN2=Patches$PatchName[A2_i]
PN=substr(PN2,1,nchar(Patches$PatchName[A2_i])-2)
PN1=paste0(PN,"_1")

#Match 'em up
A1_im=match(PN1,Patches$PatchName)



#Patches to Colony

CGall=Patches2ColonyGraphs(Patches = Patches)

VitalRate_Growth=
    data.frame(Site=head_of(CGall,E(CGall))$Site,
               DataOrError=E(CGall)$DataOrError,
               Annotator_Tail=tail_of(CGall,E(CGall))$Annotation_Round,
               Annotator_Head=head_of(CGall,E(CGall))$Annotation_Round,
               ColonyID=head_of(CGall,E(CGall))$ColonyID,
               T0_PatchName=tail_of(CGall,E(CGall))$name,
               T1_PatchName=head_of(CGall,E(CGall))$name,
               EdgeLabels=paste0(tail_of(CGall,E(CGall))$name,"-",head_of(CGall,E(CGall))$name),
               Spec_Code=tail_of(CGall,E(CGall))$Spec_Code,
               Genus_Code=tail_of(CGall,E(CGall))$Genus_Code,
               AreaSurveyed=0.5*length(unique(V(CGall)$Circrat)),
               StartingDate=as.Date(tail_of(CGall,E(CGall))$R_DATE,origin="1970/01/01"),
               EndingDate=as.Date(head_of(CGall,E(CGall))$R_DATE,origin="1970/01/01"),
               Interval_Years=E(CGall)$Interval_Years,
               StartingSize=10^4*tail_of(CGall,E(CGall))$Shape_Area,
               EndingSize=10^4*head_of(CGall,E(CGall))$Shape_Area,
               StartingPerimeter=10^2*tail_of(CGall,E(CGall))$Shape_Length,
               EndingPerimeter=10^2*head_of(CGall,E(CGall))$Shape_Length,
               StartingMaxDiam=10^2*tail_of(CGall,E(CGall))$Shape_Diameter,
               EndingMaxDiam=10^2*head_of(CGall,E(CGall))$Shape_Diameter,
               TransitionMagnitude=E(CGall)$AreaChange_cm2,
               TransitionRate=E(CGall)$GrowthRate_cm2_yr,
               TransitionType=E(CGall)$TRANSITION_TYPE)
VitalRate_Growth$PercentChange=100*VitalRate_Growth$TransitionMagnitude/VitalRate_Growth$StartingSize
VitalRate_Growth$Log2Ratio_Change=log2(VitalRate_Growth$EndingSize/VitalRate_Growth$StartingSize)
head(VitalRate_Growth)
head(Patches)
save(file = "1-Annotations_ColonyTransitions/Output/PatchTransitions.rdata",list=c("Patches","CGall","VitalRate_Growth"))













# Load and plot -----------------------------------------------------------


load(file = "1-Annotations_ColonyTransitions/Output/PatchTransitions.rdata")
GenusCodeLUnames=c("PMEA","PGRA","POSP","PLUT", "PLIC","MFLA","MPAT", "MOSP","POCS","POSP","PLOB","MCAP")
GenusCodeLU=c("POCS","POCS","POSP","POSP","POSP", "MOSP","MOSP", "MOSP","POCS","POSP","POSP","MOSP")
names(GenusCodeLU)=GenusCodeLUnames


#Convert Species Based PatchName to Genus Based Patch Name-T0
VitalRate_Growth$T0_PatchName=as.character(VitalRate_Growth$T0_PatchName)
VitalRate_Growth$T0_PatchName_G=paste0(substr(VitalRate_Growth$T0_PatchName,1,12),
                                       GenusCodeLU[substr(VitalRate_Growth$T0_PatchName,13,16)],
                                       substr(VitalRate_Growth$T0_PatchName,17,nchar(VitalRate_Growth$T0_PatchName)))
RM_i=which(substr(VitalRate_Growth$T0_PatchName,1,4)%in%c("RECR","MORT"))
VitalRate_Growth$T0_PatchName_G[RM_i]=paste0(substr(VitalRate_Growth$T0_PatchName[RM_i],1,17),
                                       GenusCodeLU[substr(VitalRate_Growth$T0_PatchName[RM_i],18,21)],
                                       substr(VitalRate_Growth$T0_PatchName[RM_i],22,nchar(VitalRate_Growth$T0_PatchName[RM_i])))

#Convert Species Based PatchName to Genus Based Patch Name-T1
VitalRate_Growth$T1_PatchName=as.character(VitalRate_Growth$T1_PatchName)
VitalRate_Growth$T1_PatchName_G=paste0(substr(VitalRate_Growth$T1_PatchName,1,12),
                                       GenusCodeLU[substr(VitalRate_Growth$T1_PatchName,13,16)],
                                       substr(VitalRate_Growth$T1_PatchName,17,nchar(VitalRate_Growth$T1_PatchName)))
RM_i=which(substr(VitalRate_Growth$T1_PatchName,1,4)%in%c("RECR","MORT"))
VitalRate_Growth$T1_PatchName_G[RM_i]=paste0(substr(VitalRate_Growth$T1_PatchName[RM_i],1,17),
                                             GenusCodeLU[substr(VitalRate_Growth$T1_PatchName[RM_i],18,21)],
                                             substr(VitalRate_Growth$T1_PatchName[RM_i],22,nchar(VitalRate_Growth$T1_PatchName[RM_i])))

VRG_GS=subset(VitalRate_Growth,TransitionType%in%c("GROWTH","SHRINK"))
VRG_GSerr=subset(VRG_GS,DataOrError=="ERROR")


