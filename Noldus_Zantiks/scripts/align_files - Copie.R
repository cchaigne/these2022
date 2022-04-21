
# BE SURE TO EMPTY THE 'NOT_ALIGNED' FOLDER BEFORE RUNNING THIS SCRIPT#

rm(list=ls())  # clean R environment

# 0.1 LIBRAIRIES ----
library(readxl)
library(writexl)
library(dplyr)
library(tcltk)
library(svDialogs)

# 0.2 SORT OR NOT BY GENOTYPE ----
# ask if you want to sort by genotype
#if yes asks the genotypes and arenas numbers
trier <- askYesNo("Do you want to sort by genotype belle fleur ?")
while(is.na(trier)==TRUE){
  trier <- askYesNo("Do you want to sort by genotype ? 
                  I would like an answer")
}
if(trier==TRUE){
  rentrertrier<-askYesNo("do you want to enter the genotypes ?")
  while(is.na(rentrertrier)==TRUE){
    trier <- askYesNo("do you want to enter the genotypes bonne mère ? 
                  I would like an answer")
  }
}

n=FALSE  # initialisation of wt
m=FALSE  # initialisation of mutants
if(rentrertrier==TRUE){  # select wt and mutants by hand
  while(n!=TRUE|m!=TRUE|is.na(n)==TRUE|is.na(m)==TRUE){
    r=1
    while(r>0){
      wt <- tk_select.list(c(1:96),preselect=NULL,multiple=TRUE,title="select les wt je vous prie")
      mut <- tk_select.list(c(1:96),preselect=NULL,multiple=TRUE,title="select les mutants je vous prie")
      r=sum(wt%in%mut)  # verifer que wt different de mutant
  }
    n <- askYesNo(paste("is there",length(wt),"wt, tabernacle ?"))  # check number of wt
    m <- askYesNo(paste("is there",length(mut),"mutants, saperlipopette ?"))  # check number of mutants
  }
}


# 1. DEFINE FOLDERS ----
setwd("Y://Clair//for_Elise//align_files//not_aligned")
dossier_sortie<-"Y://Clair//for_Elise//align_files//aligned"


# 2. LOADER FILES ----
nom_fichiers <-list.files(pattern = "^[^~]")  # load files without temporary files
nom_fichiers
fichiers <-list()
for (i in 1:length(nom_fichiers))
{
  fichiers[[i]] <- read_excel(nom_fichiers[i],col_names = FALSE)
}
names(fichiers) <- sub('\\.xlsx$', '', nom_fichiers)


# 3. CLEAN FILES ----
# detect time, arena, genotype and distance columns
coltemps <- as.numeric(which(sapply(fichiers[[1]], function(x) any(x == "Start-0:01:00"))))
colarene <- coltemps-1
coldistance <- as.numeric(which(sapply(fichiers[[1]], function(x) any(x == "Distance moved"))))
colgeno <- as.numeric(which(sapply(fichiers[[1]], function(x) any(x == "Independent Variable"))))
# time : change "Start-0:01:00" en "0:00:00-0:01:00"
for(i in 1:length(nom_fichiers)){
  fichiers[[i]][,coltemps][fichiers[[i]][,coltemps]=="Start-0:01:00"]<-"0:00:00-0:01:00"
}
# distance : change "-" en "0"  
for(i in 1:length(nom_fichiers)){
  fichiers[[i]][,coldistance][fichiers[[i]][,coldistance]=="-"]<-"0"
}
# remove useless line
ligne<-which(fichiers[[1]][coldistance]=="mm")
for(i in 1:length(nom_fichiers)){
  fichiers[[i]]<-fichiers[[i]][-c(1:ligne),]
}
rm(ligne)


# 4. TABLEAU CROISE DYNAMIQUE ----

# detect arenas in each file and keep those who are present in all the files
arenas1 <- unlist(unique(fichiers[[1]][, colarene]))
for (i in 1:(length(fichiers)-1)){
  arenas <- vector()
  arenas2 <- unlist(unique(fichiers[[i+1]][, colarene]))
  for (j in 1:length(arenas1)){
    if (arenas1[j] %in% arenas2){
      arenas <- c(arenas, arenas1[j])
    }
  }
  arenas1 <- arenas
}

# tableau croise dynamique
tableaux <-list()
for(i in 1:length(nom_fichiers)){
  temps <- unique(fichiers[[i]][,coltemps])  # extract time
  tableaux[[i]] <- cbind(temps)  # crete data frame with time column
  for(j in 1:length(arenas1)){
    data <- as.data.frame(filter(fichiers[[i]],fichiers[[i]][,colarene]==arenas1[j])%>%select(coltemps,coldistance))
    data[,2] <- as.numeric(data[,2])
    tableaux[[i]] <- merge(tableaux[[i]],data,by.x=1,by.y=1,all=T,sort=TRUE)
  }
  colnames(tableaux[[i]]) <-c("temps",arenas1)  # add correct column names
  # add 0 for time if below 10:00:00
  for(j in 1:nrow(tableaux[[i]])){
    if(nchar(tableaux[[i]][j,1])<17){
      tableaux[[i]][j,1]<-paste("0",tableaux[[i]][j,1],sep="")
    }
  }
  # order time (10:00:00 after 09:00:00 and not after 01:00:00)
  tableaux[[i]]<-tableaux[[i]][order(tableaux[[i]]$temps),]
}


# 5. ALIGN FILES ----
alignement<-data.frame()
for(i in 1:length(nom_fichiers)){
  alignement<-rbind(alignement,tableaux[[i]])
}


# 5.BIS SORT BY GENOTYPE ----
if(rentrertrier==TRUE){  # id genotypes if they were not entered by hand
  colgeno<-coldistance-1
  arenes<-which(fichiers[[1]][,coltemps]=="0:00:00-0:01:00")
  geno<-fichiers[[1]][arenes,colgeno]
  d<-cbind(fichiers[[1]][arenes,colarene],geno)
  colnames(d)<-c("arenes","genotype")
  wt<-as.numeric(unlist(filter(d,genotype=="wt")%>%select(arenes)))
  mut<-as.numeric(unlist(filter(d,genotype=="mut")%>%select(arenes)))
}
if(trier==TRUE){
  alignementtrie<-data.frame(temps=alignement[,1],wt=alignement[,(wt+1)],mut=alignement[,(mut+1)])
  alignementtrie<-cbind(alignementtrie,moyennewt=rowMeans(alignementtrie[,c(2:(length(wt)+1))]),moyennemut=rowMeans(alignementtrie[,c((length(wt)+2):ncol(alignementtrie))]))
}


# 6. EXPORTER ALIGN FILE ----

# demander a nommer le fichier d'alignement
nomfichier <- dlgInput("How do you want to name the aligned file ma douce ?", Sys.info()["nomfichier"])$res

if(trier==TRUE){
  write_xlsx(alignementtrie, path = paste(dossier_sortie,"//",nomfichier,".xlsx",sep=""), col_names = TRUE)
}else{
  write_xlsx(alignement, path = paste(dossier_sortie,"//",nomfichier,".xlsx",sep=""), col_names = TRUE)
  }
