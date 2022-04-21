
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
trier<-askYesNo("Do you want to sort by genotype belle fleur ?")
while(is.na(trier)==TRUE){
  trier<-askYesNo("Do you want to sort by genotype ? 
                  I would like an answer")
}
if(trier==TRUE){
  rentrertrier<-askYesNo("do you want to enter the genotypes ?")
  while(is.na(rentrertrier)==TRUE){
    trier<-askYesNo("do you want to enter the genotypes bonne mère ? 
                  I would like an answer")
  }
}

n=FALSE  # initialisation des valeurs wt
m=FALSE  # initialisation des valeurs mutants
if(rentrertrier==TRUE){  # boucle pour selectionner a la main les wt et mutants
  while(n!=TRUE|m!=TRUE|is.na(n)==TRUE|is.na(m)==TRUE){
    r=1
    while(r>0){
      wt<-tk_select.list(c(1:96),preselect=NULL,multiple=TRUE,title="select les wt je vous prie")
      mut<-tk_select.list(c(1:96),preselect=NULL,multiple=TRUE,title="select les mutants je vous prie")
      r=sum(wt%in%mut)  # verifer que wt different de mutant
  }
    n<-askYesNo(paste("is there",length(wt),"wt, tabernacle ?"))  # verifier le nombre de wt
    m<-askYesNo(paste("is there",length(mut),"mutants, saperlipopette ?"))  # verifier nombre de mutants
  }
}


# 1. DEFINE FOLDERS ----
setwd("Y://Clair//for_Elise//align_files//not_aligned")
dossier_sortie<-"Y://Clair//for_Elise//align_files//aligned"


# 2. LOADER LES FICHIERS ----
nom_fichiers<-list.files(pattern = "^[^~]")  # sans prendre fichiers temporaires
nom_fichiers
fichiers<-list()
for (i in 1:length(nom_fichiers))
{
  fichiers[[i]]<-read_excel(nom_fichiers[i],col_names = FALSE)
}
names(fichiers)<-sub('\\.xlsx$', '', nom_fichiers)


# 3. CLEAN FILES ----
# detect time, arena, genotype and distance columns
coltemps<-as.numeric(which(sapply(fichiers[[1]], function(x) any(x == "Start-0:01:00"))))
colarene<-coltemps-1
coldistance<-as.numeric(which(sapply(fichiers[[1]], function(x) any(x == "Distance moved"))))
colgeno<-as.numeric(which(sapply(fichiers[[1]], function(x) any(x == "Independent Variable"))))
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
narene<-nrow(unique(fichiers[[1]][,colarene]))  # nombre d'arenes
tableaux<-list()
for(i in 1:length(nom_fichiers)){
  temps<-unique(fichiers[[i]][,coltemps])  # extraire temps
  tableaux[[i]]<-cbind(temps)  # creer data frame avec colonne temps
  for(j in 1:narene){
    data<-as.data.frame(filter(fichiers[[i]],fichiers[[i]][,colarene]==j)%>%select(coltemps,coldistance))
    data[,2]<-as.numeric(data[,2])
    tableaux[[i]]<-merge(tableaux[[i]],data,by.x=1,by.y=1,all=T,sort=TRUE)
  }
  colnames(tableaux[[i]])<-c("temps",unlist(unique(fichiers[[1]][,colarene])))  # mettre bon nom de colonne
  # rajouter un 0 au temps si en dessous de 10:00:00
  for(j in 1:nrow(tableaux[[i]])){
    if(nchar(tableaux[[i]][j,1])<17){
      tableaux[[i]][j,1]<-paste("0",tableaux[[i]][j,1],sep="")
    }
  }
  # ordonner comme il faut les tableaux (pour que le temps 10:00:00 vienne après 09:00:00 et pas après 01:00:00)
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
