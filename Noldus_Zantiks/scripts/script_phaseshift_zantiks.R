rm(list=ls())  # clean R studio environment
#  Ctrl+Enter to run line by line
#  Run all the lines in libraries and functions sections to load them

# fichier
# fichier2 correspond to fichier with the 300 seconds autoreference gaps filled

# 0.1 LIBRARIES ---- 
library(readxl)
library(dplyr)
library(ggplot2)
library(plotly)
library(data.table)
library(writexl)

# 0.2 FUNCTIONS ----
# Do the sum over n rows
n.colsum = function(df, n){
  # df = data frame to do the sum on
  # n = number of rows to sum
  aggregate(x = df,
            by = list(gl(ceiling(nrow(df)/n), n)[1:nrow(df)]),
            FUN = sum)
}
# Do the mean over n rows
n.colmean<-function(df,n){
  # df = data frame to do the mean on
  # n = number of rows to mean
  aggregate(x = df,
            by = list(gl(ceiling(nrow(df)/n), n)[1:nrow(df)]),
            FUN = mean)
}
# Identify day/night and light transitions (return row numbers in a vector) from the column VARIABLE...5
transitions<-function(fichier){
  # fichier = zantiks file containing the VARIABLE...5 column specifying DAY/NIGHT/PULSE etc...
  v<-vector()
  for (i in 2:nrow(fichier)){
    if(fichier$VARIABLE...5[i-1]!=fichier$VARIABLE...5[i]){
      v<-c(v,i)
    }
  }
  return(v)
}
# Add transitions to the plot (red line represent night to day transitions, purple : day to night, grey : light pulse)
lignestrans<-function(trans,n){
  # the first transition should be a day to night one, the light pulse should correspond to the 5th transition
  # trans = row numbers in the zantiks file (data in seconds) corresponding to the transitions (identified by the transitions function)
  # n = 60 for the minute plot
  # n = 600 for the 10minute plot
  # n = 3600 for the hour plot
  for(i in 1:length(trans)){
    if((i==5)|(i==6)){  # light pulse
      abline(v=trans[i]/n,col="grey")}
    else{
      if((i%%2==0)){  # night to day
        abline(v=trans[i]/n,col="purple")
      }
      else{  # day to night
        abline(v=trans[i]/n,col="red")
      }
    }
  }
}
# Detect autoreference moments (gaps of 300 seconds every hour) 
autoreferences<-function(fichier){
  # fichier = zantiks file
  fichier$TIME <- as.numeric(fichier$TIME)  # transform time column from character to numeric
  autoref <- c()
  for (i in 1:(nrow(fichier)-1)){
    if (fichier$TIME[i+1]-fichier$TIME[i]>60) {  # detect gaps of more than 60 seconds in the TIME column
      autoref <- c(autoref, i)
    }
    return(autoref)
  }
}

# 1. DEFINE FOLDERS ----
setwd("Y://ELISE DATA_17Gb a trier//EXCELL+stats//excell//zebrafish locomotor activity//Zantiks_Clair//phaseshift_zantiks//data")  # folder containing zantiks files
dossier_sortie <-"Y://ELISE DATA_17Gb a trier//EXCELL+stats//excell//zebrafish locomotor activity//Zantiks_Clair//phaseshift_zantiks//resultats"  # folder to put transformed files in (meaned by 10minutes for example)

# 2. LOAD FILE ----
# can take several minutes
fichier <- read_excel("yourfile.xlsx")  # change "yourfile.xlsx" by the name of your file (don't forget the quotation marks)

# 3. CLEAN FILE AND FILL AUTOREFENCE GAPS ----
fichier <- fichier[-(nrow(fichier)), ]  # remove the last line containing NAs
# the next lines are used to make another data frame, fichier2 (in seconds) from fichier (also in seconds), where the 300 seconds autoreference gaps are filed
# (filled with the data from the previous 300 seconds)
autoref <- autoreferences(fichier)  # detect the lines preceding an autoreference
# next line to use if your zantiks script begins with an autoreference just after lightoff (like in phaseshift_3 and phaseshift_4)
fichier2 <- fichier[c(1:300), ]  # copy the first 300 seconds of fichier
fichier2 <- rbind(fichier2, fichier[c(1:autoref[1]), ])  # copy data until the first autoreference
for (i in c(1:((length(autoref))-1))){
  fichier2 <- rbind(fichier2, fichier[c((autoref[i]-299):autoref[i]), ])  # copy and add the data of the previous 300 seconds before the ith autoref to fill the gap
  fichier2 <- rbind(fichier2, fichier[c((autoref[i]+1):autoref[(i+1)]), ])  # add the data from the ith autoref to the next autoref
}
fichier2 <- rbind(fichier2, fichier[c((autoref[length(autoref)]-299):autoref[length(autoref)]), ])  # copy and add the data of the previous 300 seconds before the last autoref to fill the gap
fichier2 <- rbind(fichier2, fichier[c((autoref[length(autoref)]+1):nrow(fichier)), ])  # add the data from the last autoref to the last row

# 4. SPECIFICATION OF GENOTYPES ----
# specify the corresponding arena numbers to the conditions/genotypes you want to make plots from
bad <- c(1, 2, 3, 4, 5)  # arenas to discard from plots
good <- c(1:96)[!(c(1:96)%in%bad)]  # arenas to keep 
pig <- c(92,96)  # arenas containing pigmented larvae
nopig <- notpig<-c(1:96)[!(c(1:96)%in%bad)&!(c(1:96)%in%pig)]  # good arenas not containing pigmented larvae
ctrl <- c(3,18,27,34,35,52,55,73,74,84,91,94)
opn4xa <- c(2,8,16,24,25,33,36,50,57,61,67,77,83)
lak <- c(4,46,47,54,80)
double <- c(7,30,42,45,63,69,76,85,88)

# 5. DETECTION OF TRANSITIONS ----
# detection of transitions row in fichier and fichier2
t<-transitions(fichier)  # detection of transitions row in fichier
t2<-transitions(fichier2) # detection of transitions row in fichier2
# extract row numbers corresponding to the transitions +/- 60 seconds
dtol <- c((t[5]-60):(t[5]+60))  # dark to light
ltod <- c((t[6]-60):(t[6]+60))  # light to dark ATTENTION if used phaseshift_3 use this instead :ltod<-c((t[6]):(t[6]+120))
dtol2 <- c((t2[5]-60):(t2[5]+60))
ltod2 <- c((t2[6]-60):(t2[6]+60))  # ATTENTION if used phaseshift_3 use this instead :ltod2<-c((t2[6]):(t2[6]+120))

# 6. CHECKPOINTS ----
# check that the transitions in fichier2 are spaced appropriately
for(i in c(1:(length(t2)-1))){  # should be 14, 10, 14, ~2, ~2, ~6, 14, 10, 14, 10
  print(i)
  print((t2[i+1]-t2[i])/3600)  # to have the spacing between transitions in hours
}
# check that there is indeed an autoref gap between the autoref rows detected and their respective next rows
fichier$TIME<-as.numeric(fichier$TIME)  # transform time column from character to numeric
check<-c()
for (i in 1:(length(autoref))){
  d<-fichier[c((autoref[i]-10):(autoref[i]+10)), ]
  check<-c(check,(d$TIME[12]-d$TIME[11]))
}
View(check)  # should be around 301 seconds each time
# check that the differences in the number of rows between fichier and fichier2 are correct
# for 106 autoref (phaseshift_3 and phaseshift_4) we should have (number of autoref * length of autoref + length of the first autoref) 106*300+300=32100 seconds of differences
nrow(fichier2)-nrow(fichier)  # ok 32100 secondes


# 7. TRANSFORM FILE (sum over n min, mean over n min per min) ----
# sum over 1 minute
fichiermin<-n.colsum(fichier[,c(8:ncol(fichier))],n=60)  # our first arena is in the 8th column
fichiermin2<-n.colsum(fichier2[,c(8:ncol(fichier2))],n=60)
# mean over 10 minutes (per minutes)
fichier10min<-n.colmean(fichiermin[,c(2:ncol(fichiermin))],n=10)  # our first arena is in the 2nd column
fichier10min2<-n.colmean(fichiermin2[,c(2:ncol(fichiermin2))],n=10)
# threshold of 480 (if data>480, replace by 480) in the 10 minutes file
fichier10minseuil<-as.data.frame(fichier10min[,1])
for(i in (2:97)){
  fichier10minseuil[,i]<-replace(fichier10min[,i],fichier10min[,i]>480,480)
}
# moyenne glissante
fichier10minseuilglisse<-as.data.frame(fichier10min[,1])
for(i in (2:97)){
  fichier10minseuilglisse[,i]<-frollmean(fichier10min[,i],10,align="center")
}

# 8. PLOTS ----
# 8.1 plot only one condition/genotype : 
# replace fichieraploter by your file (fichiermin, fichiermin2, fichier10min, fichier10min2...)
# replace genotype by your condition/genotype (good, ctrl, opn4xa...)
# replace "monTitre" by what you want your title to be
# replace maxdey by what you want your maximum in y to be (remove ylim = c(0,maxdey) if you don't want to specify this)
plot(rowMeans(fichieraploter[, (genotype+1)]), type = "l", main = "monTitre", ylim = c(0,maxdey))
lignestrans(t, n)  # to add transitions to the plot, replace n by the number of seconds of your file (for fichier10min replace by 600)

# 8.2 plot two conditions/genotypes
# replace fichieraploter by your file (fichiermin, fichiermin2, fichier10min, fichier10min2...)
# replace genotype1 and genotype 2 by your conditions/genotypes (good, ctrl, opn4xa...)
# replace "monTitre" by what you want your title to be
# replace maxdey by what you want your maximum in y to be (remove ylim = c(0,maxdey) if you don't want to specify this)
plot(rowMeans(fichieraploter[, (genotype1+1)]), type = "l", main = "monTitre", ylim = c(0,maxdey)) +
  lines(rowMeans(fichieraploter[, (genotype2+1)]), type = "l", col = "red")
lignestrans(t, n)  # to add transitions to the plot, replace n by the number of seconds of your file (for fichier10min replace by 600)

# 8.3 plot each arena
# replace genotype by your condition/genotype (good, ctrl, opn4xa...)
# replace fichieraploter by your file (fichiermin, fichier10min, ...)
for(i in (1:length(genotype))){
  plot((fichieraploter[,(genotype[i]+1)]),type="l",main=i)
}

# 8.4 plot during the light pulse transitions +/- 60 seconds
# plot only one condition/genotype
# replace transitiontoplot by the transition you want to plot (ltod  for light to dark, dtol for light to dark)
# replace genotype by your condition/genotype (good, ctrl, opn4xa...)
plot(rowMeans(fichier[transitiontoplot,(genotype+7)]),type="l")
abline(v=61,col="grey",lty=2)  # add a grey line at 61 seconds
# plot two conditions/genotypes
# replace transitiontoplot by the transition you want to plot (ltod  for light to dark, dtol for light to dark)
# replace genotype1 and genotype 2 by your conditions/genotypes (good, ctrl, opn4xa...)
plot(rowMeans(fichier[transitiontoplot,(genotype1+7)]),type="l")+
  lines(rowMeans(fichieraploter[transitiontoplot,(genotype2+7)]),type="l",col="red")
abline(v=61,col="grey",lty=2)  # add a grey line at 61 seconds

# 9. EXPORT FILES ----
# replace fichiertoexport to which file you want to export (fichiermin, fichier10min, etc...)
# replace "nametogive.xlsx" by the name you want to give to this file (don't forget the quotes)
write_xlsx(fichiertoexport, path = paste(dossier_sortie,"//","nametogive.xlsx",sep=""), col_names = TRUE)
