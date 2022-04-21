# export : 2000 x 800 (jpeg)
# 600 x 724


#### LIBRARIES ####
library(dplyr)
library(ggplot2)
library(plotly)#interactive plot
library(readxl)
library(behavr)
library(ggetho)
library(rstatix)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # palette colorblindfriendly


setwd("Y://ELISE DATA_17Gb a trier//EXCELL+stats//excell//zebrafish locomotor activity//Clair_bilans//juin2021//1.data//2.stats_mmdeb")

##### FUNCTIONS ######
  #1.transformer data pour ggplot // transform data for ggplot
transfo_data<-function(fichier, lignes, colwt, nwt, colmut, nmut, minutes,conditionwt,conditionmut){
    #fichier : fichier à transformer // file to transform
    #lignes : nombre de lignes du fichier à considérer // number of lines of the file to consider
    #colwt : numéro de la première colonne wt (les wt doivent être à la suite les uns des autres) // number of the first column containing wt (wt all needs to be next each other)
    #nwt : nombre de wt // number of wt
    #colmut : numéro de la première colonne mut (les mut doivent être à la suite les uns des autres) // number of the first column containing mutants (mutants all needs to be next each other)
    #nmut : nombre de mut // number of wt
    #minutes : intervalle de temps en minutes (DM10, DM60 ?) // time bin in minutes
    #conditionwt : nom de la condition wt (controle, wt, non pigmente) // name of the wild type condition
    #conditionmut: nom de la condition mut(mutant, opn4xa-/-, pigmente, etc) // name of the mutant condition
  wt<-print(paste(nwt,conditionwt,sep=" "))
  mut<-print(paste(nmut,conditionmut,sep=" "))
  data<-as.data.frame(cbind(time=seq(from=0,length.out=lignes,by=mins(minutes)), arena=colnames(fichier[colwt]), genotype=wt, distance=unlist(fichier[c(1:lignes),colwt],use.names=FALSE)), stringsAsFactors=FALSE)
  for (i in 1:(nwt-1)){
    n<-i+colwt
    data<-as.data.frame(rbind(data, as.data.frame(cbind(time=seq(from=0,length.out=lignes,by=mins(minutes)), arena=colnames(fichier[n]), genotype=wt, distance=unlist(fichier[c(1:lignes),n], use.names=FALSE)))),stringsAsFactors=FALSE)
  }
  
  for (i in 0:(nmut-1)){
    n<-i+colmut
    data<-as.data.frame(rbind(data, as.data.frame(cbind(time=seq(from=0,length.out=lignes,by=mins(minutes)), arena=colnames(fichier[n]), genotype=mut, distance=unlist(fichier[c(1:lignes),n], use.names=FALSE)))),stringsAsFactors=FALSE)
  }
  return(data)
}
  #2. créer fichier metadata pour rethomics // create metadata file for rethomics
create_meta<-function(data,nwt,nmut,conditionwt,conditionmut){
  wt<-print(paste(nwt,conditionwt,sep=" "))
  mut<-print(paste(nmut,conditionmut,sep=" "))
  data_meta<-data.table::data.table(arena=unique(data$arena), genotype=c(rep(wt,nwt),rep(mut,nmut)),key="arena")
  return(data_meta)
}
  #3. créer fichier behavr pour rethomics // create behavr file for rethomics
create_behavr<-function(data,data_meta){
  data_behavr<-behavr(as.data.table(data,key="arena"),data_meta)
  return(data_behavr)
}
  #4. créer le plot avec rethomics // create plot with rethomics (by default error bars are standard errors) 
plot_ggetho<-function(data,titre,phaseLD,couleurs,couleurwt="#0072B2",couleurmutant="#E69F00"){
    #data : fichier behavr // behavr file
    #nwt : nombre de wt // number of wt
    #nmut : nombre de mutants // number of mutants
    #titre : titre (entre guillemets) // title (between " ")
    #phaseLD : en secondes, temps avant la premiere transition nuit/jour // in seconds, time before the first night/dark transition
    #couleurs : vecteur qui contient les couleurs du bandeau (c("grey","black")) par exemple // vector containing colors of the horizontal bands
  intervalles<-c(phaseLD,phaseLD+hours(14),phaseLD+hours(14)+hours(10),phaseLD+hours(14)+hours(10)+hours(14),phaseLD+hours(14)+hours(10)+hours(14)+hours(10),phaseLD+hours(14)+hours(10)+hours(14)+hours(10)+hours(14),phaseLD+hours(14)+hours(10)+hours(14)+hours(10)+hours(14)+hours(10),phaseLD+hours(14)+hours(10)+hours(14)+hours(10)+hours(14)+hours(10)+hours(14),phaseLD+hours(14)+hours(10)+hours(14)+hours(10)+hours(14)+hours(10)+hours(14)+hours(10))
  p<-ggetho(data,aes(x=as.numeric(time),y=as.numeric(distance),color=factor(genotype,levels=c(wt,mut))))+
    scale_color_manual(values=c(couleurwt,couleurmutant))+
    scale_fill_manual(values=c(couleurwt,couleurmutant))+
    stat_pop_etho()+
    labs(title=titre, x="Time", y="Mean distance travelled (mm/min over 10min)")+
    scale_x_discrete(expand=c(0,0))+
    geom_vline(xintercept=intervalles,color="darkgrey",linetype="dashed",lwd=1.4)+
    stat_ld_annotations(l_duration=hours(14), phase=(phaseLD),period=hours(24),ld_colours=couleurs)+
    theme_minimal(base_line_size=1.4,base_rect_size=1)+theme(legend.title=element_blank())+geom_vline(xintercept=0)
  return(p)
}

  #5. plot period mean +/- sd
dataperiod <- function(data, col1, col2){  # create data frame for the plot
  df <- as.data.frame(cbind(period=c(unlist(data[,col1]), unlist(data[,col2])), genotype=c(rep("wt", nrow(data)), rep("mut", nrow(data)))),stringsAsFactors = FALSE)
  df <- as.data.frame(cbind(period=df$period[!is.na(df$period)], genotype=df$genotype[!is.na(df$period)]))
  df$genotype <- factor(df$genotype,levels=c("wt","mut"))
  df$period <- as.numeric(df$period)
  return(df)
}
periodplot<-function(data,titre,colwt="#0072B2",colmut="#E69F00",wt="wt",mut="mut"){
  df.summary <- data %>%
    group_by(genotype) %>%
    summarise(
      sd = sd(period),
      period = mean(period)
    )
  df.summary
  p<-ggplot(data, aes(genotype, as.numeric(period),col=genotype)) +
    geom_jitter(
      position = position_jitter(0.2), color = "darkgray",size=4, alpha=0.4
    ) + 
    geom_pointrange(
      aes(ymin = (period-sd), ymax = (period+sd)),size = 2.5,shape=15,alpha=0.8,data = df.summary
    )+ylim(20,30)+scale_x_discrete(limits=c(wt,mut))+
    labs(title=titre, y="period")+scale_color_manual(values=c(colwt,colmut))+theme_light(base_line_size=2.5,base_rect_size = 1.4)+scale_x_discrete(breaks = NULL)
  return(p)
}



##### 1. COMPORTEMENT OPN4XA #####

  ### LD (48): M6(11), M7(21), M8(18)###
#1. tranformer data pour ggplot et rethomics
LD_opn4xa<-read_excel("LDopn4xa_M6-M7-M8.xlsx")
LD_opn4xa_data<-transfo_data(LD_opn4xa,lignes=438,colwt=3,nwt=48,colmut=51,nmut=48, minutes = 10,conditionwt="opn4xa+/+",conditionmut="opn4xa-/-")
#2. plot rethomics
LD_opn4xa_meta<-create_meta(LD_opn4xa_data,48,48,conditionwt="opn4xa+/+",conditionmut="opn4xa-/-")
LD_opn4xa_behavr<-create_behavr(LD_opn4xa_data,LD_opn4xa_meta)
LD_opn4xa_behavr$distance[is.na(LD_opn4xa_behavr$distance)]<-0
wt<-print(paste(48,"opn4xa+/+",sep=" "))
mut<-print(paste(48,"opn4xa-/-",sep=" "))
plot_ggetho(LD_opn4xa_behavr,titre="LD opn4xa",phaseLD=3000,couleurs=c("white","black"))

  ### DD (65) : DD1(18), DD2(22), DD3(25) 
DD_opn4xa<-read_excel("DDopn4xa_DD1-DD2-DD3.xlsx")
DD_opn4xa_data<-transfo_data(DD_opn4xa,lignes=583,colwt=3,nwt=65,colmut=68,nmut=65, minutes = 10,conditionwt="opn4xa+/+",conditionmut="opn4xa-/-")
DD_opn4xa_meta<-create_meta(DD_opn4xa_data,65,65,conditionwt="opn4xa+/+",conditionmut="opn4xa-/-")
DD_opn4xa_behavr<-create_behavr(DD_opn4xa_data,DD_opn4xa_meta)
DD_opn4xa_behavr$distance[is.na(DD_opn4xa_behavr$distance)]<-0
wt<-print(paste(65,"opn4xa+/+",sep=" "))
mut<-print(paste(65,"opn4xa-/-",sep=" "))
plot_ggetho(DD_opn4xa_behavr,titre="DD opn4xa",phaseLD=1800, couleurs=c("grey","black"))

  ### LL (66) : LL5(17), LL8(14), LL9(19), LL11(16) ###
LL_opn4xa<-read_excel("LLopn4xa_LL5-LL8-LL9-LL11.xlsx")
LL_opn4xa_data<-transfo_data(LL_opn4xa,lignes=436,colwt=3,nwt=66,colmut=69,nmut=66, minutes = 10,conditionwt="opn4xa+/+",conditionmut="opn4xa-/-")
LL_opn4xa_meta<-create_meta(LL_opn4xa_data,66,66,conditionwt="opn4xa+/+",conditionmut="opn4xa-/-")
LL_opn4xa_behavr<-create_behavr(LL_opn4xa_data,LL_opn4xa_meta)
LL_opn4xa_behavr$distance[is.na(LL_opn4xa_behavr$distance)]<-0
wt<-print(paste(66,"opn4xa+/+",sep=" "))
mut<-print(paste(66,"opn4xa-/-",sep=" "))
plot_ggetho(LL_opn4xa_behavr,titre="LL opn4xa",phaseLD=1800, couleurs=c("white","grey92"))


  ### PS (58):  PS1(12), PS2(23), PS5(23) ###
PS_opn4xa<-read_excel("PSopn4xa_PS1-PS2-PS5.xlsx")
PS_opn4xa_data<-transfo_data(PS_opn4xa,lignes=583,colwt=3,nwt=58,colmut=61,nmut=58, minutes = 10,conditionwt="opn4xa+/+",conditionmut="opn4xa-/-")
PS_opn4xa_meta<-create_meta(PS_opn4xa_data,58,58,conditionwt="opn4xa+/+",conditionmut="opn4xa-/-")
PS_opn4xa_behavr<-create_behavr(PS_opn4xa_data,PS_opn4xa_meta)
PS_opn4xa_behavr$distance[is.na(PS_opn4xa_behavr$distance)]<-0
wt<-print(paste(58,"opn4xa+/+",sep=" "))
mut<-print(paste(58,"opn4xa-/-",sep=" "))
phaseshift<-1800+hours(14)+hours(10)+hours(14)+hours(2)
plot_ggetho(PS_opn4xa_behavr,titre="PS opn4xa",phaseLD=1800, couleurs=c("grey","black"))+stat_ld_annotations(x_limits=(c(phaseshift,(phaseshift+hours(2)))),ld_colours=c("white","white"))+
  geom_vline(xintercept=c((1800+hours(14)+hours(10)+hours(14)+hours(2)),(1800+hours(14)+hours(10)+hours(14)+hours(2))+hours(2)),colour="deeppink",linetype="dashed",lwd=1.4,alpha=0.8)


##### 2. COMPORTEMENT LAKRITZ #####

  ### LD (55): LD2(25), LD3(15), LD4(15) ###
LD_lak<-read_excel("LDlak_LD2-LD3-LD4.xlsx")
LD_lak_data<-transfo_data(LD_lak,lignes=438,colwt=3,nwt=55,colmut=58,nmut=55, minutes = 10,conditionwt="ctrl",conditionmut="lak-/-")
LD_lak_meta<-create_meta(LD_lak_data,55,55,conditionwt="ctrl",conditionmut="lak-/-")
LD_lak_behavr<-create_behavr(LD_lak_data,LD_lak_meta)
LD_lak_behavr$distance[is.na(LD_lak_behavr$distance)]<-0
wt<-print(paste(55,"ctrl",sep=" "))
mut<-print(paste(55,"lak-/-",sep=" "))
plot_ggetho(LD_lak_behavr,titre="LD lak",phaseLD=3000, couleurs=c("white","black"),couleurmut="#CC79A7")

  ### DD (48):  DD1(18), DD2(30) ###
DD_lak<-read_excel("DDlak_DD1-DD2.xlsx")
DD_lak_data<-transfo_data(DD_lak,lignes=583,colwt=3,nwt=48,colmut=51,nmut=48,minutes=10,conditionwt="ctrl",conditionmut="lak-/-")
DD_lak_meta<-create_meta(DD_lak_data,48,48,conditionwt="ctrl",conditionmut="lak-/-")
DD_lak_behavr<-create_behavr(DD_lak_data,DD_lak_meta)
DD_lak_behavr$distance[is.na(DD_lak_behavr$distance)]<-0
wt<-print(paste(48,"ctrl",sep=" "))
mut<-print(paste(48,"lak-/-",sep=" "))
plot_ggetho(DD_lak_behavr,titre="DD lakritz",phaseLD=1800, couleurs=c("grey","black"),couleurmut="#CC79A7")

  ### LL(81): LL2(12), LL3(32),LL4(37)  ###
LL_lak<-read_excel("LLlak_LL2-LL3-LL4.xlsx")
LL_lak_data<-transfo_data(LL_lak,lignes=436,colwt=3,nwt=81,colmut=84,nmut=81, minutes = 10,conditionwt="ctrl",conditionmut="lak-/-")
LL_lak_meta<-create_meta(LL_lak_data,81,81,conditionwt="ctrl",conditionmut="lak-/-")
LL_lak_behavr<-create_behavr(LL_lak_data,LL_lak_meta)
LL_lak_behavr$distance[is.na(LL_lak_behavr$distance)]<-0
wt<-print(paste(81,"ctrl",sep=" "))
mut<-print(paste(81,"lak-/-",sep=" "))
plot_ggetho(LL_lak_behavr,titre="LL lak",phaseLD=1800, couleurs=c("white","grey92"),couleurmut="#CC79A7")

  ### PS (62): PS1(12), PS2(25), PS4(25)###
PS_lak<-read_excel("PSlak_PS1-PS2-PS4.xlsx")
PS_lak_data<-transfo_data(PS_lak,lignes=583,colwt=3,nwt=62,colmut=65,nmut=62, minutes = 10,conditionwt="ctrl",conditionmut="lak-/-")
PS_lak_meta<-create_meta(PS_lak_data,62,62,conditionwt="ctrl",conditionmut="lak-/-")
PS_lak_behavr<-create_behavr(PS_lak_data,PS_lak_meta)
PS_lak_behavr$distance[is.na(PS_lak_behavr$distance)]<-0
wt<-print(paste(62,"ctrl",sep=" "))
mut<-print(paste(62,"lak-/-",sep=" "))
phaseshift<-1800+hours(14)+hours(10)+hours(14)+hours(2)
plot_ggetho(PS_lak_behavr,titre="PS1 lak",phaseLD=1800, couleurs=c("grey","black"),couleurmut="#CC79A7")+stat_ld_annotations(x_limits=(c(phaseshift,(phaseshift+hours(2)))),ld_colours=c("white","white"))+
  geom_vline(xintercept=c((1800+hours(14)+hours(10)+hours(14)+hours(2)),(1800+hours(14)+hours(10)+hours(14)+hours(2))+hours(2)),colour="deeppink",linetype="dashed",lwd=1.4,alpha=0.8)

##### 3. COMPORTEMENT PHASE SHIFT WT #####

#WT opn4xa 
PS_wt<-read_excel("DD_PS_opn4xa.xlsx")
PS_wt_data<-transfo_data(PS_wt,lignes=583,colwt=3,nwt=65,colmut=68,nmut=58, minutes = 10,conditionwt="ctrl",conditionmut="phase_shift")
PS_wt_meta<-create_meta(PS_wt_data,65,58,conditionwt="ctrl",conditionmut="phase_shift")
PS_wt_behavr<-create_behavr(PS_wt_data,PS_wt_meta)
PS_wt_behavr$distance[is.na(PS_wt_behavr$distance)]<-0
wt<-print(paste(65,"ctrl",sep=" "))
mut<-print(paste(58,"phase_shift",sep=" "))
phaseshift<-1800+hours(14)+hours(10)+hours(14)+hours(2)
plot_ggetho(PS_wt_behavr,titre="PS wt opn4xa",phaseLD=1800, couleurs=c("grey","black"),couleurmut="#56B4E9")+stat_ld_annotations(x_limits=(c(phaseshift,(phaseshift+hours(2)))),ld_colours=c("white","white"))+
  geom_vline(xintercept=c((1800+hours(14)+hours(10)+hours(14)+hours(2)),(1800+hours(14)+hours(10)+hours(14)+hours(2))+hours(2)),colour="deeppink",linetype="dashed",lwd=1.4,alpha=0.8)

#WT lakritz 
PS_wt<-read_excel("DD_PS_lak.xlsx")
PS_wt_data<-transfo_data(PS_wt,lignes=583,colwt=3,nwt=48,colmut=51,nmut=62, minutes = 10,conditionwt="ctrl",conditionmut="phase_shift")
PS_wt_meta<-create_meta(PS_wt_data,48,62,conditionwt="ctrl",conditionmut="phase_shift")
PS_wt_behavr<-create_behavr(PS_wt_data,PS_wt_meta)
PS_wt_behavr$distance[is.na(PS_wt_behavr$distance)]<-0
wt<-print(paste(48,"ctrl",sep=" "))
mut<-print(paste(62,"phase_shift",sep=" "))
phaseshift<-1800+hours(14)+hours(10)+hours(14)+hours(2)
plot_ggetho(PS_wt_behavr,titre="PS wt lakrtiz",phaseLD=1800, couleurs=c("grey","black"),couleurmut="#56B4E9")+stat_ld_annotations(x_limits=(c(phaseshift,(phaseshift+hours(2)))),ld_colours=c("white","white"))+
  geom_vline(xintercept=c((1800+hours(14)+hours(10)+hours(14)+hours(2)),(1800+hours(14)+hours(10)+hours(14)+hours(2))+hours(2)),colour="deeppink",linetype="dashed",lwd=1.4,alpha=0.8)


##### 4. COMPORTEMENT DOUBLE #####

### PS (12):  PS1(6), PS3 (9) ###
PS_double<-read_excel("PSdouble_1-2-3.xlsx")
PS_double_data<-transfo_data(PS_double,lignes=575,colwt=3,nwt=27,colmut=30,nmut=27, minutes = 10,conditionwt="lak-/-",conditionmut="double")
PS_double_meta<-create_meta(PS_double_data,27,27,conditionwt="lak-/-",conditionmut="double")
PS_double_behavr<-create_behavr(PS_double_data,PS_double_meta)
PS_double_behavr$distance[is.na(PS_double_behavr$distance)]<-0
wt<-print(paste(27,"lak-/-",sep=" "))
mut<-print(paste(27,"double",sep=" "))
phaseshift<-1800+hours(14)+hours(10)+hours(14)+hours(2)
plot_ggetho(PS_double_behavr,titre="PS double",phaseLD=1800, couleurs=c("grey","black"),couleurwt="#CC79A7", couleurmut="#999999")+stat_ld_annotations(x_limits=(c(phaseshift,(phaseshift+hours(2)))),ld_colours=c("white","white"))+
  geom_vline(xintercept=c((1800+hours(14)+hours(10)+hours(14)+hours(2)),(1800+hours(14)+hours(10)+hours(14)+hours(2))+hours(2)),colour="deeppink",linetype="dashed",lwd=1.4,alpha=0.8)


##### 5. PERIODS #####
periods <- read_excel("period_detrend.xlsx")

periods_LLopn4xa <- dataperiod(periods, col1=1, col2=2)
periods_LLlak <- dataperiod(periods, col1=3, col2=4)
periods_DDopn4xa <- dataperiod(periods, col1=5, col2=6)
periods_DDlak <- dataperiod(periods, col1=7, col2=8)

periodplot(periods_LLopn4xa, "LL opn4xa")
periodplot(periods_DDopn4xa, "DD opn4xa")
periodplot(periods_DDlak, "DD lak",colmut="#CC79A7")
periodplot(periods_LLlak, "LL lak",colmut="#CC79A7")

wilcox.test(as.numeric(period)~genotype, data=periods_LLopn4xa)
wilcox.test(as.numeric(period)~genotype, data=periods_DDopn4xa)
wilcox.test(as.numeric(period)~genotype, data=periods_LLlak)
wilcox.test(as.numeric(period)~genotype, data=periods_DDlak)


# comparaison LL opn4xa

df <- as.data.frame(cbind(period=c(unlist(periods[,9]), unlist(periods[,10]), unlist(periods[,11]), unlist(periods[,12]), unlist(periods[,13]), unlist(periods[,14]), unlist(periods[,15]), unlist(periods[,16])), genotype=c(rep("LL5wt", nrow(periods)), rep("LL5mut", nrow(periods)), rep("LL8wt", nrow(periods)), rep("LL8mut", nrow(periods)), rep("LL9wt", nrow(periods)), rep("LL9mut", nrow(periods)), rep("LL11wt", nrow(periods)), rep("LL11mut", nrow(periods)))),stringsAsFactors = FALSE)
df <- as.data.frame(cbind(period=df$period[!is.na(df$period)], genotype=df$genotype[!is.na(df$period)]))
df$period <- as.numeric(df$period)
df$genotype <- factor(df$genotype, levels=c("LL5wt", "LL5mut", "LL8wt", "LL8mut", "LL9wt", "LL9mut", "LL11wt", "LL11mut"))


df.summary <- df %>%
  group_by(genotype) %>%
  summarise(
    premquartile = summary(period)[2],
    troisquartile = summary(period)[5],
    period = median(period)
  )
df.summary
p<-ggplot(df, aes(genotype, as.numeric(period),col=genotype)) +
  geom_jitter(
    position = position_jitter(0.2), color = "darkgray",size=4, alpha=0.4
  ) + 
  geom_pointrange(
    aes(ymin = premquartile, ymax = troisquartile),size = 2.5,shape=15,alpha=0.8,data = df.summary
  )+ylim(20,30)+scale_x_discrete(limits=c("LL5wt", "LL5mut", "LL8wt", "LL8mut", "LL9wt", "LL9mut", "LL11wt", "LL11mut"))+
  labs(title="comparaison LL opn4xa", y="period")+theme_light(base_line_size=2.5,base_rect_size = 1.4)+scale_x_discrete(breaks = NULL)
p


##### 7. DDT ####
DDT<-read_excel("DDT_comparaison.xlsx")

DDT_data_opn4xaDD<-transfo_data(DDT,lignes=583,colwt=4,nwt=65,colmut=239,nmut=59, minutes = 10,conditionwt="DDopn4xa",conditionmut="DDT")
DDT_data_opn4xaPS<-transfo_data(DDT,lignes=583,colwt=69,nwt=58,colmut=239,nmut=59, minutes = 10,conditionwt="PSopn4xa",conditionmut="DDT")
DDT_data_lakDD<-transfo_data(DDT,lignes=583,colwt=128,nwt=48,colmut=239,nmut=59, minutes = 10,conditionwt="DDlak",conditionmut="DDT")
DDT_data_lakPS<-transfo_data(DDT,lignes=583,colwt=176,nwt=62,colmut=239,nmut=59, minutes = 10,conditionwt="PSlak",conditionmut="DDT")

DDT_opn4xaDD_meta<-create_meta(DDT_data_opn4xaDD,65,59,conditionwt="DDopn4xa",conditionmut="DDT")
DDT_opn4xaPS_meta<-create_meta(DDT_data_opn4xaPS,58,59,conditionwt="PSopn4xa",conditionmut="DDT")
DDT_lakDD_meta<-create_meta(DDT_data_lakDD,48,59,conditionwt="DDlak",conditionmut="DDT")
DDT_lakPS_meta<-create_meta(DDT_data_lakPS,62,59,conditionwt="PSlak",conditionmut="DDT")

DDT_opn4xaDD_behavr<-create_behavr(DDT_data_opn4xaDD,DDT_opn4xaDD_meta)
DDT_opn4xaDD_behavr$distance[is.na(DDT_opn4xaDD_behavr$distance)]<-0
DDT_opn4xaPS_behavr<-create_behavr(DDT_data_opn4xaPS,DDT_opn4xaPS_meta)
DDT_opn4xaPS_behavr$distance[is.na(DDT_opn4xaPS_behavr$distance)]<-0
DDT_lakDD_behavr<-create_behavr(DDT_data_lakDD,DDT_lakDD_meta)
DDT_lakDD_behavr$distance[is.na(DDT_lakDD_behavr$distance)]<-0
DDT_lakPS_behavr<-create_behavr(DDT_data_lakPS,DDT_lakPS_meta)
DDT_lakPS_behavr$distance[is.na(DDT_lakPS_behavr$distance)]<-0

# DDopn4xa vs DDT
wt<-print(paste(65,"DDopn4xa",sep=" "))
mut<-print(paste(59,"DDT",sep=" "))
phaseshift<-1800+hours(14)+hours(10)+hours(14)+hours(2)
plot_ggetho(DDT_opn4xaDD_behavr,titre="DDopn4xa vs DDT",phaseLD=1800, couleurs=c("grey","black"),couleurmut="deeppink")+stat_ld_annotations(x_limits=(c(phaseshift,(phaseshift+hours(2)))),ld_colours=c("white","black"))+
  geom_vline(xintercept=c((1800+hours(14)+hours(10)+hours(14)+hours(2)),(1800+hours(14)+hours(10)+hours(14)+hours(2))+hours(2)),colour="deeppink",linetype="dashed",lwd=1.4,alpha=0.8)
# PS opn4xa vs DDT
wt<-print(paste(58,"PSopn4xa",sep=" "))
mut<-print(paste(59,"DDT",sep=" "))
phaseshift<-1800+hours(14)+hours(10)+hours(14)+hours(2)
plot_ggetho(DDT_opn4xaPS_behavr,titre="PSopn4xa vs DDT",phaseLD=1800, couleurs=c("grey","black"), couleurwt="#56B4E9", couleurmut="deeppink")+stat_ld_annotations(x_limits=(c(phaseshift,(phaseshift+hours(2)))),ld_colours=c("white","black"))+
  geom_vline(xintercept=c((1800+hours(14)+hours(10)+hours(14)+hours(2)),(1800+hours(14)+hours(10)+hours(14)+hours(2))+hours(2)),colour="deeppink",linetype="dashed",lwd=1.4,alpha=0.8)

# DDlak vs DDT
wt<-print(paste(48,"DDlak",sep=" "))
mut<-print(paste(59,"DDT",sep=" "))
phaseshift<-1800+hours(14)+hours(10)+hours(14)+hours(2)
plot_ggetho(DDT_lakDD_behavr,titre="DDlak vs DDT",phaseLD=1800, couleurs=c("grey","black"), couleurmut="deeppink")+stat_ld_annotations(x_limits=(c(phaseshift,(phaseshift+hours(2)))),ld_colours=c("white","black"))+
  geom_vline(xintercept=c((1800+hours(14)+hours(10)+hours(14)+hours(2)),(1800+hours(14)+hours(10)+hours(14)+hours(2))+hours(2)),colour="deeppink",linetype="dashed",lwd=1.4,alpha=0.8)
# PSlak vs DDT
wt<-print(paste(62,"PSlak",sep=" "))
mut<-print(paste(59,"DDT",sep=" "))
phaseshift<-1800+hours(14)+hours(10)+hours(14)+hours(2)
plot_ggetho(DDT_lakPS_behavr,titre="PSlak vs DDT",phaseLD=1800, couleurs=c("grey","black"), couleurwt="#56B4E9", couleurmut="deeppink")+stat_ld_annotations(x_limits=(c(phaseshift,(phaseshift+hours(2)))),ld_colours=c("white","black"))+
  geom_vline(xintercept=c((1800+hours(14)+hours(10)+hours(14)+hours(2)),(1800+hours(14)+hours(10)+hours(14)+hours(2))+hours(2)),colour="deeppink",linetype="dashed",lwd=1.4,alpha=0.8)



#### 8. GRAPHES COMPARAISONS EXPERIENCE LL ####

### OPN4XA ###
# LL2 opn4xa (12wt, 20mut, attention NO MDM)
LL2 <- read_excel("DM10_LL2_div10.xlsx")
LL2_data <- transfo_data(LL2, lignes=436, colwt=55, nwt=12, colmut=35, nmut=20, minutes=10, conditionwt="opn4xa+/+", conditionmut="opn4xa-/-")
LL2_meta <- create_meta(LL2_data, 12, 20, conditionwt="opn4xa+/+", conditionmut="opn4xa-/-")
LL2_behavr <- create_behavr(LL2_data, LL2_meta)
LL2_behavr$distance[is.na(LL2_behavr$distance)]<-0
wt<-print(paste(12,"opn4xa+/+",sep=" "))
mut<-print(paste(20,"opn4xa-/-",sep=" "))
plot_ggetho(LL2_behavr,titre="LL2 opn4xa",phaseLD=1800, couleurs=c("white","grey92"))

# LL opn4xa par manip : LL5(17), LL8(14), LL9(19), LL11(16)
LL_opn4xa<-read_excel("opn4xa_LL_LD.xlsx")
#LL5
LL5_data<-transfo_data(LL_opn4xa,lignes=436,colwt=3,nwt=17,colmut=69,nmut=17, minutes = 10,conditionwt="opn4xa+/+",conditionmut="opn4xa-/-")
LL5_meta <- create_meta(LL5_data, 17, 17, conditionwt="opn4xa+/+", conditionmut="opn4xa-/-")
LL5_behavr <- create_behavr(LL5_data, LL5_meta)
LL5_behavr$distance[is.na(LL5_behavr$distance)]<-0
wt<-print(paste(17,"opn4xa+/+",sep=" "))
mut<-print(paste(17,"opn4xa-/-",sep=" "))
plot_ggetho(LL5_behavr,titre="LL5 opn4xa",phaseLD=1800, couleurs=c("white","grey92"))
#LL8
LL8_data<-transfo_data(LL_opn4xa,lignes=436,colwt=20,nwt=14,colmut=86,nmut=14, minutes = 10,conditionwt="opn4xa+/+",conditionmut="opn4xa-/-")
LL8_meta <- create_meta(LL8_data, 14, 14, conditionwt="opn4xa+/+", conditionmut="opn4xa-/-")
LL8_behavr <- create_behavr(LL8_data, LL8_meta)
LL8_behavr$distance[is.na(LL8_behavr$distance)]<-0
wt<-print(paste(14,"opn4xa+/+",sep=" "))
mut<-print(paste(14,"opn4xa-/-",sep=" "))
plot_ggetho(LL8_behavr,titre="LL8 opn4xa",phaseLD=1800, couleurs=c("white","grey92"))
#LL9
LL9_data<-transfo_data(LL_opn4xa,lignes=436,colwt=34,nwt=19,colmut=100,nmut=19, minutes = 10,conditionwt="opn4xa+/+",conditionmut="opn4xa-/-")
LL9_meta <- create_meta(LL9_data, 19, 19, conditionwt="opn4xa+/+", conditionmut="opn4xa-/-")
LL9_behavr <- create_behavr(LL9_data, LL9_meta)
LL9_behavr$distance[is.na(LL9_behavr$distance)]<-0
wt<-print(paste(19,"opn4xa+/+",sep=" "))
mut<-print(paste(19,"opn4xa-/-",sep=" "))
plot_ggetho(LL9_behavr,titre="LL9 opn4xa",phaseLD=1800, couleurs=c("white","grey92"))
#LL11
LL11_data<-transfo_data(LL_opn4xa,lignes=436,colwt=53,nwt=16,colmut=119,nmut=16, minutes = 10,conditionwt="opn4xa+/+",conditionmut="opn4xa-/-")
LL11_meta <- create_meta(LL11_data, 16, 16, conditionwt="opn4xa+/+", conditionmut="opn4xa-/-")
LL11_behavr <- create_behavr(LL11_data, LL11_meta)
LL11_behavr$distance[is.na(LL11_behavr$distance)]<-0
wt<-print(paste(16,"opn4xa+/+",sep=" "))
mut<-print(paste(16,"opn4xa-/-",sep=" "))
plot_ggetho(LL11_behavr,titre="LL11 opn4xa",phaseLD=1800, couleurs=c("white","grey92"))

# LDwt vs LLwt
LDwt_LLwt_data<-transfo_data(LL_opn4xa,lignes=436,colwt=137,nwt=48,colmut=3,nmut=66, minutes = 10,conditionwt="LD opn4xa+/+",conditionmut="LL opn4xa+/+")
LDwt_LLwt_meta <- create_meta(LDwt_LLwt_data, 48, 66, conditionwt="LD opn4xa+/+", conditionmut="LL opn4xa+/+")
LDwt_LLwt_behavr <- create_behavr(LDwt_LLwt_data, LDwt_LLwt_meta)
LDwt_LLwt_behavr$distance[is.na(LDwt_LLwt_behavr$distance)]<-0
wt<-print(paste(48,"LD opn4xa+/+",sep=" "))
mut<-print(paste(66,"LL opn4xa+/+",sep=" "))
plot_ggetho(LDwt_LLwt_behavr,titre="LDwt vs LLwt opn4xa",phaseLD=1800, couleurs=c("white","grey92"), couleurmut="#56B4E9")
# LDmut vs LLwt
LDwt_LLwt_data<-transfo_data(LL_opn4xa,lignes=436,colwt=185,nwt=48,colmut=3,nmut=66, minutes = 10,conditionwt="LD opn4xa-/-",conditionmut="LL opn4xa+/+")
LDwt_LLwt_meta <- create_meta(LDwt_LLwt_data, 48, 66, conditionwt="LD opn4xa-/-", conditionmut="LL opn4xa+/+")
LDwt_LLwt_behavr <- create_behavr(LDwt_LLwt_data, LDwt_LLwt_meta)
LDwt_LLwt_behavr$distance[is.na(LDwt_LLwt_behavr$distance)]<-0
wt<-print(paste(48,"LD opn4xa-/-",sep=" "))
mut<-print(paste(66,"LL opn4xa+/+",sep=" "))
plot_ggetho(LDwt_LLwt_behavr,titre="LDmut vs LLwt opn4xa",phaseLD=1800, couleurs=c("white","grey92"), couleurmut="#56B4E9", couleurwt="#E69F00")
# LDwt vs LLmut
LDwt_LLwt_data<-transfo_data(LL_opn4xa,lignes=436,colwt=137,nwt=48,colmut=69,nmut=66, minutes = 10,conditionwt="LD opn4xa+/+",conditionmut="LL opn4xa-/-")
LDwt_LLwt_meta <- create_meta(LDwt_LLwt_data, 48, 66, conditionwt="LD opn4xa+/+", conditionmut="LL opn4xa-/-")
LDwt_LLwt_behavr <- create_behavr(LDwt_LLwt_data, LDwt_LLwt_meta)
LDwt_LLwt_behavr$distance[is.na(LDwt_LLwt_behavr$distance)]<-0
wt<-print(paste(48,"LD opn4xa+/+",sep=" "))
mut<-print(paste(66,"LL opn4xa-/-",sep=" "))
plot_ggetho(LDwt_LLwt_behavr,titre="LDwt vs LLmut opn4xa",phaseLD=1800, couleurs=c("white","grey92"), couleurmut="#F0E442")
# LDmut vs LLmut
LDwt_LLwt_data<-transfo_data(LL_opn4xa,lignes=436,colwt=185,nwt=48,colmut=69,nmut=66, minutes = 10,conditionwt="LD opn4xa-/-",conditionmut="LL opn4xa-/-")
LDwt_LLwt_meta <- create_meta(LDwt_LLwt_data, 48, 66, conditionwt="LD opn4xa-/-", conditionmut="LL opn4xa-/-")
LDwt_LLwt_behavr <- create_behavr(LDwt_LLwt_data, LDwt_LLwt_meta)
LDwt_LLwt_behavr$distance[is.na(LDwt_LLwt_behavr$distance)]<-0
wt<-print(paste(48,"LD opn4xa-/-",sep=" "))
mut<-print(paste(66,"LL opn4xa-/-",sep=" "))
plot_ggetho(LDwt_LLwt_behavr,titre="LDmut vs LLmut opn4xa",phaseLD=1800, couleurs=c("white","grey92"), couleurmut="#F0E442", couleurwt="#E69F00")

### LAKRITZ ###
LL_lak <- read_excel("lak_LL_LD.xlsx")

# LL lak par manip : LL2(12), LL3(32),LL4(37)
#LL2
LL2_data <- transfo_data(LL_lak, lignes=436, colwt=3, nwt=12, colmut=84, nmut=12, minutes=10, conditionwt="ctrl", conditionmut="lak-/-")
LL2_meta <- create_meta(LL2_data, 12, 12, conditionwt="ctrl", conditionmut="lak-/-")
LL2_behavr <- create_behavr(LL2_data, LL2_meta)
LL2_behavr$distance[is.na(LL2_behavr$distance)]<-0
wt<-print(paste(12,"ctrl",sep=" "))
mut<-print(paste(12,"lak-/-",sep=" "))
plot_ggetho(LL2_behavr,titre="LL2 lak",phaseLD=1800, couleurs=c("white","grey92"),couleurmut="#CC79A7")
#LL3
LL3_data <- transfo_data(LL_lak, lignes=436, colwt=15, nwt=32, colmut=96, nmut=32, minutes=10, conditionwt="ctrl", conditionmut="lak-/-")
LL3_meta <- create_meta(LL3_data, 32, 32, conditionwt="ctrl", conditionmut="lak-/-")
LL3_behavr <- create_behavr(LL3_data, LL3_meta)
LL3_behavr$distance[is.na(LL3_behavr$distance)]<-0
wt<-print(paste(32,"ctrl",sep=" "))
mut<-print(paste(32,"lak-/-",sep=" "))
plot_ggetho(LL3_behavr,titre="LL3 lak",phaseLD=1800, couleurs=c("white","grey92"),couleurmut="#CC79A7")
#LL4
LL4_data <- transfo_data(LL_lak, lignes=436, colwt=47, nwt=37, colmut=128, nmut=37, minutes=10, conditionwt="ctrl", conditionmut="lak-/-")
LL4_meta <- create_meta(LL4_data, 37, 37, conditionwt="ctrl", conditionmut="lak-/-")
LL4_behavr <- create_behavr(LL4_data, LL4_meta)
LL4_behavr$distance[is.na(LL4_behavr$distance)]<-0
wt<-print(paste(37,"ctrl",sep=" "))
mut<-print(paste(37,"lak-/-",sep=" "))
plot_ggetho(LL4_behavr,titre="LL4 lak",phaseLD=1800, couleurs=c("white","grey92"),couleurmut="#CC79A7")

# LD lak par manip : LD2(25), LD3(15), LD4(15)
#LD2
LD2_data <- transfo_data(LL_lak, lignes=436, colwt=167, nwt=25, colmut=222, nmut=25, minutes=10, conditionwt="ctrl", conditionmut="lak-/-")
LD2_meta <- create_meta(LD2_data, 25, 25, conditionwt="ctrl", conditionmut="lak-/-")
LD2_behavr <- create_behavr(LD2_data, LD2_meta)
LD2_behavr$distance[is.na(LD2_behavr$distance)]<-0
wt<-print(paste(25,"ctrl",sep=" "))
mut<-print(paste(25,"lak-/-",sep=" "))
plot_ggetho(LD2_behavr,titre="LD2 lak",phaseLD=1800, couleurs=c("white","black"),v)
#LD3
LD3_data <- transfo_data(LL_lak, lignes=436, colwt=192, nwt=15, colmut=247, nmut=15, minutes=10, conditionwt="ctrl", conditionmut="lak-/-")
LD3_meta <- create_meta(LD3_data, 15, 15, conditionwt="ctrl", conditionmut="lak-/-")
LD3_behavr <- create_behavr(LD3_data, LD3_meta)
LD3_behavr$distance[is.na(LD3_behavr$distance)]<-0
wt<-print(paste(15,"ctrl",sep=" "))
mut<-print(paste(15,"lak-/-",sep=" "))
plot_ggetho(LD3_behavr,titre="LD3 lak",phaseLD=1800, couleurs=c("white","black"),couleurmut="#CC79A7")
#LD4
LD4_data <- transfo_data(LL_lak, lignes=436, colwt=207, nwt=15, colmut=262, nmut=15, minutes=10, conditionwt="ctrl", conditionmut="lak-/-")
LD4_meta <- create_meta(LD4_data, 15, 15, conditionwt="ctrl", conditionmut="lak-/-")
LD4_behavr <- create_behavr(LD4_data, LD4_meta)
LD4_behavr$distance[is.na(LD4_behavr$distance)]<-0
wt<-print(paste(15,"ctrl",sep=" "))
mut<-print(paste(15,"lak-/-",sep=" "))
plot_ggetho(LD4_behavr,titre="LD4 lak",phaseLD=1800, couleurs=c("white","black"),couleurmut="#CC79A7")

# LDwt vs LLwt
LDvsLL_data <- transfo_data(LL_lak,lignes=436,colwt=167,nwt=55,colmut=3,nmut=81, minutes = 10,conditionwt="LD ctrl",conditionmut="LL ctrl")
LDvsLL_meta <- create_meta(LDvsLL_data, 55, 81, conditionwt="LD ctrl", conditionmut="LL ctrl")
LDvsLL_behavr <- create_behavr(LDvsLL_data, LDvsLL_meta)
LDvsLL_behavr$distance[is.na(LDvsLL_behavr$distance)]<-0
wt<-print(paste(55,"LD ctrl",sep=" "))
mut<-print(paste(81,"LL ctrl",sep=" "))
plot_ggetho(LDvsLL_behavr,titre="LDctrl vs LLctrl",phaseLD=1800, couleurs=c("white","black"),couleurmut="#56B4E9")
# LDmut vs LLwt

# LDwt vs LLmut
LDvsLL_data <- transfo_data(LL_lak,lignes=436,colwt=167,nwt=55,colmut=84,nmut=81, minutes = 10,conditionwt="LD ctrl",conditionmut="LL mut")
LDvsLL_meta <- create_meta(LDvsLL_data, 55, 81, conditionwt="LD ctrl", conditionmut="LL mut")
LDvsLL_behavr <- create_behavr(LDvsLL_data, LDvsLL_meta)
LDvsLL_behavr$distance[is.na(LDvsLL_behavr$distance)]<-0
wt<-print(paste(55,"LD ctrl",sep=" "))
mut<-print(paste(81,"LL mut",sep=" "))
plot_ggetho(LDvsLL_behavr,titre="LDctrl vs LLmut",phaseLD=1800, couleurs=c("white","black"),couleurmut="#CC79A7")
# LDmut vs LLmut
