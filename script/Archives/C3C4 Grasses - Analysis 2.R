################################
### C3C4 Grasses CALSTATE LA ###
########## Analysis 2 ##########
################################

### libraries ####
library(ggplot2)
library(agricolae)
library(plyr)
library(multcomp)
library(corrplot)
library(factoextra)
library(readxl)
### PV, data loading, variable calculation, means ####
#PVdata<-read.table("C:/Users/mnnad/OneDrive/Escritorio/R datasets/C3C4G PV R.txt",dec=".",header=T)
#on Marion's PC
PVdata<-read_xlsx(path = "C3C4 Grasses - PV.xlsx", sheet = "C3C4G PV R")

PVdata$Metabolism<-factor(PVdata$Metabolism)
PVdata$Species<-factor(PVdata$Species,levels=c("Avsa","Hovu","Trae","Pegl","Zema","Chga"))
# variable calculation
PVdata$T.LV<-PVdata$T.LA*PVdata$T.LT/10
PVdata$tlp.LV<-PVdata$tlp.LA*PVdata$tlp.LT/10
PVdata$D.LV<-PVdata$D.LA*PVdata$D.LT/10
PVdata$tlp.PLT<-100-PVdata$tlp.LT/PVdata$T.LT*100
PVdata$D.PLT<-100-PVdata$D.LT/PVdata$T.LT*100
PVdata$tlp.PLA<-100-PVdata$tlp.LA/PVdata$T.LA*100
PVdata$D.PLA<-100-PVdata$D.LA/PVdata$T.LA*100
PVdata$tlp.PLV<-100-PVdata$tlp.LV/PVdata$T.LV*100
PVdata$D.PLV<-100-PVdata$D.LV/PVdata$T.LV*100
PVdata$SWCm<-(PVdata$TW-PVdata$DW)/PVdata$DW
PVdata$LMA<-PVdata$DW/1000/(PVdata$T.LA/10000)
PVdata$LD<-PVdata$DW/1000/PVdata$T.LV
PVdata$LDMC<-PVdata$DW/PVdata$TW
PVdata$SWCpv<-PVdata$SWpv/PVdata$DW
PVdata$cft.mass<-PVdata$cft*PVdata$SWpv/PVdata$DW*1000
PVdata$cft.area<-PVdata$cft.mass/PVdata$LMA
PVdata$ctlp.mass<-PVdata$ctlp*PVdata$SWpv/PVdata$DW*1000
PVdata$ctlp.area<-PVdata$ctlp.mass/PVdata$LMA
PVdata$em<-PVdata$em.*(1-PVdata$af)
# means
m.PVdata<-ddply(PVdata,c("Species","Metabolism"),summarise,
            m.SWCpv=mean(SWCpv,na.rm=TRUE),n.SWCpv=sum(!is.na(SWCpv)),se.SWCpv=sd(SWCpv,na.rm=T)/sqrt(n.SWCpv),
            m.rwc.tlp=mean(rwc.tlp,na.rm=TRUE),n.rwc.tlp=sum(!is.na(rwc.tlp)),se.rwc.tlp=sd(rwc.tlp,na.rm=T)/sqrt(n.rwc.tlp),
            m.osp.tlp=mean(osp.tlp,na.rm=TRUE),n.osp.tlp=sum(!is.na(osp.tlp)),se.osp.tlp=sd(osp.tlp,na.rm=T)/sqrt(n.osp.tlp),
            m.osp.ft=mean(osp.ft,na.rm=TRUE),n.osp.ft=sum(!is.na(osp.ft)),se.osp.ft=sd(osp.ft,na.rm=T)/sqrt(n.osp.ft),
            m.em.=mean(em.,na.rm=TRUE),n.em.=sum(!is.na(em.)),se.em.=sd(em.,na.rm=T)/sqrt(n.em.),
            m.em=mean(em,na.rm=TRUE),n.em=sum(!is.na(em)),se.em=sd(em,na.rm=T)/sqrt(n.em),
            m.af=mean(af,na.rm=TRUE),n.af=sum(!is.na(af)),se.af=sd(af,na.rm=T)/sqrt(n.af),
            m.cft=mean(cft,na.rm=TRUE),n.cft=sum(!is.na(cft)),se.cft=sd(cft,na.rm=T)/sqrt(n.cft),
            m.ctlp=mean(ctlp,na.rm=TRUE),n.ctlp=sum(!is.na(ctlp)),se.ctlp=sd(ctlp,na.rm=T)/sqrt(n.ctlp),
            m.cft.area=mean(cft.area,na.rm=TRUE),n.cft.area=sum(!is.na(cft.area)),se.cft.area=sd(cft.area,na.rm=T)/sqrt(n.cft.area),
            m.cft.mass=mean(cft.mass,na.rm=TRUE),n.cft.mass=sum(!is.na(cft.mass)),se.cft.mass=sd(cft.mass,na.rm=T)/sqrt(n.cft.mass),
            m.ctlp.area=mean(ctlp.area,na.rm=TRUE),n.ctlp.area=sum(!is.na(ctlp.area)),se.ctlp.area=sd(ctlp.area,na.rm=T)/sqrt(n.ctlp.area),
            m.ctlp.mass=mean(ctlp.mass,na.rm=TRUE),n.ctlp.mass=sum(!is.na(ctlp.mass)),se.ctlp.mass=sd(ctlp.mass,na.rm=T)/sqrt(n.ctlp.mass),
            m.LMA=mean(LMA,na.rm=TRUE),n.LMA=sum(!is.na(LMA)),se.LMA=sd(LMA,na.rm=T)/sqrt(n.LMA),
            m.T.LT=mean(T.LT,na.rm=TRUE),n.T.LT=sum(!is.na(T.LT)),se.T.LT=sd(T.LT,na.rm=T)/sqrt(n.T.LT),
            m.LD=mean(LD,na.rm=TRUE),n.LD=sum(!is.na(LD)),se.LD=sd(LD,na.rm=T)/sqrt(n.LD),
            m.tlp.PLT=mean(tlp.PLT,na.rm=TRUE),n.tlp.PLT=sum(!is.na(tlp.PLT)),se.tlp.PLT=sd(tlp.PLT,na.rm=T)/sqrt(n.tlp.PLT),
            m.D.PLT=mean(D.PLT,na.rm=TRUE),n.D.PLT=sum(!is.na(D.PLT)),se.D.PLT=sd(D.PLT,na.rm=T)/sqrt(n.D.PLT),
            m.tlp.PLA=mean(tlp.PLA,na.rm=TRUE),n.tlp.PLA=sum(!is.na(tlp.PLA)),se.tlp.PLA=sd(tlp.PLA,na.rm=T)/sqrt(n.tlp.PLA),
            m.D.PLA=mean(D.PLA,na.rm=TRUE),n.D.PLA=sum(!is.na(D.PLA)),se.D.PLA=sd(D.PLA,na.rm=T)/sqrt(n.D.PLA),
            m.tlp.PLV=mean(tlp.PLV,na.rm=TRUE),n.tlp.PLV=sum(!is.na(tlp.PLV)),se.tlp.PLV=sd(tlp.PLV,na.rm=T)/sqrt(n.tlp.PLV),
            m.D.PLV=mean(D.PLV,na.rm=TRUE),n.D.PLV=sum(!is.na(D.PLV)),se.D.PLV=sd(D.PLV,na.rm=T)/sqrt(n.D.PLV))

### GEx, data loading, variable calculation, means ####
GExdata<-read.table("C:/Users/mnnad/OneDrive/Escritorio/R datasets/C3C4G GEx R.txt",dec=".",header=T) #Miquel

GExdata<-read_xlsx(path = "C3C4 Grasses - GEx.xlsx", sheet = "C3C4G GEx R")#on Marion's PC

GExdata$Dataset<-factor(GExdata$Dataset)
GExdata$Metabolism<-factor(GExdata$Metabolism)
GExdata$Species<-factor(GExdata$Species,levels=c("Avsa","Hovu","Trae","Pegl","Zema","Chga"))
GExdata$Individual<-factor(GExdata$Individual)
GExdata$Leaf<-factor(GExdata$Leaf)
GExdata$Time.min<-factor(GExdata$Time.min)
# variable calculation
GExdata$RWC<-(GExdata$FW-GExdata$DW)/(GExdata$RW-GExdata$DW)*100
# subsets
GExdata.0<-subset(GExdata,Time.min=="0")
GExdata.midd<-subset(GExdata,Dataset=="Midday")
GExdata.dehy<-subset(GExdata,Dataset=="Dehydration")
# means GEx midday
m.GExdata.midd<-ddply(GExdata.midd,c("Species","Metabolism"),summarise,
            m.An.area=mean(An.area,na.rm=TRUE),n.An.area=sum(!is.na(An.area)),se.An.area=sd(An.area,na.rm=T)/sqrt(n.An.area),
            m.gsw.area=mean(gsw.area,na.rm=TRUE),n.gsw.area=sum(!is.na(gsw.area)),se.gsw.area=sd(gsw.area,na.rm=T)/sqrt(n.gsw.area),
            m.Ci=mean(Ci,na.rm=TRUE),n.Ci=sum(!is.na(Ci)),se.Ci=sd(Ci,na.rm=T)/sqrt(n.Ci),
            m.WUEi=mean(WUEi,na.rm=TRUE),n.WUEi=sum(!is.na(WUEi)),se.WUEi=sd(WUEi,na.rm=T)/sqrt(n.WUEi),
            m.phiPSII=mean(phiPSII,na.rm=TRUE),n.phiPSII=sum(!is.na(phiPSII)),se.phiPSII=sd(phiPSII,na.rm=T)/sqrt(n.phiPSII))

### SRC, data loading, variable calculation, outliers, means ####
SRCdata<-read.table("C:/Users/mnnad/OneDrive/Escritorio/R datasets/C3C4G SRC R.txt",dec=".",header=T) #Miquel

SRCdata <- read_xlsx(path = "C3C4 Grasses - OV monitoring, SRC.xlsx", sheet = "C3C4G SRC") #marion

SRCdata$Metabolism<-factor(SRCdata$Metabolism)
SRCdata$Species<-factor(SRCdata$Species,levels=c("Avsa","Hovu","Trae","Pegl","Zema","Chga"))
# variable calculation
SRCdata$T.LT<-rowMeans(SRCdata[c("T.LT1","T.LT2","T.LT3","T.LT4","T.LT5","T.LT6")],na.rm=T)
SRCdata$F.LT<-rowMeans(SRCdata[c("F.LT1","F.LT2","F.LT3","F.LT4","F.LT5","F.LT6")],na.rm=T)
SRCdata$D.LT<-rowMeans(SRCdata[c("D.LT1","D.LT2","D.LT3","D.LT4","D.LT5","D.LT6")],na.rm=T)
SRCdata$T.LW1<-rowMeans(SRCdata[c("T.LW1","T.LW2","T.LW3")],na.rm=T)
SRCdata$F.LW1<-rowMeans(SRCdata[c("F.LW1","F.LW2","F.LW3")],na.rm=T)
SRCdata$T.LW2<-rowMeans(SRCdata[c("T.LW1","T.LW4","T.LW5")],na.rm=T)
SRCdata$F.LW2<-rowMeans(SRCdata[c("F.LW1","F.LW4","F.LW5")],na.rm=T)
SRCdata$D.LW<-rowMeans(SRCdata[c("D.LW1","D.LW2","D.LW3","D.LW4","D.LW5")],na.rm=T)
SRCdata$T.LV<-SRCdata$T.LA*SRCdata$T.LT/10
SRCdata$F.LV<-SRCdata$F.LA*SRCdata$F.LT/10
SRCdata$D.LV<-SRCdata$D.LA*SRCdata$D.LT/10
SRCdata$F.PLT<-100-SRCdata$F.LT/SRCdata$T.LT*100
SRCdata$D.PLT<-100-SRCdata$D.LT/SRCdata$T.LT*100
SRCdata$F.PLA<-100-SRCdata$F.LA/SRCdata$T.LA*100
SRCdata$D.PLA<-100-SRCdata$D.LA/SRCdata$T.LA*100
SRCdata$F.PLW1<-100-SRCdata$F.LW1/SRCdata$T.LW1*100
SRCdata$D.PLW1<-100-SRCdata$D.LW/SRCdata$T.LW1*100
SRCdata$F.PLLe<-100-SRCdata$F.LLe/SRCdata$T.LLe*100
SRCdata$D.PLLe<-100-SRCdata$D.LLe/SRCdata$T.LLe*100
SRCdata$F.PLV<-100-SRCdata$F.LV/SRCdata$T.LV*100
SRCdata$D.PLV<-100-SRCdata$D.LV/SRCdata$T.LV*100
SRCdata$T.SWC<-(SRCdata$TW-SRCdata$DW)/SRCdata$DW
SRCdata$F.SWC<-(SRCdata$FW-SRCdata$DW)/SRCdata$DW
SRCdata$R.SWC<-(SRCdata$RW-SRCdata$DW)/SRCdata$DW
SRCdata$F.RWC<-(SRCdata$FW-SRCdata$DW)/(SRCdata$TW-SRCdata$DW)*100
SRCdata$F.1_LWi.LWhy<-1-SRCdata$FW/SRCdata$TW
SRCdata$LMA<-SRCdata$DW/1000/(SRCdata$T.LA/10000)
SRCdata$LD<-SRCdata$DW/1000/SRCdata$T.LV
SRCdata$LDMC<-SRCdata$DW/SRCdata$TW
SRCdata$PLRC<-100-SRCdata$R.SWC/SRCdata$T.SWC*100
## SRC, outliers
# PLRC of -20% (Hovu)
SRCdata$PLRC[!SRCdata$PLRC>-20]<-NA
# Avsa
#SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==17,5:84]<-NA
#SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==20,5:84]<-NA
#SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==21,5:84]<-NA
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==22,5:84]<-NA # FvFm
#SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==24,5:84]<-NA
#SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==25,5:84]<-NA
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==26,5:84]<-NA # FvFm
#SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==29,5:84]<-NA
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==30,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==31,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==32,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==43,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==44,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==45,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==46,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==47,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==49,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==51,5:84]<-NA # too small leaf (!)
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==52,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==53,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==54,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==55,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==56,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==57,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==59,5:84]<-NA # PLA
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==17,84]<-NA # only PLRC
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==20,84]<-NA # only PLRC
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==21,84]<-NA # only PLRC
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==24,84]<-NA # only PLRC
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==25,84]<-NA # only PLRC
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==24,74]<-NA # only D.PLLe - to check photos!
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==27,74]<-NA # only D.PLLe - to check photos!
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==28,74]<-NA # only D.PLLe - to check photos!
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==65,74]<-NA # only D.PLLe - to check photos!
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==66,74]<-NA # only D.PLLe - to check photos!
SRCdata[SRCdata$Species=="Avsa" & SRCdata$Leaf==48,70]<-NA # only D.PLA - to check photos!
# Chga
SRCdata[SRCdata$Species=="Chga" & SRCdata$Leaf==1,5:84]<-NA # PLT - to check
SRCdata[SRCdata$Species=="Chga" & SRCdata$Leaf==6,5:84]<-NA # PLA - to check photos!
SRCdata[SRCdata$Species=="Chga" & SRCdata$Leaf==9,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Chga" & SRCdata$Leaf==21,5:84]<-NA # PLT - to check
SRCdata[SRCdata$Species=="Chga" & SRCdata$Leaf==22,5:84]<-NA # LMA vs LDMC
SRCdata[SRCdata$Species=="Chga" & SRCdata$Leaf==24,74]<-NA # only D.PLLe - to check photos!
SRCdata[SRCdata$Species=="Chga" & SRCdata$Leaf==25,5:84]<-NA # LMA vs LDMC
SRCdata[SRCdata$Species=="Chga" & SRCdata$Leaf==30,5:84]<-NA # PLT - to check
SRCdata[SRCdata$Species=="Chga" & SRCdata$Leaf==31,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Chga" & SRCdata$Leaf==32,5:84]<-NA # PLRC
SRCdata[SRCdata$Species=="Chga" & SRCdata$Leaf==48,5:84]<-NA # FvFm
# Hovu
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==3,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==4,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==5,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==6,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==7,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==8,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==15,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==17,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==19,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==49,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==50,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==56,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==55,84]<-NA # only PLRC
SRCdata[SRCdata$Species=="Hovu" & SRCdata$Leaf==14,c(69,75)]<-NA # only PLA, PLV - to check photos!
# Pegl 
SRCdata[SRCdata$Species=="Pegl" & SRCdata$Leaf==35,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Pegl" & SRCdata$Leaf==38,74]<-NA # only D.PLLe - to check photos!
SRCdata[SRCdata$Species=="Pegl" & SRCdata$Leaf==42,74]<-NA # only D.PLLe - to check photos!
# Trae
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==8,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==9,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==10,5:84]<-NA # FvFm
#SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==18,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==20,5:84]<-NA # FvFm
#SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==23,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==28,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==29,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==30,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==31,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==32,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==45,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==46,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==47,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==49,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==50,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==52,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==53,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Trae" & SRCdata$Leaf==42,44]<-NA # only R FvFm - to check!
# Zema
SRCdata[SRCdata$Species=="Zema" & SRCdata$Leaf==8,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Zema" & SRCdata$Leaf==10,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Zema" & SRCdata$Leaf==15,5:84]<-NA # LDMC vs LD
SRCdata[SRCdata$Species=="Zema" & SRCdata$Leaf==20,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Zema" & SRCdata$Leaf==23,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Zema" & SRCdata$Leaf==34,5:84]<-NA # FvFm
SRCdata[SRCdata$Species=="Zema" & SRCdata$Leaf==24,84]<-NA # only PLRC
SRCdata[SRCdata$Species=="Zema" & SRCdata$Leaf==6,c(70,76)]<-NA # only DPLA, DPLV - to check photos!
## means 
m.SRCdata<-ddply(SRCdata,c("Species","Metabolism"),summarise,
                m.LMA=mean(LMA,na.rm=TRUE),n.LMA=sum(!is.na(LMA)),se.LMA=sd(LMA,na.rm=T)/sqrt(n.LMA),
                m.T.LT=mean(T.LT,na.rm=TRUE),n.T.LT=sum(!is.na(T.LT)),se.T.LT=sd(T.LT,na.rm=T)/sqrt(n.T.LT),
                m.LD=mean(LD,na.rm=TRUE),n.LD=sum(!is.na(LD)),se.LD=sd(LD,na.rm=T)/sqrt(n.LD),
                m.LDMC=mean(LDMC,na.rm=TRUE),n.LDMC=sum(!is.na(LDMC)),se.LDMC=sd(LDMC,na.rm=T)/sqrt(n.LDMC),
                m.D.PLT=mean(D.PLT,na.rm=TRUE),n.D.PLT=sum(!is.na(D.PLT)),se.D.PLT=sd(D.PLT,na.rm=T)/sqrt(n.D.PLT),
                m.D.PLA=mean(D.PLA,na.rm=TRUE),n.D.PLA=sum(!is.na(D.PLA)),se.D.PLA=sd(D.PLA,na.rm=T)/sqrt(n.D.PLA),
                m.D.PLV=mean(D.PLV,na.rm=TRUE),n.D.PLV=sum(!is.na(D.PLV)),se.D.PLV=sd(D.PLV,na.rm=T)/sqrt(n.D.PLV),
                m.D.PLLe=mean(D.PLLe,na.rm=TRUE),n.D.PLLe=sum(!is.na(D.PLLe)),se.D.PLLe=sd(D.PLLe,na.rm=T)/sqrt(n.D.PLLe),
                m.D.PLW1=mean(D.PLW1,na.rm=TRUE),n.D.PLW1=sum(!is.na(D.PLW1)),se.D.PLW1=sd(D.PLW1,na.rm=T)/sqrt(n.D.PLW1))

### X50, data loading, subsets ####
X50data<-read.table("C:/Users/mnnad/OneDrive/Escritorio/R datasets/C3C4G X50 R.txt",dec=".",header=T)
X50data$WaterStatus<-factor(X50data$WaterStatus)
X50data$Species<-factor(X50data$Species,levels=c("Avsa","Hovu","Trae","Pegl","Zema","Chga"))
X50data$Metabolism<-factor(X50data$Metabolism)
# subsets
X50data.RWC<-subset(X50data,WaterStatus=="RWC")
## PV, stat analysis ####
lm1<-lm(SWCpv~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(rwc.tlp~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(osp.tlp~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(osp.ft~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(em.~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(em~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(af~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(cft~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(ctlp~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(cft.area~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(ctlp.area~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(cft.mass~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(ctlp.mass~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(LMA~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(T.LT~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(LD~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(tlp.PLT~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(D.PLT~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(tlp.PLA~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(D.PLA~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(tlp.PLV~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(D.PLV~Species,data=PVdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups

## GEx (midday), stat analysis ####
lm1<-lm(An.area~Species,data=GExdata.midd)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(gsw.area~Species,data=GExdata.midd)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(Ci~Species,data=GExdata.midd)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(WUEi~Species,data=GExdata.midd)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
lm1<-lm(phiPSII~Species,data=GExdata.midd)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
# Midday water potential check ####
g1<-ggplot(data=GExdata.0[GExdata.0$Dataset=="Midday",],aes(x=Species,y=F.wp.leaf,fill=Metabolism))+
  geom_boxplot(lwd=0.2,fatten=5)+
  geom_point(data=GExdata.0[GExdata.0$Dataset=="Dehydration",],aes(x=Species,y=F.wp.leaf),size=4,shape=21,fill="yellow")+
  annotate("point",y=-1.02,x=1,fill="purple",size=5,alpha=0.5,shape=21)+
  annotate("point",y=-1.00,x=2,fill="purple",size=5,alpha=0.5,shape=21)+
  annotate("point",y=-1.00,x=3,fill="purple",size=5,alpha=0.5,shape=21)+
  annotate("point",y=-0.77,x=4,fill="purple",size=5,alpha=0.5,shape=21)+
  annotate("point",y=-0.74,x=5,fill="purple",size=5,alpha=0.5,shape=21)+
  annotate("point",y=-0.84,x=6,fill="purple",size=5,alpha=0.5,shape=21)+
  #scale_fill_manual(values=c("gray99","palegreen","green4"))+
  ylim(-2,0)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13,face="italic"),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste(Psi[leaf]," (MPa)")))
## GEx (dehydration) #### 
# quick graph for calculating RWC50
g1<-ggplot()+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Chga",],aes(x=RWC,y=An.area),size=4,shape=21,fill="red")+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Chga",],aes(x=RWC,y=An.area),color="red",alpha=0.5,fill="tomato1")+
  xlim(0,100)+
  annotate("text",x=10,y=0,label=expression(paste(italic("Chga"))))
# graphs per species ####
# Avsa
g1<-ggplot()+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_hline(yintercept=50,color="darkgrey",lty=2)+
  geom_hline(yintercept=100,color="darkgrey",lty=2)+
  geom_vline(xintercept=73.02,color="red",lty=2,lwd=1)+
  geom_vline(xintercept=77.76,color="blue",lty=2,lwd=1)+
  geom_vline(xintercept=46.12,color="green",lty=2,lwd=1)+
  geom_vline(xintercept=89.5,color="black",lty=2,lwd=1)+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Avsa",],aes(x=RWC,y=phiPSII/0.364*100),size=3,shape=21,fill="green",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Avsa",],aes(x=RWC,y=phiPSII/0.364*100),color="green",alpha=0.5,fill="palegreen")+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Avsa",],aes(x=RWC,y=An.area/25.20*100),size=3,shape=21,fill="red",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Avsa",],aes(x=RWC,y=An.area/25.20*100),color="red",alpha=0.5,fill="tomato1")+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Avsa",],aes(x=RWC,y=gsw.area/0.417*100),size=3,shape=21,fill="blue",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Avsa",],aes(x=RWC,y=gsw.area/0.417*100),color="blue",alpha=0.5,fill="skyblue2")+
  ylim(-25,150)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_text(size=14),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("relative A"[n],", g"[sw],", ",Phi["PSII"]," (%)")))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=10,y=90,label=expression(paste(italic("Avsa"))))
# Hovu
g1<-ggplot()+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_hline(yintercept=50,color="darkgrey",lty=2)+
  geom_hline(yintercept=100,color="darkgrey",lty=2)+
  geom_vline(xintercept=75.78,color="red",lty=2,lwd=1)+
  geom_vline(xintercept=83.95,color="blue",lty=2,lwd=1)+
  geom_vline(xintercept=66.30,color="green",lty=2,lwd=1)+
  geom_vline(xintercept=90.7,color="black",lty=2,lwd=1)+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Hovu",],aes(x=RWC,y=phiPSII/0.385*100),size=3,shape=21,fill="green",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Hovu",],aes(x=RWC,y=phiPSII/0.385*100),color="green",alpha=0.5,fill="palegreen")+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Hovu",],aes(x=RWC,y=An.area/28.44*100),size=3,shape=21,fill="red",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Hovu",],aes(x=RWC,y=An.area/28.44*100),color="red",alpha=0.5,fill="tomato1")+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Hovu",],aes(x=RWC,y=gsw.area/0.501*100),size=3,shape=21,fill="blue",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Hovu",],aes(x=RWC,y=gsw.area/0.501*100),color="blue",alpha=0.5,fill="skyblue2")+
  ylim(-25,150)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_text(size=14),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("relative A"[n],", g"[sw],", ",Phi["PSII"]," (%)")))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=10,y=90,label=expression(paste(italic("Hovu"))))
# Trae
g1<-ggplot()+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_hline(yintercept=50,color="darkgrey",lty=2)+
  geom_hline(yintercept=100,color="darkgrey",lty=2)+
  geom_vline(xintercept=73.37,color="red",lty=2,lwd=1)+
  geom_vline(xintercept=81.31,color="blue",lty=2,lwd=1)+
  geom_vline(xintercept=61.15,color="green",lty=2,lwd=1)+
  geom_vline(xintercept=90.0,color="black",lty=2,lwd=1)+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Trae",],aes(x=RWC,y=phiPSII/0.358*100),size=3,shape=21,fill="green",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Trae",],aes(x=RWC,y=phiPSII/0.358*100),color="green",alpha=0.5,fill="palegreen")+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Trae",],aes(x=RWC,y=An.area/32.91*100),size=3,shape=21,fill="red",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Trae",],aes(x=RWC,y=An.area/32.91*100),color="red",alpha=0.5,fill="tomato1")+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Trae",],aes(x=RWC,y=gsw.area/0.699*100),size=3,shape=21,fill="blue",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Trae",],aes(x=RWC,y=gsw.area/0.699*100),color="blue",alpha=0.5,fill="skyblue2")+
  ylim(-25,150)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_text(size=14),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("relative A"[n],", g"[sw],", ",Phi["PSII"]," (%)")))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=10,y=90,label=expression(paste(italic("Trae"))))
# Pegl
g1<-ggplot()+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_hline(yintercept=50,color="darkgrey",lty=2)+
  geom_hline(yintercept=100,color="darkgrey",lty=2)+
  geom_vline(xintercept=86.94,color="red",lty=2,lwd=1)+
  geom_vline(xintercept=86.16,color="blue",lty=2,lwd=1)+
  geom_vline(xintercept=83.79,color="green",lty=2,lwd=1)+
  geom_vline(xintercept=91.2,color="black",lty=2,lwd=1)+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Pegl",],aes(x=RWC,y=phiPSII/0.350*100),size=3,shape=21,fill="green",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Pegl",],aes(x=RWC,y=phiPSII/0.350*100),color="green",alpha=0.5,fill="palegreen")+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Pegl",],aes(x=RWC,y=An.area/40.28*100),size=3,shape=21,fill="red",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Pegl",],aes(x=RWC,y=An.area/40.28*100),color="red",alpha=0.5,fill="tomato1")+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Pegl",],aes(x=RWC,y=gsw.area/0.230*100),size=3,shape=21,fill="blue",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Pegl",],aes(x=RWC,y=gsw.area/0.230*100),color="blue",alpha=0.5,fill="skyblue2")+
  ylim(-25,150)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_text(size=14),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("relative A"[n],", g"[sw],", ",Phi["PSII"]," (%)")))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=10,y=90,label=expression(paste(italic("Pegl"))))
# Chga
g1<-ggplot()+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_hline(yintercept=50,color="darkgrey",lty=2)+
  geom_hline(yintercept=100,color="darkgrey",lty=2)+
  geom_vline(xintercept=80.90,color="red",lty=2,lwd=1)+
  geom_vline(xintercept=82.80,color="blue",lty=2,lwd=1)+
  geom_vline(xintercept=77.47,color="green",lty=2,lwd=1)+
  geom_vline(xintercept=93.1,color="black",lty=2,lwd=1)+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Chga",],aes(x=RWC,y=phiPSII/0.252*100),size=3,shape=21,fill="green",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Chga",],aes(x=RWC,y=phiPSII/0.252*100),color="green",alpha=0.5,fill="palegreen")+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Chga",],aes(x=RWC,y=An.area/26.75*100),size=3,shape=21,fill="red",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Chga",],aes(x=RWC,y=An.area/26.75*100),color="red",alpha=0.5,fill="tomato1")+
  geom_point(data=GExdata.dehy[GExdata.dehy$Species=="Chga",],aes(x=RWC,y=gsw.area/0.154*100),size=3,shape=21,fill="blue",alpha=0.5)+
  geom_smooth(data=GExdata.dehy[GExdata.dehy$Species=="Chga",],aes(x=RWC,y=gsw.area/0.154*100),color="blue",alpha=0.5,fill="skyblue2")+
  ylim(-25,150)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_text(size=14),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("relative A"[n],", g"[sw],", ",Phi["PSII"]," (%)")))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=10,y=90,label=expression(paste(italic("Chga"))))

## SRC, boxplots, stat analysis ####
# LMA
g1<-ggplot(SRCdata,aes(x=Species,y=LMA,fill=Metabolism))+
  geom_boxplot(lwd=0.2,fatten=5)+
  #geom_point(data=GExdata.0[GExdata.0$Dataset=="Dehydration",],aes(x=Species,y=LMA),size=4,shape=21,fill="yellow")+
  ylim(0,50)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13,face="italic"),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("LMA (g m"^-2,")")))+
  annotate("text",x=1,y=7,label=expression(paste(italic(bar("x"))," = 31.1")),size=4)+
  annotate("text",x=1,y=4,label=expression(paste(italic("n")," = 49")),size=4)+
  annotate("text",x=2,y=7,label=expression(paste("21.2")),size=4)+
  annotate("text",x=2,y=4,label=expression(paste("48")),size=4)+
  annotate("text",x=3,y=7,label=expression(paste("27.4")),size=4)+
  annotate("text",x=3,y=4,label=expression(paste("38")),size=4)+
  annotate("text",x=4,y=7,label=expression(paste("24.7")),size=4)+
  annotate("text",x=4,y=4,label=expression(paste("41")),size=4)+
  annotate("text",x=5,y=7,label=expression(paste("20.6")),size=4)+
  annotate("text",x=5,y=4,label=expression(paste("42")),size=4)+
  annotate("text",x=6,y=7,label=expression(paste("27.5")),size=4)+
  annotate("text",x=6,y=4,label=expression(paste("42")),size=4)+
  annotate("text",x=1,y=max(SRCdata[SRCdata$Species=="Avsa","LMA"],na.rm=T)+4,label=expression(paste(bold("a"))),size=5)+
  annotate("text",x=2,y=max(SRCdata[SRCdata$Species=="Hovu","LMA"],na.rm=T)+4,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=3,y=max(SRCdata[SRCdata$Species=="Trae","LMA"],na.rm=T)+4,label=expression(paste(bold("b"))),size=5)+
  annotate("text",x=4,y=max(SRCdata[SRCdata$Species=="Pegl","LMA"],na.rm=T)+4,label=expression(paste(bold("b"))),size=5)+
  annotate("text",x=5,y=max(SRCdata[SRCdata$Species=="Zema","LMA"],na.rm=T)+4,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=6,y=max(SRCdata[SRCdata$Species=="Chga","LMA"],na.rm=T)+4,label=expression(paste(bold("b"))),size=5)
# stat
lm1<-lm(LMA~Species,data=SRCdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
# LT
g1<-ggplot(SRCdata,aes(x=Species,y=T.LT,fill=Metabolism))+
  geom_boxplot(lwd=0.2,fatten=5)+
  #geom_point(data=GExdata.0[GExdata.0$Dataset=="Dehydration",],aes(x=Species,y=LMA),size=4,shape=21,fill="yellow")+
  ylim(0,0.4)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13,face="italic"),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("LT (mm)")))+
  annotate("text",x=1,y=0.056,label=expression(paste(italic(bar("x"))," = 0.22")),size=4)+
  annotate("text",x=1,y=0.032,label=expression(paste(italic("n")," = 49")),size=4)+
  annotate("text",x=2,y=0.056,label=expression(paste("0.19")),size=4)+
  annotate("text",x=2,y=0.032,label=expression(paste("48")),size=4)+
  annotate("text",x=3,y=0.056,label=expression(paste("0.17")),size=4)+
  annotate("text",x=3,y=0.032,label=expression(paste("38")),size=4)+
  annotate("text",x=4,y=0.056,label=expression(paste("0.13")),size=4)+
  annotate("text",x=4,y=0.032,label=expression(paste("41")),size=4)+
  annotate("text",x=5,y=0.056,label=expression(paste("0.16")),size=4)+
  annotate("text",x=5,y=0.032,label=expression(paste("42")),size=4)+
  annotate("text",x=6,y=0.056,label=expression(paste("0.12")),size=4)+
  annotate("text",x=6,y=0.032,label=expression(paste("42")),size=4)+
  annotate("text",x=1,y=max(SRCdata[SRCdata$Species=="Avsa","T.LT"],na.rm=T)+0.032,label=expression(paste(bold("a"))),size=5)+
  annotate("text",x=2,y=max(SRCdata[SRCdata$Species=="Hovu","T.LT"],na.rm=T)+0.032,label=expression(paste(bold("b"))),size=5)+
  annotate("text",x=3,y=max(SRCdata[SRCdata$Species=="Trae","T.LT"],na.rm=T)+0.032,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=4,y=max(SRCdata[SRCdata$Species=="Pegl","T.LT"],na.rm=T)+0.032,label=expression(paste(bold("d"))),size=5)+
  annotate("text",x=5,y=max(SRCdata[SRCdata$Species=="Zema","T.LT"],na.rm=T)+0.032,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=6,y=max(SRCdata[SRCdata$Species=="Chga","T.LT"],na.rm=T)+0.032,label=expression(paste(bold("e"))),size=5)
# stat
lm1<-lm(T.LT~Species,data=SRCdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
# LD
g1<-ggplot(SRCdata,aes(x=Species,y=LD,fill=Metabolism))+
  geom_boxplot(lwd=0.2,fatten=5)+
  #geom_point(data=GExdata.0[GExdata.0$Dataset=="Dehydration",],aes(x=Species,y=LMA),size=4,shape=21,fill="yellow")+
  ylim(0,0.4)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13,face="italic"),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("LD (g cm"^-3,")")))+
  annotate("text",x=1,y=0.040,label=expression(paste(italic(bar("x"))," = 0.14")),size=4)+
  annotate("text",x=1,y=0.016,label=expression(paste(italic("n")," = 49")),size=4)+
  annotate("text",x=2,y=0.040,label=expression(paste("0.11")),size=4)+
  annotate("text",x=2,y=0.016,label=expression(paste("48")),size=4)+
  annotate("text",x=3,y=0.040,label=expression(paste("0.17")),size=4)+
  annotate("text",x=3,y=0.016,label=expression(paste("38")),size=4)+
  annotate("text",x=4,y=0.040,label=expression(paste("0.19")),size=4)+
  annotate("text",x=4,y=0.016,label=expression(paste("41")),size=4)+
  annotate("text",x=5,y=0.040,label=expression(paste("0.13")),size=4)+
  annotate("text",x=5,y=0.016,label=expression(paste("42")),size=4)+
  annotate("text",x=6,y=0.040,label=expression(paste("0.24")),size=4)+
  annotate("text",x=6,y=0.016,label=expression(paste("42")),size=4)+
  annotate("text",x=1,y=max(SRCdata[SRCdata$Species=="Avsa","LD"],na.rm=T)+0.032,label=expression(paste(bold("d"))),size=5)+
  annotate("text",x=2,y=max(SRCdata[SRCdata$Species=="Hovu","LD"],na.rm=T)+0.032,label=expression(paste(bold("e"))),size=5)+
  annotate("text",x=3,y=max(SRCdata[SRCdata$Species=="Trae","LD"],na.rm=T)+0.032,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=4,y=max(SRCdata[SRCdata$Species=="Pegl","LD"],na.rm=T)+0.032,label=expression(paste(bold("b"))),size=5)+
  annotate("text",x=5,y=max(SRCdata[SRCdata$Species=="Zema","LD"],na.rm=T)+0.032,label=expression(paste(bold("d"))),size=5)+
  annotate("text",x=6,y=max(SRCdata[SRCdata$Species=="Chga","LD"],na.rm=T)+0.032,label=expression(paste(bold("a"))),size=5)
# stat
lm1<-lm(LD~Species,data=SRCdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
# LDMC
g1<-ggplot(SRCdata,aes(x=Species,y=LDMC,fill=Metabolism))+
  geom_boxplot(lwd=0.2,fatten=5)+
  #geom_point(data=GExdata.0[GExdata.0$Dataset=="Dehydration",],aes(x=Species,y=LMA),size=4,shape=21,fill="yellow")+
  ylim(0,0.3)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13,face="italic"),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("LDMC (g g"^-1,")")))+
  annotate("text",x=1,y=0.024,label=expression(paste(italic(bar("x"))," = 0.13")),size=4)+
  annotate("text",x=1,y=0,label=expression(paste(italic("n")," = 44")),size=4)+
  annotate("text",x=2,y=0.024,label=expression(paste("0.11")),size=4)+
  annotate("text",x=2,y=0,label=expression(paste("47")),size=4)+
  annotate("text",x=3,y=0.024,label=expression(paste("0.16")),size=4)+
  annotate("text",x=3,y=0,label=expression(paste("38")),size=4)+
  annotate("text",x=4,y=0.024,label=expression(paste("0.11")),size=4)+
  annotate("text",x=4,y=0,label=expression(paste("41")),size=4)+
  annotate("text",x=5,y=0.024,label=expression(paste("0.10")),size=4)+
  annotate("text",x=5,y=0,label=expression(paste("42")),size=4)+
  annotate("text",x=6,y=0.024,label=expression(paste("0.16")),size=4)+
  annotate("text",x=6,y=0,label=expression(paste("42")),size=4)+
  annotate("text",x=1,y=max(SRCdata[SRCdata$Species=="Avsa","LDMC"],na.rm=T)+0.024,label=expression(paste(bold("b"))),size=5)+
  annotate("text",x=2,y=max(SRCdata[SRCdata$Species=="Hovu","LDMC"],na.rm=T)+0.024,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=3,y=max(SRCdata[SRCdata$Species=="Trae","LDMC"],na.rm=T)+0.024,label=expression(paste(bold("a"))),size=5)+
  annotate("text",x=4,y=max(SRCdata[SRCdata$Species=="Pegl","LDMC"],na.rm=T)+0.024,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=5,y=max(SRCdata[SRCdata$Species=="Zema","LDMC"],na.rm=T)+0.024,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=6,y=max(SRCdata[SRCdata$Species=="Chga","LDMC"],na.rm=T)+0.024,label=expression(paste(bold("a"))),size=5)
# stat
lm1<-lm(LDMC~Species,data=SRCdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
# PLTdry
g1<-ggplot(SRCdata,aes(x=Species,y=D.PLT,fill=Metabolism))+
  geom_boxplot(lwd=0.2,fatten=5)+
  #geom_point(data=GExdata.0[GExdata.0$Dataset=="Dehydration",],aes(x=Species,y=LMA),size=4,shape=21,fill="yellow")+
  ylim(0,100)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13,face="italic"),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("PLT"["dry"]," (%)")))+
  annotate("text",x=1,y=10,label=expression(paste(italic(bar("x"))," = 46.1")),size=4)+
  annotate("text",x=1,y=4,label=expression(paste(italic("n")," = 49")),size=4)+
  annotate("text",x=2,y=10,label=expression(paste("69.8")),size=4)+
  annotate("text",x=2,y=4,label=expression(paste("48")),size=4)+
  annotate("text",x=3,y=10,label=expression(paste("45.8")),size=4)+
  annotate("text",x=3,y=4,label=expression(paste("38")),size=4)+
  annotate("text",x=4,y=10,label=expression(paste("56.5")),size=4)+
  annotate("text",x=4,y=4,label=expression(paste("41")),size=4)+
  annotate("text",x=5,y=10,label=expression(paste("69.2")),size=4)+
  annotate("text",x=5,y=4,label=expression(paste("42")),size=4)+
  annotate("text",x=6,y=10,label=expression(paste("46.4")),size=4)+
  annotate("text",x=6,y=4,label=expression(paste("42")),size=4)+
  annotate("text",x=1,y=max(SRCdata[SRCdata$Species=="Avsa","D.PLT"],na.rm=T)+8,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=2,y=max(SRCdata[SRCdata$Species=="Hovu","D.PLT"],na.rm=T)+8,label=expression(paste(bold("a"))),size=5)+
  annotate("text",x=3,y=max(SRCdata[SRCdata$Species=="Trae","D.PLT"],na.rm=T)+8,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=4,y=max(SRCdata[SRCdata$Species=="Pegl","D.PLT"],na.rm=T)+8,label=expression(paste(bold("b"))),size=5)+
  annotate("text",x=5,y=max(SRCdata[SRCdata$Species=="Zema","D.PLT"],na.rm=T)+8,label=expression(paste(bold("a"))),size=5)+
  annotate("text",x=6,y=max(SRCdata[SRCdata$Species=="Chga","D.PLT"],na.rm=T)+8,label=expression(paste(bold("c"))),size=5)
# stat
lm1<-lm(D.PLT~Species,data=SRCdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
# PLAdry
g1<-ggplot(SRCdata,aes(x=Species,y=D.PLA,fill=Metabolism))+
  geom_boxplot(lwd=0.2,fatten=5)+
  #geom_point(data=GExdata.0[GExdata.0$Dataset=="Dehydration",],aes(x=Species,y=LMA),size=4,shape=21,fill="yellow")+
  ylim(10,90)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13,face="italic"),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("PLA"["dry"]," (%)")))+
  annotate("text",x=1,y=21.2,label=expression(paste(italic(bar("x"))," = 49.4")),size=4)+
  annotate("text",x=1,y=16.4,label=expression(paste(italic("n")," = 48")),size=4)+
  annotate("text",x=2,y=21.2,label=expression(paste("43.7")),size=4)+
  annotate("text",x=2,y=16.4,label=expression(paste("46")),size=4)+
  annotate("text",x=3,y=21.2,label=expression(paste("45.8")),size=4)+
  annotate("text",x=3,y=16.4,label=expression(paste("38")),size=4)+
  annotate("text",x=4,y=21.2,label=expression(paste("60.2")),size=4)+
  annotate("text",x=4,y=16.4,label=expression(paste("41")),size=4)+
  annotate("text",x=5,y=21.2,label=expression(paste("53.0")),size=4)+
  annotate("text",x=5,y=16.4,label=expression(paste("42")),size=4)+
  annotate("text",x=6,y=21.2,label=expression(paste("63.1")),size=4)+
  annotate("text",x=6,y=16.4,label=expression(paste("42")),size=4)+
  annotate("text",x=1,y=max(SRCdata[SRCdata$Species=="Avsa","D.PLA"],na.rm=T)+6.4,label=expression(paste(bold("bc"))),size=5)+
  annotate("text",x=2,y=max(SRCdata[SRCdata$Species=="Hovu","D.PLA"],na.rm=T)+6.4,label=expression(paste(bold("d"))),size=5)+
  annotate("text",x=3,y=max(SRCdata[SRCdata$Species=="Trae","D.PLA"],na.rm=T)+6.4,label=expression(paste(bold("cd"))),size=5)+
  annotate("text",x=4,y=max(SRCdata[SRCdata$Species=="Pegl","D.PLA"],na.rm=T)+6.4,label=expression(paste(bold("a"))),size=5)+
  annotate("text",x=5,y=max(SRCdata[SRCdata$Species=="Zema","D.PLA"],na.rm=T)+6.4,label=expression(paste(bold("b"))),size=5)+
  annotate("text",x=6,y=max(SRCdata[SRCdata$Species=="Chga","D.PLA"],na.rm=T)+6.4,label=expression(paste(bold("a"))),size=5)
# stat
lm1<-lm(D.PLA~Species,data=SRCdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
# PLLedry
g1<-ggplot(SRCdata,aes(x=Species,y=D.PLLe,fill=Metabolism))+
  geom_boxplot(lwd=0.2,fatten=5)+
  #geom_point(data=GExdata.0[GExdata.0$Dataset=="Dehydration",],aes(x=Species,y=LMA),size=4,shape=21,fill="yellow")+
  ylim(-15,20)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13,face="italic"),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("PLLe"["dry"]," (%)")))+
  annotate("text",x=1,y=-6,label=expression(paste(italic(bar("x"))," = 1.9")),size=4)+
  annotate("text",x=1,y=-8,label=expression(paste(italic("n")," = 41")),size=4)+
  annotate("text",x=2,y=-6,label=expression(paste("2.0")),size=4)+
  annotate("text",x=2,y=-8,label=expression(paste("46")),size=4)+
  annotate("text",x=3,y=-6,label=expression(paste("1.4")),size=4)+
  annotate("text",x=3,y=-8,label=expression(paste("38")),size=4)+
  annotate("text",x=4,y=-6,label=expression(paste("2.9")),size=4)+
  annotate("text",x=4,y=-8,label=expression(paste("39")),size=4)+
  annotate("text",x=5,y=-6,label=expression(paste("3.2")),size=4)+
  annotate("text",x=5,y=-8,label=expression(paste("42")),size=4)+
  annotate("text",x=6,y=-6,label=expression(paste("3.4")),size=4)+
  annotate("text",x=6,y=-8,label=expression(paste("41")),size=4)+
  annotate("text",x=1,y=max(SRCdata[SRCdata$Species=="Avsa","D.PLLe"],na.rm=T)+3,label=expression(paste(bold("bc"))),size=5)+
  annotate("text",x=2,y=max(SRCdata[SRCdata$Species=="Hovu","D.PLLe"],na.rm=T)+3,label=expression(paste(bold("bc"))),size=5)+
  annotate("text",x=3,y=max(SRCdata[SRCdata$Species=="Trae","D.PLLe"],na.rm=T)+3,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=4,y=max(SRCdata[SRCdata$Species=="Pegl","D.PLLe"],na.rm=T)+3,label=expression(paste(bold("ab"))),size=5)+
  annotate("text",x=5,y=max(SRCdata[SRCdata$Species=="Zema","D.PLLe"],na.rm=T)+3,label=expression(paste(bold("a"))),size=5)+
  annotate("text",x=6,y=max(SRCdata[SRCdata$Species=="Chga","D.PLLe"],na.rm=T)+3,label=expression(paste(bold("a"))),size=5)
# stat
lm1<-lm(D.PLLe~Species,data=SRCdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
# PLWidry (PLW1dry)
g1<-ggplot(SRCdata,aes(x=Species,y=D.PLW1,fill=Metabolism))+
  geom_boxplot(lwd=0.2,fatten=5)+
  #geom_point(data=GExdata.0[GExdata.0$Dataset=="Dehydration",],aes(x=Species,y=LMA),size=4,shape=21,fill="yellow")+
  ylim(0,100)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13,face="italic"),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("PLWi"["dry"]," (%)")))+
  annotate("text",x=1,y=14,label=expression(paste(italic(bar("x"))," = 49.0")),size=4)+
  annotate("text",x=1,y=8,label=expression(paste(italic("n")," = 49")),size=4)+
  annotate("text",x=2,y=14,label=expression(paste("44.9")),size=4)+
  annotate("text",x=2,y=8,label=expression(paste("46")),size=4)+
  annotate("text",x=3,y=14,label=expression(paste("46.8")),size=4)+
  annotate("text",x=3,y=8,label=expression(paste("38")),size=4)+
  annotate("text",x=4,y=14,label=expression(paste("57.9")),size=4)+
  annotate("text",x=4,y=8,label=expression(paste("41")),size=4)+
  annotate("text",x=5,y=14,label=expression(paste("56.3")),size=4)+
  annotate("text",x=5,y=8,label=expression(paste("42")),size=4)+
  annotate("text",x=6,y=14,label=expression(paste("61.6")),size=4)+
  annotate("text",x=6,y=8,label=expression(paste("42")),size=4)+
  annotate("text",x=1,y=max(SRCdata[SRCdata$Species=="Avsa","D.PLW1"],na.rm=T)+8,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=2,y=max(SRCdata[SRCdata$Species=="Hovu","D.PLW1"],na.rm=T)+8,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=3,y=max(SRCdata[SRCdata$Species=="Trae","D.PLW1"],na.rm=T)+8,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=4,y=max(SRCdata[SRCdata$Species=="Pegl","D.PLW1"],na.rm=T)+8,label=expression(paste(bold("ab"))),size=5)+
  annotate("text",x=5,y=max(SRCdata[SRCdata$Species=="Zema","D.PLW1"],na.rm=T)+8,label=expression(paste(bold("b"))),size=5)+
  annotate("text",x=6,y=max(SRCdata[SRCdata$Species=="Chga","D.PLW1"],na.rm=T)+8,label=expression(paste(bold("a"))),size=5)
# stat
lm1<-lm(D.PLW1~Species,data=SRCdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups
# PLV (PLW1dry)
g1<-ggplot(SRCdata,aes(x=Species,y=D.PLV,fill=Metabolism))+
  geom_boxplot(lwd=0.2,fatten=5)+
  #geom_point(data=GExdata.0[GExdata.0$Dataset=="Dehydration",],aes(x=Species,y=LMA),size=4,shape=21,fill="yellow")+
  ylim(40,110)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=13),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5,size=13,face="italic"),
        legend.title=element_blank(),
        legend.text=element_text(size=13))+
  ylab(expression(paste("PLV"["dry"]," (%)")))+
  annotate("text",x=1,y=50,label=expression(paste(italic(bar("x"))," = 71.6")),size=4)+
  annotate("text",x=1,y=45,label=expression(paste(italic("n")," = 49")),size=4)+
  annotate("text",x=2,y=50,label=expression(paste("82.5")),size=4)+
  annotate("text",x=2,y=45,label=expression(paste("46")),size=4)+
  annotate("text",x=3,y=50,label=expression(paste("70.5")),size=4)+
  annotate("text",x=3,y=45,label=expression(paste("38")),size=4)+
  annotate("text",x=4,y=50,label=expression(paste("82.7")),size=4)+
  annotate("text",x=4,y=45,label=expression(paste("41")),size=4)+
  annotate("text",x=5,y=50,label=expression(paste("85.3")),size=4)+
  annotate("text",x=5,y=45,label=expression(paste("41")),size=4)+
  annotate("text",x=6,y=50,label=expression(paste("80.2")),size=4)+
  annotate("text",x=6,y=45,label=expression(paste("42")),size=4)+
  annotate("text",x=1,y=max(SRCdata[SRCdata$Species=="Avsa","D.PLV"],na.rm=T)+8,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=2,y=max(SRCdata[SRCdata$Species=="Hovu","D.PLV"],na.rm=T)+8,label=expression(paste(bold("ab"))),size=5)+
  annotate("text",x=3,y=max(SRCdata[SRCdata$Species=="Trae","D.PLV"],na.rm=T)+8,label=expression(paste(bold("c"))),size=5)+
  annotate("text",x=4,y=max(SRCdata[SRCdata$Species=="Pegl","D.PLV"],na.rm=T)+8,label=expression(paste(bold("ab"))),size=5)+
  annotate("text",x=5,y=max(SRCdata[SRCdata$Species=="Zema","D.PLV"],na.rm=T)+8,label=expression(paste(bold("a"))),size=5)+
  annotate("text",x=6,y=max(SRCdata[SRCdata$Species=="Chga","D.PLV"],na.rm=T)+8,label=expression(paste(bold("b"))),size=5)
# stat
lm1<-lm(D.PLV~Species,data=SRCdata)
anova(lm1)
tuk1<-HSD.test(lm1,"Species",group=T)
tuk1$groups

## SRC, scatterplots, stat analysis ####
# LMA, components ####
# LMA vs LT
g1<-ggplot(SRCdata,aes(x=T.LT,y=LMA,colour=Metabolism,shape=Species,lty=Species))+
  geom_point(size=2,aes(alpha=Species))+
  scale_alpha_manual(values=c(1,0.6,0.3,1,0.6,0.3))+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  geom_smooth(size=1.25,alpha=0.8,method=lm,se=F)+
  scale_linetype_manual(values=c(1,2,6,1,2,6))+
  ylim(7,45)+
  xlim(0.06,0.31)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  ylab(expression(paste("LMA (g m"^-2,")")))+
  xlab(expression(paste("LT (mm)")))
# stat
lm1<-lm(LMA~T.LT*Species,data=SRCdata)
anova(lm1)
summary(lm1)
#ph<-glht(lm1,linfct=mcp(Species="Tukey"))
#summary(ph)
# LMA vs LD
g1<-ggplot(SRCdata,aes(x=LD,y=LMA,colour=Metabolism,shape=Species,lty=Species))+
  geom_point(size=2,aes(alpha=Species))+
  scale_alpha_manual(values=c(1,0.6,0.3,1,0.6,0.3))+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  geom_smooth(size=1.25,alpha=0.8,method=lm,se=F)+
  scale_linetype_manual(values=c(1,2,6,1,2,6))+
  ylim(7,45)+
  xlim(0.04,0.35)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  ylab(expression(paste("LMA (g m"^-2,")")))+
  xlab(expression(paste("LD (g cm"^-3,")")))
# stat
lm1<-lm(LMA~LD*Species,data=SRCdata)
anova(lm1)
summary(lm1)
#ph<-glht(lm1,linfct=mcp(Species="Tukey"))
#summary(ph)
# LMA vs LDMC
g1<-ggplot(SRCdata,aes(x=LDMC,y=LMA,colour=Metabolism,shape=Species,lty=Species))+
  geom_point(size=2,aes(alpha=Species))+
  scale_alpha_manual(values=c(1,0.6,0.3,1,0.6,0.3))+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  geom_smooth(size=1.25,alpha=0.8,method=lm,se=F)+
  scale_linetype_manual(values=c(1,2,6,1,2,6))+
  ylim(7,45)+
  xlim(0.04,0.23)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  ylab(expression(paste("LMA (g m"^-2,")")))+
  xlab(expression(paste("LDMC (g g"^-1,")")))
# stat
lm1<-lm(LMA~LDMC*Species,data=SRCdata)
anova(lm1)
summary(lm1)
#ph<-glht(lm1,linfct=mcp(Species="Tukey"))
#summary(ph)
# LD vs LDMC
g1<-ggplot(SRCdata,aes(x=LDMC,y=LD,colour=Metabolism,shape=Species,lty=Species))+
  geom_abline(intercept=0,slope=1,colour="darkgrey")+
  geom_point(size=2,aes(alpha=Species))+
  scale_alpha_manual(values=c(1,0.6,0.3,1,0.6,0.3))+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  geom_smooth(size=1.25,alpha=0.8,method=lm,se=F)+
  scale_linetype_manual(values=c(1,2,6,1,2,6))+
  ylim(0.04,0.35)+
  xlim(0.04,0.23)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  ylab(expression(paste("LD (g cm"^-3,")")))+
  xlab(expression(paste("LDMC (g g"^-1,")")))
# stat
lm1<-lm(LD~LDMC*Species,data=SRCdata)
anova(lm1)
summary(lm1)
#ph<-glht(lm1,linfct=mcp(Species="Tukey"))
#summary(ph)
# LT vs LD
g1<-ggplot(SRCdata,aes(x=LD,y=T.LT,colour=Metabolism,shape=Species,lty=Species))+
  geom_point(size=2,aes(alpha=Species))+
  scale_alpha_manual(values=c(1,0.6,0.3,1,0.6,0.3))+
  geom_smooth(size=1.25,alpha=0.8,method=lm,se=F)+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  scale_linetype_manual(values=c(1,2,6,1,2,6))+
  ylim(0.06,0.31)+
  xlim(0.04,0.35)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  ylab(expression(paste("LT (mm)")))+
  xlab(expression(paste("LD (g cm"^-3,")")))
# stat
lm1<-lm(T.LT~LD*Species,data=SRCdata)
anova(lm1)
summary(lm1)
#ph<-glht(lm1,linfct=mcp(Species="Tukey"))
#summary(ph)
# LT vs LDMC
g1<-ggplot(SRCdata,aes(x=LDMC,y=T.LT,colour=Metabolism,shape=Species,lty=Species))+
  geom_point(size=2,aes(alpha=Species))+
  scale_alpha_manual(values=c(1,0.6,0.3,1,0.6,0.3))+
  geom_smooth(size=1.25,alpha=0.8,method=lm,se=F)+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  scale_linetype_manual(values=c(1,2,6,1,2,6))+
  ylim(0.06,0.31)+
  xlim(0.04,0.23)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  ylab(expression(paste("LT (mm)")))+
  xlab(expression(paste("LDMC (g g"^-1,")")))
# stat
lm1<-lm(T.LT~LDMC*Species,data=SRCdata)
anova(lm1)
summary(lm1)
#ph<-glht(lm1,linfct=mcp(Species="Tukey"))
#summary(ph)

# Shrinkage ####
# PLV vs PLT dry
g1<-ggplot(SRCdata,aes(x=D.PLT,y=D.PLV,colour=Metabolism,shape=Species,lty=Species))+
  geom_point(size=2,aes(alpha=Species))+
  scale_alpha_manual(values=c(1,0.6,0.3,1,0.6,0.3))+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  geom_smooth(size=1.25,alpha=0.8,method=lm,se=F)+
  scale_linetype_manual(values=c(1,2,6,1,2,6))+
  ylim(50,99)+
  xlim(10,90)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  xlab(expression(paste("PLT"[dry]," (%)")))+
  ylab(expression(paste("PLV"[dry]," (%)")))
# stat
lm1<-lm(D.PLV~D.PLT*Species,data=SRCdata)
anova(lm1)
summary(lm1)
#ph<-glht(lm1,linfct=mcp(Species="Tukey"))
#summary(ph)
# PLV vs PLA dry
g1<-ggplot(SRCdata,aes(x=D.PLA,y=D.PLV,colour=Metabolism,shape=Species,lty=Species))+
  geom_point(size=2,aes(alpha=Species))+
  scale_alpha_manual(values=c(1,0.6,0.3,1,0.6,0.3))+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  geom_smooth(size=1.25,alpha=0.8,method=lm,se=F)+
  scale_linetype_manual(values=c(1,2,6,1,2,6))+
  ylim(50,99)+
  xlim(20,80)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  xlab(expression(paste("PLA"[dry]," (%)")))+
  ylab(expression(paste("PLV"[dry]," (%)")))
# stat
lm1<-lm(D.PLV~D.PLA*Species,data=SRCdata)
anova(lm1)
summary(lm1)
#ph<-glht(lm1,linfct=mcp(Species="Tukey"))
#summary(ph)
g1<-ggplot(SRCdata,aes(x=T.LT,y=D.PLT,colour=Metabolism,shape=Species,lty=Species))+
  geom_point(size=2,aes(alpha=Species))+
  scale_alpha_manual(values=c(1,0.6,0.3,1,0.6,0.3))+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  geom_smooth(size=1.25,alpha=0.8,method=lm,se=F)+
  scale_linetype_manual(values=c(1,2,6,1,2,6))+
  ylim(10,90)+
  xlim(0.06,0.31)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  ylab(expression(paste("PLT"[dry]," (%)")))+
  xlab(expression(paste("LT (mm)")))
# stat
lm1<-lm(D.PLT~T.LT*Species,data=SRCdata)
anova(lm1)
summary(lm1)
#ph<-glht(lm1,linfct=mcp(Species="Tukey"))
#summary(ph)
g1<-ggplot(SRCdata,aes(x=T.LW1,y=D.PLW1,colour=Metabolism,shape=Species,lty=Species))+
  geom_point(size=2,aes(alpha=Species))+
  scale_alpha_manual(values=c(1,0.6,0.3,1,0.6,0.3))+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  geom_smooth(size=1.25,alpha=0.8,method=lm,se=F)+
  scale_linetype_manual(values=c(1,2,6,1,2,6))+
  ylim(15,80)+
  xlim(0.1,4.25)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  ylab(expression(paste("PLWi"[dry]," (%)")))+
  xlab(expression(paste("LWi (cm)")))
# stat
lm1<-lm(D.PLW1~T.LW1*Species,data=SRCdata)
anova(lm1)
summary(lm1)
#ph<-glht(lm1,linfct=mcp(Species="Tukey"))
#summary(ph)
g1<-ggplot(SRCdata,aes(x=LD,y=D.PLV,colour=Metabolism,shape=Species,lty=Species))+
  geom_point(size=2,aes(alpha=Species))+
  scale_alpha_manual(values=c(1,0.6,0.3,1,0.6,0.3))+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  geom_smooth(size=1.25,alpha=0.8,method=lm,se=F)+
  scale_linetype_manual(values=c(1,2,6,1,2,6))+
  ylim(50,98)+
  xlim(0.04,0.35)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  ylab(expression(paste("PLV"[dry]," (%)")))+
  xlab(expression(paste("LD (g cm"^-3,")")))
# stat
lm1<-lm(D.PLV~LD*Species,data=SRCdata)
anova(lm1)
summary(lm1)
#ph<-glht(lm1,linfct=mcp(Species="Tukey"))
#summary(ph)
g1<-ggplot(SRCdata,aes(x=LDMC,y=D.PLV,colour=Metabolism,shape=Species,lty=Species))+
  geom_point(size=2,aes(alpha=Species))+
  scale_alpha_manual(values=c(1,0.6,0.3,1,0.6,0.3))+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  geom_smooth(size=1.25,alpha=0.8,method=lm,se=F)+
  scale_linetype_manual(values=c(1,2,6,1,2,6))+
  ylim(50,98)+
  xlim(0.04,0.23)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  ylab(expression(paste("PLV"[dry]," (%)")))+
  xlab(expression(paste("LDMC (g g"^-1,")")))
# stat
lm1<-lm(D.PLV~LDMC*Species,data=SRCdata)
anova(lm1)
summary(lm1)
#ph<-glht(lm1,linfct=mcp(Species="Tukey"))
#summary(ph)
## SRC, dehydration/rehydration analysis ####
# all species, facet ####
# FvFm
g1<-ggplot(data=SRCdata,aes(x=F.FvFm,y=R.FvFm))+
  geom_abline(intercept=0,slope=1,colour="darkgrey")+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_hline(yintercept=0.6,color="orange",lty=2)+
  geom_hline(yintercept=0.35,color="red",lty=2)+
  geom_point(size=3,shape=21,fill="black",alpha=0.5)+
  geom_smooth(aes(x=F.FvFm,y=R.FvFm),color="black",fill="grey")+
  scale_y_continuous(limits=c(-0.25,0.8),breaks=c(0,0.25,0.5,0.75))+
  scale_x_continuous(limits=c(0,0.8),breaks=c(0,0.25,0.5,0.7))+
  #xlim(0,1)+
  facet_wrap(~Species,ncol=3)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        strip.text=element_text(size=13,face="italic"))+
  ylab(expression(paste("rehydration F"[v],"/F"[m])))+
  xlab(expression(paste("dehydration F"[v],"/F"[m])))
# water potential - 'mean'
g1<-ggplot(data=SRCdata,aes(x=F.RWC,y=WPleaf2))+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_point(size=3,shape=21,fill="black",alpha=0.5)+
  geom_smooth(aes(x=F.RWC,y=WPleaf2),color="black",fill="grey")+
  #scale_y_continuous(limits=c(-0.25,0.8),breaks=c(0,0.25,0.5,0.75))+
  scale_x_continuous(limits=c(0,100))+
  #xlim(0,1)+
  facet_wrap(~Species,ncol=3)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        strip.text=element_text(size=13,face="italic"))+
  ylab(expression(paste(Psi["leaf"]," (MPa)")))+
  xlab(expression(paste("RWC (%)")))
# data per species, FvFm, PLRC ####
# Avsa
g1<-ggplot(SRCdata,aes(x=F.RWC,y=F.FvFm))+
  geom_rect(data=data.frame(xmin=16,xmax=75,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_hline(yintercept=mean(SRCdata[SRCdata$Species=="Avsa","T.FvFm"],na.rm=T),color="green1",lty=2)+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_vline(xintercept=89.5,color="black",lty=2)+
  geom_vline(xintercept=77.8,color="darkgrey",lty=2)+
  geom_vline(xintercept=46.1,color="black",lty=3)+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_hline(yintercept=0.6,color="orange",lty=2)+
  geom_hline(yintercept=0.35,color="red",lty=2)+
  geom_hline(yintercept=0.1,color="skyblue1",lty=2)+
  geom_hline(yintercept=0.25,color="skyblue2",lty=2)+
  geom_hline(yintercept=0.5,color="skyblue3",lty=2)+
  geom_point(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=F.FvFm),size=2,shape=21,alpha=0.5,fill="green1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=F.FvFm),color="green1",fill="palegreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=R.FvFm),size=2,shape=21,alpha=0.5,fill="darkgoldenrod1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=R.FvFm),color="darkgoldenrod1",fill="lightgoldenrod1")+
  geom_point(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=PLRC/100),size=2,shape=21,alpha=0.5,fill="lightskyblue1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=PLRC/100),color="royalblue1",fill="lightskyblue1")+
  ylim(-0.25,1)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=9))+
  ylab(expression(paste("PLRC; dehy, rehy F"[v],"/F"[m])))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=15,y=0.97,label=expression(paste(italic("Avsa"))))
# Hovu
g1<-ggplot(SRCdata,aes(x=F.RWC,y=F.FvFm))+
  geom_rect(data=data.frame(xmin=22,xmax=74,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
    alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_rect(data=data.frame(xmin=23,xmax=60,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
    alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_rect(data=data.frame(xmin=19,xmax=74,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
    alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_hline(yintercept=mean(SRCdata[SRCdata$Species=="Hovu","T.FvFm"],na.rm=T),color="green1",lty=2)+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_vline(xintercept=90.7,color="black",lty=2)+
  geom_vline(xintercept=84.0,color="darkgrey",lty=2)+
  geom_vline(xintercept=66.3,color="black",lty=3)+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_hline(yintercept=0.6,color="orange",lty=2)+
  geom_hline(yintercept=0.35,color="red",lty=2)+
  geom_hline(yintercept=0.1,color="skyblue1",lty=2)+
  geom_hline(yintercept=0.25,color="skyblue2",lty=2)+
  geom_hline(yintercept=0.5,color="skyblue3",lty=2)+
  geom_point(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=F.FvFm),size=2,shape=21,alpha=0.5,fill="green1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=F.FvFm),color="green1",fill="palegreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=R.FvFm),size=2,shape=21,alpha=0.5,fill="darkgoldenrod1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=R.FvFm),color="darkgoldenrod1",fill="lightgoldenrod1")+
  geom_point(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=PLRC/100),size=2,shape=21,alpha=0.5,fill="lightskyblue1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=PLRC/100),color="royalblue1",fill="lightskyblue1")+
  ylim(-0.25,1)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=9))+
  ylab(expression(paste("PLRC; dehy, rehy F"[v],"/F"[m])))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=15,y=0.97,label=expression(paste(italic("Hovu"))))
# Trae
g1<-ggplot(SRCdata,aes(x=F.RWC,y=F.FvFm))+
  geom_rect(data=data.frame(xmin=11.5,xmax=61.5,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_rect(data=data.frame(xmin=23,xmax=70,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_rect(data=data.frame(xmin=44,xmax=72,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_hline(yintercept=mean(SRCdata[SRCdata$Species=="Trae","T.FvFm"],na.rm=T),color="green1",lty=2)+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_vline(xintercept=90.0,color="black",lty=2)+
  geom_vline(xintercept=81.3,color="darkgrey",lty=2)+
  geom_vline(xintercept=61.2,color="black",lty=3)+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_hline(yintercept=0.6,color="orange",lty=2)+
  geom_hline(yintercept=0.35,color="red",lty=2)+
  geom_hline(yintercept=0.1,color="skyblue1",lty=2)+
  geom_hline(yintercept=0.25,color="skyblue2",lty=2)+
  geom_hline(yintercept=0.5,color="skyblue3",lty=2)+
  geom_point(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=F.FvFm),size=2,shape=21,alpha=0.5,fill="green1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=F.FvFm),color="green1",fill="palegreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=R.FvFm),size=2,shape=21,alpha=0.5,fill="darkgoldenrod1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=R.FvFm),color="darkgoldenrod1",fill="lightgoldenrod1")+
  geom_point(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=PLRC/100),size=2,shape=21,alpha=0.5,fill="lightskyblue1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=PLRC/100),color="royalblue1",fill="lightskyblue1")+
  ylim(-0.25,1)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=9))+
  ylab(expression(paste("PLRC; dehy, rehy F"[v],"/F"[m])))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=15,y=0.97,label=expression(paste(italic("Trae"))))
# Pegl
g1<-ggplot(SRCdata,aes(x=F.RWC,y=F.FvFm))+
  geom_rect(data=data.frame(xmin=7,xmax=9,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_rect(data=data.frame(xmin=20,xmax=30,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_hline(yintercept=mean(SRCdata[SRCdata$Species=="Pegl","T.FvFm"],na.rm=T),color="green1",lty=2)+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_vline(xintercept=91.2,color="black",lty=2)+
  geom_vline(xintercept=86.2,color="darkgrey",lty=2)+
  geom_vline(xintercept=83.8,color="black",lty=3)+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_hline(yintercept=0.6,color="orange",lty=2)+
  geom_hline(yintercept=0.35,color="red",lty=2)+
  geom_hline(yintercept=0.1,color="skyblue1",lty=2)+
  geom_hline(yintercept=0.25,color="skyblue2",lty=2)+
  geom_hline(yintercept=0.5,color="skyblue3",lty=2)+
  geom_point(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=F.FvFm),size=2,shape=21,alpha=0.5,fill="green1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=F.FvFm),color="green1",fill="palegreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=R.FvFm),size=2,shape=21,alpha=0.5,fill="darkgoldenrod1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=R.FvFm),color="darkgoldenrod1",fill="lightgoldenrod1")+
  geom_point(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=PLRC/100),size=2,shape=21,alpha=0.5,fill="lightskyblue1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=PLRC/100),color="royalblue1",fill="lightskyblue1")+
  ylim(-0.25,1)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=9))+
  ylab(expression(paste("PLRC; dehy, rehy F"[v],"/F"[m])))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=15,y=0.97,label=expression(paste(italic("Pegl"))))
# Zema
g1<-ggplot(SRCdata,aes(x=F.RWC,y=F.FvFm))+
  geom_rect(data=data.frame(xmin=4,xmax=5,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_rect(data=data.frame(xmin=11,xmax=35,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_hline(yintercept=mean(SRCdata[SRCdata$Species=="Zema","T.FvFm"],na.rm=T),color="green1",lty=2)+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_vline(xintercept=91.8,color="black",lty=2)+
  #geom_vline(xintercept=86.2,color="darkgrey",lty=2)+
  #geom_vline(xintercept=83.8,color="black",lty=3)+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_hline(yintercept=0.6,color="orange",lty=2)+
  geom_hline(yintercept=0.35,color="red",lty=2)+
  geom_hline(yintercept=0.1,color="skyblue1",lty=2)+
  geom_hline(yintercept=0.25,color="skyblue2",lty=2)+
  geom_hline(yintercept=0.5,color="skyblue3",lty=2)+
  geom_point(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=F.FvFm),size=2,shape=21,alpha=0.5,fill="green1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=F.FvFm),color="green1",fill="palegreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=R.FvFm),size=2,shape=21,alpha=0.5,fill="darkgoldenrod1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=R.FvFm),color="darkgoldenrod1",fill="lightgoldenrod1")+
  geom_point(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=PLRC/100),size=2,shape=21,alpha=0.5,fill="lightskyblue1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=PLRC/100),color="royalblue1",fill="lightskyblue1")+
  ylim(-0.25,1)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=9))+
  ylab(expression(paste("PLRC; dehy, rehy F"[v],"/F"[m])))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=15,y=0.97,label=expression(paste(italic("Zema"))))
# Chga
g1<-ggplot(SRCdata,aes(x=F.RWC,y=F.FvFm))+
  geom_rect(data=data.frame(xmin=9,xmax=10,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_rect(data=data.frame(xmin=9,xmax=10,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_rect(data=data.frame(xmin=7,xmax=40,ymin=0,ymax=1),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            alpha=0.05,fill="royalblue",inherit.aes=F)+
  geom_hline(yintercept=mean(SRCdata[SRCdata$Species=="Chga","T.FvFm"],na.rm=T),color="green1",lty=2)+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_vline(xintercept=93.1,color="black",lty=2)+
  geom_vline(xintercept=82.8,color="darkgrey",lty=2)+
  geom_vline(xintercept=77.5,color="black",lty=3)+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_hline(yintercept=0.6,color="orange",lty=2)+
  geom_hline(yintercept=0.35,color="red",lty=2)+
  geom_hline(yintercept=0.1,color="skyblue1",lty=2)+
  geom_hline(yintercept=0.25,color="skyblue2",lty=2)+
  geom_hline(yintercept=0.5,color="skyblue3",lty=2)+
  geom_point(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=F.FvFm),size=2,shape=21,alpha=0.5,fill="green1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=F.FvFm),color="green1",fill="palegreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=R.FvFm),size=2,shape=21,alpha=0.5,fill="darkgoldenrod1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=R.FvFm),color="darkgoldenrod1",fill="lightgoldenrod1")+
  geom_point(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=PLRC/100),size=2,shape=21,alpha=0.5,fill="lightskyblue1")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=PLRC/100),color="royalblue1",fill="lightskyblue1")+
  ylim(-0.25,1)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=9))+
  ylab(expression(paste("PLRC; dehy, rehy F"[v],"/F"[m])))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=15,y=0.97,label=expression(paste(italic("Chga"))))

# data per species, Shrinkage ####
# Avsa
g1<-ggplot()+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_vline(xintercept=89.5,color="black",lty=2)+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_hline(yintercept=10,color="grey",lty=2)+
  geom_hline(yintercept=25,color="darkgrey",lty=2)+
  geom_hline(yintercept=50,color="black",lty=2)+
  geom_point(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=F.PLV),size=2,shape=21,alpha=0.5,fill="forestgreen")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=F.PLV),color="forestgreen",fill="forestgreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=F.PLT),size=2,shape=21,alpha=0.5,fill="chartreuse")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=F.PLT),color="chartreuse3",fill="chartreuse")+
  geom_point(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=F.PLA),size=2,shape=21,alpha=0.5,fill="palegreen")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=F.PLA),color="chartreuse2",fill="palegreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=F.PLW1),size=2,shape=21,alpha=0.5,fill="white",colour="chartreuse2")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Avsa",],aes(x=F.RWC,y=F.PLW1),color="chartreuse2",fill="palegreen",lty=2)+
  ylim(-10,100)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=9))+
  ylab(expression(paste("Shrinkage; PLV,PLT,PLA,PLWi (%)")))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=75,y=97,label=expression(paste(italic("Avsa"))))
# Hovu
g1<-ggplot()+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_vline(xintercept=90.7,color="black",lty=2)+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_hline(yintercept=10,color="grey",lty=2)+
  geom_hline(yintercept=25,color="darkgrey",lty=2)+
  geom_hline(yintercept=50,color="black",lty=2)+
  geom_point(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=F.PLV),size=2,shape=21,alpha=0.5,fill="forestgreen")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=F.PLV),color="forestgreen",fill="forestgreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=F.PLT),size=2,shape=21,alpha=0.5,fill="chartreuse")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=F.PLT),color="chartreuse3",fill="chartreuse")+
  geom_point(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=F.PLA),size=2,shape=21,alpha=0.5,fill="palegreen")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=F.PLA),color="chartreuse2",fill="palegreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=F.PLW1),size=2,shape=21,alpha=0.5,fill="white",colour="chartreuse2")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Hovu",],aes(x=F.RWC,y=F.PLW1),color="chartreuse2",fill="palegreen",lty=2)+
  ylim(-10,100)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=9))+
  ylab(expression(paste("Shrinkage; PLV,PLT,PLA,PLWi (%)")))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=75,y=97,label=expression(paste(italic("Hovu"))))
# Trae
g1<-ggplot()+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_vline(xintercept=90.0,color="black",lty=2)+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_hline(yintercept=10,color="grey",lty=2)+
  geom_hline(yintercept=25,color="darkgrey",lty=2)+
  geom_hline(yintercept=50,color="black",lty=2)+
  geom_point(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=F.PLV),size=2,shape=21,alpha=0.5,fill="forestgreen")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=F.PLV),color="forestgreen",fill="forestgreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=F.PLT),size=2,shape=21,alpha=0.5,fill="chartreuse")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=F.PLT),color="chartreuse3",fill="chartreuse")+
  geom_point(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=F.PLA),size=2,shape=21,alpha=0.5,fill="palegreen")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=F.PLA),color="chartreuse2",fill="palegreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=F.PLW1),size=2,shape=21,alpha=0.5,fill="white",colour="chartreuse2")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Trae",],aes(x=F.RWC,y=F.PLW1),color="chartreuse2",fill="palegreen",lty=2)+
  ylim(-10,100)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=9))+
  ylab(expression(paste("Shrinkage; PLV,PLT,PLA,PLWi (%)")))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=75,y=97,label=expression(paste(italic("Trae"))))
# Pegl
g1<-ggplot()+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_vline(xintercept=91.2,color="black",lty=2)+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_hline(yintercept=10,color="grey",lty=2)+
  geom_hline(yintercept=25,color="darkgrey",lty=2)+
  geom_hline(yintercept=50,color="black",lty=2)+
  geom_point(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=F.PLV),size=2,shape=21,alpha=0.5,fill="forestgreen")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=F.PLV),color="forestgreen",fill="forestgreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=F.PLT),size=2,shape=21,alpha=0.5,fill="chartreuse")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=F.PLT),color="chartreuse3",fill="chartreuse")+
  geom_point(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=F.PLA),size=2,shape=21,alpha=0.5,fill="palegreen")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=F.PLA),color="chartreuse2",fill="palegreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=F.PLW1),size=2,shape=21,alpha=0.5,fill="white",colour="chartreuse2")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Pegl",],aes(x=F.RWC,y=F.PLW1),color="chartreuse2",fill="palegreen",lty=2)+
  ylim(-10,100)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=9))+
  ylab(expression(paste("Shrinkage; PLV,PLT,PLA,PLWi (%)")))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=75,y=97,label=expression(paste(italic("Pegl"))))
# Zema
g1<-ggplot()+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_vline(xintercept=91.8,color="black",lty=2)+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_hline(yintercept=10,color="grey",lty=2)+
  geom_hline(yintercept=25,color="darkgrey",lty=2)+
  geom_hline(yintercept=50,color="black",lty=2)+
  geom_point(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=F.PLV),size=2,shape=21,alpha=0.5,fill="forestgreen")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=F.PLV),color="forestgreen",fill="forestgreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=F.PLT),size=2,shape=21,alpha=0.5,fill="chartreuse")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=F.PLT),color="chartreuse3",fill="chartreuse")+
  geom_point(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=F.PLA),size=2,shape=21,alpha=0.5,fill="palegreen")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=F.PLA),color="chartreuse2",fill="palegreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=F.PLW1),size=2,shape=21,alpha=0.5,fill="white",colour="chartreuse2")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Zema",],aes(x=F.RWC,y=F.PLW1),color="chartreuse2",fill="palegreen",lty=2)+
  ylim(-10,100)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=9))+
  ylab(expression(paste("Shrinkage; PLV,PLT,PLA,PLWi (%)")))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=75,y=97,label=expression(paste(italic("Zema"))))
# Chga
g1<-ggplot()+
  geom_vline(xintercept=0,color="darkgrey")+
  geom_vline(xintercept=93.1,color="black",lty=2)+
  geom_hline(yintercept=0,color="darkgrey")+
  geom_hline(yintercept=10,color="grey",lty=2)+
  geom_hline(yintercept=25,color="darkgrey",lty=2)+
  geom_hline(yintercept=50,color="black",lty=2)+
  geom_point(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=F.PLV),size=2,shape=21,alpha=0.5,fill="forestgreen")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=F.PLV),color="forestgreen",fill="forestgreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=F.PLT),size=2,shape=21,alpha=0.5,fill="chartreuse")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=F.PLT),color="chartreuse3",fill="chartreuse")+
  geom_point(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=F.PLA),size=2,shape=21,alpha=0.5,fill="palegreen")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=F.PLA),color="chartreuse2",fill="palegreen")+
  geom_point(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=F.PLW1),size=2,shape=21,alpha=0.5,fill="white",colour="chartreuse2")+
  geom_smooth(data=SRCdata[SRCdata$Species=="Chga",],aes(x=F.RWC,y=F.PLW1),color="chartreuse2",fill="palegreen",lty=2)+
  ylim(-10,100)+
  xlim(0,100)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=9))+
  ylab(expression(paste("Shrinkage; PLV,PLT,PLA,PLWi (%)")))+
  xlab(expression(paste("RWC (%)")))+
  annotate("text",x=75,y=97,label=expression(paste(italic("Chga"))))



## X50, correlations ####
d1<-X50data.RWC[c("LMA","LT","LD","LDMC","RWCtlp",
                "An50","gsw50","phiPSII50","dehyF0.6","dehyF0.35","PLRC10","PLRC25",
                "PLT10","PLT25","PLTdry","PLA10","PLA25","PLAdry","PLV10","PLV25","PLVdry")]
m1<-cor(d1,use="pairwise.complete.obs")
g1<-corrplot.mixed(m1)
# most interesting, clear relationships
# RWC,dehyF0.35 vs RWC,PLV25
g1<-ggplot(X50data.RWC,aes(x=PLV25,y=dehyF0.35))+
  geom_abline(intercept=0,slope=1,colour="darkgrey")+
  geom_smooth(size=1,alpha=0.8,method=lm,se=T,colour="darkgrey",fill="lightgrey")+
  geom_point(size=3,aes(colour=Metabolism,shape=Species))+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  ylim(-5,55)+
  xlim(60,85)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  ylab(expression(paste("RWC"["dehyF=0.35"]," (%)")))+
  xlab(expression(paste("RWC"["PLV=25"]," (%)")))
# stat
lm1<-lm(dehyF0.35~PLV25,data=X50data.RWC)
anova(lm1)
summary(lm1)
# RWC,dehyF0.60 vs RWC,PLA25
g1<-ggplot(X50data.RWC,aes(x=PLA25,y=dehyF0.6))+
  geom_abline(intercept=0,slope=1,colour="darkgrey")+
  geom_smooth(size=1,alpha=0.8,method=lm,se=T,colour="darkgrey",fill="lightgrey")+
  geom_point(size=3,aes(colour=Metabolism,shape=Species))+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  ylim(10,65)+
  xlim(25,45)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  ylab(expression(paste("RWC"["dehyF=0.6"]," (%)")))+
  xlab(expression(paste("RWC"["PLA=25"]," (%)")))
# stat
lm1<-lm(dehyF0.6~PLA25,data=X50data.RWC)
anova(lm1)
summary(lm1)
# RWC,PLT25 vs LDMC
g1<-ggplot(X50data.RWC,aes(x=LDMC,y=PLT25))+
  geom_abline(intercept=0,slope=1,colour="darkgrey")+
  geom_smooth(size=1,alpha=0.8,method=lm,se=T,colour="darkgrey",fill="lightgrey")+
  geom_point(size=3,aes(colour=Metabolism,shape=Species))+
  scale_shape_manual(values=c(16,15,17,16,15,17))+
  ylim(40,85)+
  xlim(0.075,0.175)+
  #geom_hline(yintercept=0,lty=2)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=13),
        axis.text.x=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'cm'))+
  ylab(expression(paste("RWC"["PLT=25"]," (%)")))+
  xlab(expression(paste("LDMC (g g"^-1,")")))
# stat
lm1<-lm(PLT25~LDMC,data=X50data.RWC)
anova(lm1)
summary(lm1)

## X50, PCA ####
# complete obs
d1<-X50data.RWC[c("Species","Metabolism",
                  "dehyF0.6","dehyF0.35","PLRC10","PLRC25",
                  "PLT10","PLT25","PLTdry","PLA10","PLA25","PLAdry","PLV10","PLV25","PLVdry")]
pc1<-prcomp(scale(d1[c(3:15)]))
met<-as.factor(d1$Metabolism[1:6])
sp<-as.factor(d1$Species[1:6])
# 
pca1<-fviz_pca_biplot(pc1,label="var",col.var="black",repel=T,
                      geom.ind="point",pointsize=5,col.ind=met,invisible="quali")+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(size=1,colour="#222222"),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        legend.position="right",
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        legend.background=element_blank(),
        title=element_blank())+
  ylab(expression(paste("PC2 (28.7 %)")))+
  xlab(expression(paste("PC1 (54.1 %)")))
