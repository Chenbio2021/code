#先做DDC_path10   train的矩阵#
getwd()
setwd("/boot2/chenjiaqi/R/87drug/path10")
library(vroom)
train10<-vroom("NCI_train38.csv")
train10<-as.data.frame(train10)[-1]
colnames(train10)[9]<-c("drugname1")
colnames(train10)[10]<-c("drugname2")
#DD_maccs#
DD_maccs<-vroom("DD_maccs.csv")
DD_maccs<-as.data.frame(DD_maccs)
row.names(DD_maccs)<-DD_maccs[,1]
DD_maccs<-DD_maccs[,-1]
colnames(DD_maccs)[3:111]<-paste0("DD_maccs",sep="_",colnames(DD_maccs)[3:111])

#DD_tox#
DD_tox<-vroom("DD_tox.csv")
DD_tox<-as.data.frame(DD_tox)
row.names(DD_tox)<-DD_tox[,1]
DD_tox<-DD_tox[,-1]
colnames(DD_tox)[3:99]<-paste0("DD_tox",sep="_",colnames(DD_tox)[3:99])

#DDP10#
DDP<-vroom("DDP10.csv")
DDP<-as.data.frame(DDP)
row.names(DDP)<-DDP[,1]
DDP<-DDP[,-1]
colnames(DDP)[3:159]<-paste0("DDP",sep="_",colnames(DDP)[3:159])

#exp#
CP_EXP<-vroom("exp10.csv")
CP_EXP<-as.data.frame(CP_EXP)
row.names(CP_EXP)<-CP_EXP[,1]
colnames(CP_EXP)<-paste0("exp",sep="_",colnames(CP_EXP))
colnames(CP_EXP)[1]<-c("depmap_id")
#methy#
CP_METHY<-vroom("methy10.csv")
CP_METHY<-as.data.frame(CP_METHY)
row.names(CP_METHY)<-CP_METHY[,1]
colnames(CP_METHY)<-paste0("methy",sep="_",colnames(CP_METHY))
colnames(CP_METHY)[1]<-c("depmap_id")
#cnv#
CP_CNV<-vroom("cnv10.csv")
CP_CNV<-as.data.frame(CP_CNV)
row.names(CP_CNV)<-CP_CNV[,1]
colnames(CP_CNV)<-paste0("cnv",sep="_",colnames(CP_CNV))
colnames(CP_CNV)[1]<-c("depmap_id")
#mut#
CP_MUT<-vroom("mut10.csv")
CP_MUT<-as.data.frame(CP_MUT)
row.names(CP_MUT)<-CP_MUT[,1]
colnames(CP_MUT)<-paste0("mut",sep="_",colnames(CP_MUT))
colnames(CP_MUT)[1]<-c("depmap_id")
#RNAi#
CP_RNAi<-vroom("RNAi10.csv")
CP_RNAi<-as.data.frame(CP_RNAi)
row.names(CP_RNAi)<-CP_RNAi[,1]
colnames(CP_RNAi)<-paste0("RNAi",sep="_",colnames(CP_RNAi))
colnames(CP_RNAi)[1]<-c("depmap_id")

#input 矩阵   DDP必要#
#train10的药物特征merge#
data1<-merge(train10,DD_maccs,by=c("drugname1","drugname2"))
data2<-merge(data1,DD_tox,by=c("drugname1","drugname2"))
data3<-merge(data2,DDP,by=c("drugname1","drugname2"))
#cell line feature先合并在mergedata3
data4<-merge(CP_CNV,CP_EXP,by="depmap_id")
data5<-merge(data4,CP_METHY,by="depmap_id")
data6<-merge(data5,CP_MUT,by="depmap_id")
data7<-merge(data6,CP_RNAi,by="depmap_id")
#所有的merge#
data8<-merge(data3,data7,by="depmap_id")
write.csv(data8,"train10.csv",row.names = T)

#DDC_path10   test的矩阵#
getwd()
setwd("/boot2/chenjiaqi/R/87drug/path10")
library(vroom)
test10<-vroom("NCI_test38.csv")
test10<-as.data.frame(test10)[-1]
colnames(test10)[5]<-c("depmap_id")

#test 矩阵#
#test10的药物特征merge#
data1<-merge(test10,DD_maccs,by=c("drugname1","drugname2"))
data2<-merge(data1,DD_tox,by=c("drugname1","drugname2"))
data3<-merge(data2,DDP,by=c("drugname1","drugname2"))
#cell line feature先合并在mergedata3
data4<-merge(CP_CNV,CP_EXP,by="depmap_id")
data5<-merge(data4,CP_METHY,by="depmap_id")
data6<-merge(data5,CP_MUT,by="depmap_id")
data7<-merge(data6,CP_RNAi,by="depmap_id")
#所有的merge#
data8<-merge(data3,data7,by="depmap_id")
write.csv(data8,"test10.csv",row.names = T)






# ####################先做DDC_path157   train的矩阵###################
# getwd()
# setwd("/boot2/chenjiaqi/R/87drug/path10")
# library(vroom)
# train10<-vroom("NCI_train38.csv")
# train10<-as.data.frame(train10)[-1]
# colnames(train10)[9]<-c("drugname1")
# colnames(train10)[10]<-c("drugname2")
# #DD_maccs#
# DD_maccs<-vroom("DD_maccs.csv")
# DD_maccs<-as.data.frame(DD_maccs)
# row.names(DD_maccs)<-DD_maccs[,1]
# DD_maccs<-DD_maccs[,-1]
# colnames(DD_maccs)[3:111]<-paste0("DD_maccs",sep="_",colnames(DD_maccs)[3:111])
# 
# #DD_tox#
# DD_tox<-vroom("DD_tox.csv")
# DD_tox<-as.data.frame(DD_tox)
# row.names(DD_tox)<-DD_tox[,1]
# DD_tox<-DD_tox[,-1]
# colnames(DD_tox)[3:99]<-paste0("DD_tox",sep="_",colnames(DD_tox)[3:99])
# 
# #DDP10#
# DDP<-vroom("DDP(157).csv")
# DDP<-as.data.frame(DDP)
# row.names(DDP)<-DDP[,1]
# DDP<-DDP[,-1]
# colnames(DDP)[3:159]<-paste0("DDP",sep="_",colnames(DDP)[3:159])
# 
# #exp#
# CP_EXP<-vroom("exp10.csv")
# CP_EXP<-as.data.frame(CP_EXP)
# row.names(CP_EXP)<-CP_EXP[,1]
# colnames(CP_EXP)<-paste0("exp",sep="_",colnames(CP_EXP))
# colnames(CP_EXP)[1]<-c("depmap_id")
# #methy#
# CP_METHY<-vroom("methy10.csv")
# CP_METHY<-as.data.frame(CP_METHY)
# row.names(CP_METHY)<-CP_METHY[,1]
# colnames(CP_METHY)<-paste0("methy",sep="_",colnames(CP_METHY))
# colnames(CP_METHY)[1]<-c("depmap_id")
# #cnv#
# CP_CNV<-vroom("cnv10.csv")
# CP_CNV<-as.data.frame(CP_CNV)
# row.names(CP_CNV)<-CP_CNV[,1]
# colnames(CP_CNV)<-paste0("cnv",sep="_",colnames(CP_CNV))
# colnames(CP_CNV)[1]<-c("depmap_id")
# #mut#
# CP_MUT<-vroom("mut10.csv")
# CP_MUT<-as.data.frame(CP_MUT)
# row.names(CP_MUT)<-CP_MUT[,1]
# colnames(CP_MUT)<-paste0("mut",sep="_",colnames(CP_MUT))
# colnames(CP_MUT)[1]<-c("depmap_id")
# #RNAi#
# CP_RNAi<-vroom("RNAi10.csv")
# CP_RNAi<-as.data.frame(CP_RNAi)
# row.names(CP_RNAi)<-CP_RNAi[,1]
# colnames(CP_RNAi)<-paste0("RNAi",sep="_",colnames(CP_RNAi))
# colnames(CP_RNAi)[1]<-c("depmap_id")
# 
# #input 矩阵   DDP必要#
# #train10的药物特征merge#
# data1<-merge(train10,DD_maccs,by=c("drugname1","drugname2"))
# data2<-merge(data1,DD_tox,by=c("drugname1","drugname2"))
# data3<-merge(data2,DDP,by=c("drugname1","drugname2"))
# #cell line feature先合并在mergedata3
# data4<-merge(CP_CNV,CP_EXP,by="depmap_id")
# data5<-merge(data4,CP_METHY,by="depmap_id")
# data6<-merge(data5,CP_MUT,by="depmap_id")
# data7<-merge(data6,CP_RNAi,by="depmap_id")
# #所有的merge#
# data8<-merge(data3,data7,by="depmap_id")
# write.csv(data8,"train_DDP157.csv",row.names = T)
# 
# #DDC_path10   test的矩阵#
# getwd()
# setwd("/boot2/chenjiaqi/R/87drug/path10")
# library(vroom)
# test10<-vroom("NCI_test38.csv")
# test10<-as.data.frame(test10)[-1]
# colnames(test10)[5]<-c("depmap_id")
# 
# #test 矩阵#
# #test10的药物特征merge#
# data1<-merge(test10,DD_maccs,by=c("drugname1","drugname2"))
# data2<-merge(data1,DD_tox,by=c("drugname1","drugname2"))
# data3<-merge(data2,DDP,by=c("drugname1","drugname2"))
# #cell line feature先合并在mergedata3
# data4<-merge(CP_CNV,CP_EXP,by="depmap_id")
# data5<-merge(data4,CP_METHY,by="depmap_id")
# data6<-merge(data5,CP_MUT,by="depmap_id")
# data7<-merge(data6,CP_RNAi,by="depmap_id")
# #所有的merge#
# data8<-merge(data3,data7,by="depmap_id")
# write.csv(data8,"test_DDP157.csv",row.names = T)
