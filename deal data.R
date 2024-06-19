########药物靶点大于1######
library(vroom)
setwd("D:\\database\\NEW\\01 deal\\gene大于1")
drug<-data.frame(vroom("drug.txt",delim = "\t"))
int<-data.frame(vroom("int.txt",delim = "\t"))
drug[2,]
for (i in 1:18085) {
  for(j in 1:10141)
  
}

###bindingdb####
#############BindingDB###########
BindingDB_data<-read.table("E:\\迅雷下载\\BindingDB_All.tsv",header = F,
                           sep ="\t",quote="",stringsAsFactors = F,
                           na.strings = c("NA"),fill=T,fileEncoding = "")

colnames(BindingDB_data) <- BindingDB_data[1,]
file<-BindingDB_data[-1,c(1,9,10,11,12)]
######提取Ki<10########
rm(file1)
ki<-file[,2]
ki1<-as.numeric(ki)
index_Ki<-which(ki1<=10000)
####IC50####
IC50<-file[,3]
IC501<-as.numeric(IC50)
index_IC<-which(IC501<=10000)
#######Kd#######
Kd<-file[,4]
Kd1<-as.numeric(Kd)
index_Kd<-which(Kd1<=10000)
########EC50########
EC50<-file[,5]
EC501<-as.numeric(EC50)
index_EC50<-which(EC501<=10000)
#######取并集###########
a<-union(index_EC50,index_IC)
b<-union(index_Kd,index_Ki)
select<-union(a,b)
select
file<-BindingDB_data[c(select),]
unique(file[,6])
write.csv(file,file = "D:\\subject\\01 database\\05.26 deal\\BindingDB\\int.csv",
          sep = ",",row.names = F)
###选择人类###
homo<-file[grep(pattern="Homo sapiens",file[,8]),]
unique(homo[,6])
homo<-homo[,c(2,3,5,6,9,10,11,12,29,42)]
write.csv(homo,file = "D:\\subject\\01 database\\05.26 deal\\BindingDB\\int.csv",
          sep = ",",row.names = F)
