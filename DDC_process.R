#数据预处理#
#导入NCI数据#
library(vroom)
NCI<-vroom("D://subject//02 Drug prediction//dataset//02 article//data//ComboDrugGrowth_Nov2017.csv")

NCI<-as.data.frame(NCI)
NCI60<-NCI[,c(9,15,26,28,29)]
#删除含有na的行
NCI60<-na.omit(NCI60)
NSC<-NCI60[,c(1,2)]
#NSC1和NSC2不考虑顺序问题,小的在前，大的在后
NSC <- data.frame(t(apply(NSC, 1, function(x) {
  if (x[1] > x[2]) {
    c(x[2], x[1])
  } else {
    x
  }
})))
NCI60[,c(1,2)]<-NSC
#计算SCORE平均值#
library(dplyr)
NCI60 <- NCI60 %>%
  group_by(NSC1, NSC2, PANEL, CELLNAME) %>%
  summarize(mean_score = mean(SCORE, na.rm = TRUE))
colnames(NCI60)<-c("drug1","drug2","cancer","cell line","score")
#按照得分为10，将结果区分为敏感（>10），不敏感（<=10），敏感为1，不敏感为0
NCI1<-cut(NCI60$score,breaks = c(-Inf,10,Inf),labels = c(0,1),right = TRUE)
NCI60[,6]<-NCI1
#按照得分为10，将结果区分为敏感（>4），不敏感（<=4），敏感为1，不敏感为0
NCI2<-cut(NCI60$score,breaks = c(-Inf,4,Inf),labels = c(0,1),right = TRUE)
NCI60[,7]<-NCI2
colnames(NCI60)<-c("drug1","drug2","cancer","cell line","score","cut10","cut4")
#药物组合信息#
NCI60$drug_combination <- paste(NCI60$drug1, NCI60$drug2, sep = "-")
#####附加上药物名称#####
head(NCI60)
drugname<-read.csv("D://subject//02 Drug prediction//dataset//02 article//data//ComboCompoundNames_small.csv",header = T)
colnames(drugname)<-c("drug1","drug1name")
a<-left_join(NCI60,drugname,by="drug1")
colnames(drugname)<-c("drug2","drug2name")
b<-left_join(a,drugname,by="drug2")
NCI60<-b
#由于753082和761431是一样的，且全部名称中753082给的更少，删除这些行#
a<-grep("753082", NCI60$drug_combination)
NCI60<-NCI60[-a,]
#以上是104种drug的#
#先制作药物-靶点01矩阵，并对药物进行筛选#
#drug-gene01矩阵制作#
##在excel中完成药物靶点的匹配##
drug_gene<-vroom("D://subject//02 Drug prediction//dataset//02 article//data//drug_gene_interaction.csv")
#####104个NCI药物的01药靶矩阵######
#取基因名称（行名）#
gene_name<-unique(drug_gene$Gene)
#取药物名称（列名）#
drug_name<-unique(drug_gene$Drug)

#建立空表#
drug_gene_int<-matrix(nrow=1583,ncol=104)
row.names(drug_gene_int)<-gene_name
colnames(drug_gene_int)<-drug_name
#建立01矩阵#
for (i in 1:5155) {
  drug<-drug_gene[i,1]
  gene<-drug_gene[i,2]
  drug_gene_int[match(gene,rownames(drug_gene_int)),match(drug,colnames(drug_gene_int))]<-1
}

drug_gene_int[is.na(drug_gene_int)]<-0
write.csv(drug_gene_int,"drug_gene_01矩阵.csv")
###########对药物进行筛选##剩下87种并建立01矩阵############
library(vroom)
drug_gene<-vroom("D://subject//02 Drug prediction//dataset//02 article//data//drug_gene_01矩阵.csv")
drug_gene<-as.data.frame(drug_gene)
rownames(drug_gene)<-drug_gene[,1]
drug_gene<-drug_gene[,-1]
sum_row <- colSums(drug_gene)
data <- rbind(drug_gene, sum_row)
#只保留大于等于10的靶点药物一共87个#
selected_columns <- data[, sum_row >= 10]
drug_gene1<-selected_columns[-105,]
drug<-colnames(drug_gene1)
drug<-as.data.frame(drug)
write.csv(drug,"87drugs.csv",row.names = F)#注意把NSC id加上,在excel中直接vlookup

#87种药物的DTI#
DTI<-vroom("D://subject//02 Drug prediction//dataset//02 article//data//87drug_gene_interaction.csv")
#取基因名称（行名）#
gene_name<-unique(DTI$GENE)
#取药物名称（列名）#
drug_name<-unique(DTI$Drug)
#建立空表#
drug_gene_int<-matrix(nrow=length(gene_name),ncol=length(drug_name))
row.names(drug_gene_int)<-gene_name
colnames(drug_gene_int)<-drug_name
#建立01矩阵#
for (i in 1:nrow(DTI)) {
  drug<-DTI[i,1]
  gene<-DTI[i,2]
  drug_gene_int[match(gene,rownames(drug_gene_int)),match(drug,colnames(drug_gene_int))]<-1
}

drug_gene_int[is.na(drug_gene_int)]<-0
write.csv(drug_gene_int,"87drug_gene01.csv",row.names = T)

#####以上是104种drug的，现在需要找到其中87个drug的NCI中DDC#######
#方法：删除所有不在drug87列表中的#
drug87<-colnames(selected_columns)
NCI87 <- NCI60[NCI60$drug1name %in% drug87 & NCI60$drug2name %in% drug87, ]
drug87<-colnames(selected_columns)
NCI87 <- NCI60[NCI60$drug1name %in% drug87 & NCI60$drug2name %in% drug87, ]
write.csv(NCI87,"NCI_train.csv",row.names = T)
#测试集数据，依据87种药物-60种细胞系相互组合，减去NCI87就是测试集#
drug<-vroom("87drug//87drugs.csv",col_names = F)
drug<-as.data.frame(drug)
cell<-unique(NCI87$`cell line`)
combination <- expand.grid(drug$X2,drug$X2,cell)
colnames(combination) <- c("drug1", "drug2","cell line")
#不考虑顺序问题
class(combination$drug1)
combination$drug1 <- pmin(combination$drug1, combination$drug2)
combination$drug2 <- pmax(combination$drug1, combination$drug2)

combination <- unique(combination)
# 过滤掉Drug列和Cell列不相同的情况,这些是所有的可能#
DDC <- combination[combination$drug1 != combination$drug2, ]
write.csv(DDC,"87NCI_DDC_all.csv",row.names = T)
#训练集是NCI87，测试集用上面去掉#
train<-NCI87[,c(1,2,4)]
library(dplyr)
test<- anti_join(DDC,train,by = c("drug1", "drug2", "cell line"))
write.csv(test,"NCI_test.csv",row.names = T)

#####特征处理####
#####药物特征maccs####
drugmaccs<-vroom("D://subject//02 Drug prediction//dataset//02 article//data//drug_maccs.csv")
drugmaccs<-as.data.frame(drugmaccs)
rownames(drugmaccs)<-drugmaccs[,1]
drugmaccs<-drugmaccs[,c(-1,-2)]
#data是104种药物的，选择之前筛选的87种药物数据#
drug<-drug87
drugmaccs<-drugmaccs[which(rownames(drugmaccs)%in%drug),]
sum_row <- colSums(drugmaccs)
data <- rbind(drugmaccs, sum_row)
selected_columns <- data[, sum_row > 10]
drugmaccs<-selected_columns[-88,]
#DD_maccs矩阵#这里比之前多了两列，药物名称，在excel中做的匹配
DDC<-vroom("D://subject//02 Drug prediction//dataset//02 article//data//87NCI_DDC_all.csv")
DDC<-as.data.frame(DDC[,4:5])
DD<-unique(DDC)
DD<-cbind(DD, matrix(NA, nrow = 3741, ncol = 109))
a<-colnames(drugmaccs)
colnames(DD)<-c("drugname1","drugname2",a)
for (i in 1:nrow(DD)) {
  a <-as.character(DD[i, 1])
  b <-as.character(DD[i, 2])
  DD[i,3:111]<-drugmaccs[which(rownames(drugmaccs)==a),]+ drugmaccs[which(rownames(drugmaccs)==b),]
}
write.csv(DD,"DD_maccs.csv",row.names = T)

######药物特征tox#####
drugtox<-vroom("D://subject//02 Drug prediction//dataset//02 article//data//drugtox.csv")
drugtox<-as.data.frame(drugtox)
rownames(drugtox)<-drugtox[,1]
drugtox<-drugtox[,-1]
drugtox<-t(drugtox)
#data是104种药物的，选择之前筛选的87种药物数据#
drug<-drug87
drug<-as.data.frame(drug)
drug<-drug[,1]
setdiff(drug,rownames(drugtox))
drugtox<-drugtox[which(rownames(drugtox)%in%drug),]
sum_row <- colSums(drugtox)
data <- rbind(drugtox, sum_row)
selected_columns <- data[, sum_row > 1]
drugtox<-selected_columns[-88,]
#DD_tox矩阵#
DDC<-vroom("D://subject//02 Drug prediction//dataset//02 article//data//87NCI_DDC_all.csv")
DDC<-as.data.frame(DDC[,4:5])
DD<-unique(DDC)
DD<-cbind(DD, matrix(NA, nrow = 3741, ncol = 97))
a<-colnames(drugtox)
colnames(DD)<-c("drugname1","drugname2",a)
for (i in 1:nrow(DD)) {
  a <-as.character(DD[i, 1])
  b <-as.character(DD[i, 2])
  DD[i,3:99]<-drugtox[which(rownames(drugtox)==a),]+ drugtox[which(rownames(drugtox)==b),]
}
write.csv(DD,"DD_tox.csv",row.names = T)
#####药物通路特征#######
#先制作药物-通路01矩阵#
drug_gene<-vroom("D://subject//02 Drug prediction//dataset//02 article//data//87drug_gene01.csv")
drug_gene<-as.data.frame(drug_gene)
rownames(drug_gene)<-drug_gene[,1]
drug_gene<-drug_gene[,-1]
drug<-colnames(drug_gene)
gene<-row.names(drug_gene)
#先对kegg pathway进行初筛，筛去小于X个基因的#
pathway<-vroom("D://subject//02 Drug prediction//dataset//02 article//data//KEGGREST_WithGene.csv")
pathway<-as.data.frame(pathway[-1])
#统计非na个数，用来筛选通路#
df <- data.frame()
# 遍历每一列，计算非NA值的个数
for (col_name in colnames(pathway)) {
  non_na_count <- sum(!is.na(pathway[[col_name]]))
  df <- rbind(df, non_na_count)
}
# 将结果转置，使其成为第532行
df <- t(df)
colnames(df)<-colnames(pathway)
# 将结果添加到原始数据框的末尾
pathway <- rbind(pathway, df)
#初筛#依据通路基因个数
selected_columns1 <- pathway[, df > 10]
pathway10<-selected_columns1[-532,]
#制作pathway10的基因01矩阵#
#取通路名称#
pathway_name<-colnames(pathway10)
#取全部基因并集#
gene_name<-as.matrix(pathway10)
for (i in 1:323) {
  gene<-as.matrix(pathway10[i])
  gene_name<-union(gene,gene_name)
}
#建立空表##
pathway_gene<-matrix(nrow=8152,ncol=323)
row.names(pathway_gene)<-gene_name
colnames(pathway_gene)<-pathway_name
#建立01矩阵#
for (i in 1:323) {
  kegg_gene_1<-as.matrix(pathway10[,i])
  pathway_gene[match(kegg_gene_1,rownames(pathway_gene)),i]<-1
  pathway_gene[-match(kegg_gene_1,rownames(pathway_gene)),i]<-0
}
#之前制作过340通路和104种药物的01矩阵#
#kegg通路基因#
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("GSEABase")
BiocManager::install("GSVA")
library(dplyr)
library(GSEABase)
library(GSVA)
kegg_gmt_GSVA <- readLines('D://subject//02 Drug prediction//dataset//02 article//data//KEGGREST_WithGene.gmt')
strsplit_no_name <- function(gmt.list_layer){
  as.character(unlist(strsplit(gmt.list_layer, split = '\t',
                               fixed = T)))[-2]
}
database_list_GSVA <- lapply(kegg_gmt_GSVA, strsplit_no_name)
for (layers in 1:length(database_list_GSVA)){
  names(database_list_GSVA)[layers] <- database_list_GSVA[layers][[1]][1]
  database_list_GSVA[layers][[1]] <- database_list_GSVA[layers][[1]][-1]
}
#药物对应的基因#
#先将药物01矩阵制成列表，行为药物
library(vroom)
Drug_gene<-vroom("D://subject//02 Drug prediction//dataset//02 article//data//drug_gene_01矩阵.csv")
a<-data.frame()
for (i in 1:1583) {
  for (j in 1:105) {
    if (Drug_gene[i,j]==1){a[i,j]<-Drug_gene[i,1]}
    else {a[i,j]<-""}
  }
}
colnames(a)<-colnames(Drug_gene)
a[,1]<-Drug_gene[1:1583,1]
row.names(a)<-a[,1]
Drug_gene_1<-t(a[,-1])
list <- split(Drug_gene_1, 1:nrow(Drug_gene_1))
#计算占比#之前打算依据药物通路jaccard作为单药特征，但是看了一些文献还是制作成drug-pathway01矩阵
#对表格做循环，计算方法drug与pathway求交集，再除以并集#
head(database_list_GSVA)
head(list)
drug_pathway_score1<-data.frame()
for (i in 1:340) {
  for (j in 1:104) {
    b<-length(which(as.data.frame(database_list_GSVA[i]) != "")) ###一条pathway中基因个数
    c<-length(which(as.data.frame(list[j])[,1]!= ""))  #######药物基因个数
    d<-unlist(database_list_GSVA[[i]][which(as.data.frame(database_list_GSVA[i]) != "")])  #####pathway中的基因
    e<-unlist(list[[j]][which(as.data.frame(list[j])[,1]!= "")])   ######药物中的基因
    f<-length(intersect(d,e))
    g<-length(union(d,e))
    score<-f/g
    drug_pathway_score1[j,i]<-score
    
  }
  
}
colnames(drug_pathway_score1)<-names(database_list_GSVA)
row.names(drug_pathway_score1)<-row.names(Drug_gene_1)
write.csv(drug_pathway_score1,"drug_pathway_score1.csv")
#检查有无行列全为0的情况#
##列##
del <- c()
for (i in seq(1, ncol(drug_pathway_score1))) {
  if(sum(drug_pathway_score1[,i])==0){
    print(i)
    del <- append(del, -i)
  }
}

drug_pathway_score2 <- drug_pathway_score1[,del]
write.csv(drug_pathway_score2,"drug_pathway_score2.csv")
##行##
del <- c()
for (i in 1:104){
  if(sum(drug_pathway_score1[i,])==0){
    print(i)
    
  }
}
drug_pathway_score3 <- drug_pathway_score2[-101,]
write.csv(drug_pathway_score3,"drug_pathway_score3.csv")#最终选择score3
#drug_pathway_01矩阵(340条).csv在excel中制作#选出其中87个药物以及323条通路的
drug_path<-vroom("D://subject//02 Drug prediction//dataset//02 article//data//drug_pathway_01矩阵(340条).csv")
drug87<-colnames(drug_gene)
path10<-colnames(pathway10)
drug_path<-as.data.frame(drug_path)
rownames(drug_path)<-drug_path[,1]
drug_path<-drug_path[,-1]
#drug-pathway10#
selected_rows <- rownames(drug_path) %in% drug87
selected_cols <- colnames(drug_path) %in% path10
# 使用逻辑条件选择数据
drug_path10 <- drug_path[selected_rows, selected_cols]
sum_row <- colSums(drug_path10)
data <- rbind(drug_path10, sum_row)
selected_columns <- data[, sum_row >43]#确保有50%的药物作用到这些通路上
drug_path10<-selected_columns[-88,]
write.csv(drug_path10,"drug_path10.csv",row.names = T)
#DDP矩阵制作#
#计算pathway的相关性，直接从之前340中挑选即可#
path10<-colnames(drug_path10)
pathcor<-vroom("D://subject//02 Drug prediction//dataset//02 article//data//pathway_cor_jac.csv")
pathcor<-as.data.frame(pathcor)
rownames(pathcor)<-pathcor[,1]
pathcor<-pathcor[,-1]

path10<-pathcor[row.names(pathcor) %in% path10,colnames(pathcor) %in% path10]
write.csv(path10,"pathcor10.csv",row.names = T)
#计算DDP#
#DDP10#
drug_path10<-vroom("87drug//drug_path10.csv")
drug_path10<-as.data.frame(drug_path10)
row.names(drug_path10)<-drug_path10[,1]
drug_path10<-drug_path10[,-1]
path10<-vroom("87drug//pathcor10.csv")
path10<-as.data.frame(path10)
row.names(path10)<-path10[,1]
path10<-path10[,-1]

DDC10<-vroom("87drug//87NCI_DDC_all.csv")
DDC10<-as.data.frame(DDC10[,4:5])
DD10<-unique(DDC10)
DD10<-cbind(DD10, matrix(NA, nrow = 3741, ncol = 157))
a<-colnames(path10)
colnames(DD10)<-c("drugname1","drugname2",a)
for (i in 1:3741) {
  for (j in 3:156) {
    drug1<-DD10[i,1]
    drug2<-DD10[i,2]
    pathway<-colnames(DD10[j])
    D1<-drug_path10[match(drug1,row.names(drug_path10)),]
    D2<-drug_path10[match(drug2,row.names(drug_path10)),]
    pathway1<-path10[,match(pathway,colnames(path10))]
    a<-drug_path10[drug1,pathway]
    b<-drug_path10[drug2,pathway]
    DD10[i,j]<-sum((a*D2*pathway1)+(b*D1*pathway1))/2
    
  }
  
}
write.csv(DD10,"DDP10.csv",row.names = T)










#####细胞系通路特征#####
#之前cnv、exp、RNAi、methy、mut结果就不再重做了，前四个做的gsva，mut做的jaccard#
#5种cell line特征与NCI选择细胞系#
library(vroom)
cnv<-vroom("D://subject//02 Drug prediction//dataset//feature最终版//02 特征表格//kegg_cnv_gsva.csv")
RNAi<-vroom("D://subject//02 Drug prediction//dataset//feature最终版//02 特征表格//kegg_RNAi_gsva.csv")
exp<-vroom("D://subject//02 Drug prediction//dataset//feature最终版//02 特征表格//kegg_exp_gsva.csv")
methy<-vroom("D://subject//02 Drug prediction//dataset//feature最终版//02 特征表格//kegg_methy_gsva.csv")
mut<-vroom("D://subject//02 Drug prediction//dataset//feature最终版//02 特征表格//kegg_mut_jaccard.csv")
NCI<-vroom("D://subject//02 Drug prediction//dataset//feature最终版//02 特征表格//NCI-DDC.csv")
cnv<-colnames(cnv)[-1]
RNAi<-colnames(RNAi)[-1]
exp<-colnames(exp)[-1]
methy<-colnames(methy)[-1]
mut<-colnames(mut)[-1]
NCI<-unique(NCI$DepMap_ID)
#6组数据绘制韦恩图#
library(venn)      
library(VennDiagram) 
venn_list<-list(cnv,RNAi,exp,methy,mut,NCI)
names(venn_list) <- c("cnv","RNAi","exp","methy","mut","NCI")

pdf("venn_diagram.pdf",width = 10, height = 10)
venn(venn_list,
     zcolor='style',
     opacity = 0.3,
     box = F,
     ilcs = 2,
     sncs = 3.2,
     ellipse = F,
     borders = TRUE
)
dev.off()


inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
inter <- subset(inter, select = -..values.. )
inter <- subset(inter, select = -..set.. )
write.table(inter, "6组数据交集情况.csv", row.names = FALSE, sep = ',', quote = FALSE)


#选择pathway——38细胞系#
#cnv#
cellline<-read.csv("87drug//38cell line.csv",header = FALSE)
cnv<-vroom("D://subject//02 Drug prediction//dataset//feature最终版//02 特征表格//kegg_cnv_gsva.csv")
cnv<-as.data.frame(cnv)
rownames(cnv)<-cnv[,1]
cnv<-cnv[,-1]
pathway_name10<-colnames(path10)
cnv10<-cnv[which(rownames(cnv)%in%pathway_name10),colnames(cnv)%in%cellline$V1]
write.csv(cnv10,"cnv10.csv",row.names = T)

#RNAi#
RNAi<-vroom("D://subject//02 Drug prediction//dataset//feature最终版//02 特征表格//kegg_RNAi_gsva.csv")
RNAi<-as.data.frame(RNAi)
rownames(RNAi)<-RNAi[,1]
RNAi<-RNAi[,-1]
RNAi10<-RNAi[which(rownames(RNAi)%in%pathway_name10),colnames(RNAi)%in%cellline$V1]
write.csv(RNAi10,"RNAi10.csv",row.names = T)


#exp#
exp<-vroom("D://subject//02 Drug prediction//dataset//feature最终版//02 特征表格//kegg_exp_gsva.csv")
exp<-as.data.frame(exp)
rownames(exp)<-exp[,1]
exp<-exp[,-1]
exp10<-exp[which(rownames(exp)%in%pathway_name10),colnames(exp)%in%cellline$V1]
write.csv(exp10,"exp10.csv",row.names = T)

#methy#
methy<-vroom("D://subject//02 Drug prediction//dataset//feature最终版//02 特征表格//kegg_methy_gsva.csv")
methy<-as.data.frame(methy)
rownames(methy)<-methy[,1]
methy<-methy[,-1]
methy10<-methy[which(rownames(methy)%in%pathway_name10),colnames(methy)%in%cellline$V1]
write.csv(methy10,"methy10.csv",row.names = T)

#mut#
mut<-vroom("D://subject//02 Drug prediction//dataset//feature最终版//02 特征表格//kegg_mut_jaccard.csv")
mut<-as.data.frame(mut)
rownames(mut)<-mut[,1]
mut<-mut[,-1]
mut10<-mut[which(rownames(mut)%in%pathway_name10),colnames(mut)%in%cellline$V1]
write.csv(mut10,"mut10.csv",row.names = T)











#####细胞系筛选60仅剩38，这38个具有5种cell line特征########
#后面只使用了38个细胞系，对train和test数据进行筛选#
train<-vroom("87drug//NCI_train.csv")
train<-as.data.frame(train)
train<-train[,-1]
cellline
train<-train[which(train$depmap_id%in%cellline$V1),]
write.csv(train,"NCI_train38.csv",row.names = T)

test<-vroom("87drug//NCI_test.csv")
test<-as.data.frame(test)
test<-test[,-1]

test<-test[which(test$Depmap_ID%in%cellline$V1),]
write.csv(test,"NCI_test38.csv",row.names = T)




#后面只使用了38个细胞系，对train和test数据进行筛选#
train<-vroom("87drug//NCI_train.csv")
train<-as.data.frame(train)
train<-train[,-1]
cellline
train<-train[which(train$depmap_id%in%cellline$V1),]
write.csv(train,"NCI_train38.csv",row.names = T)

test<-vroom("87drug//NCI_test.csv")
test<-as.data.frame(test)
test<-test[,-1]

test<-test[which(test$Depmap_ID%in%cellline$V1),]
write.csv(test,"NCI_test38.csv",row.names = T)


