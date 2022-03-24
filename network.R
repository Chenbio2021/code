install.packages("BiocManager")      ######为使用下载功能，安装相关的包
BiocManager::install("TCGAbiolinks")
install.packages("SummarizedExperiment")
library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)
TCGAbiolinks:::getGDCprojects()$project_id
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")
#乳腺癌、
getGDCprojects()$project_id
CancerProject <- "TCGA-BRCA"
# 设置数据存放路径
DataDirectory <- paste0("./",gsub("-","_",CancerProject)) # paste0：字符串连接，sep=“”, gsub(“-","_",CancerProject)  替换"-"为"_"
# 数据名称
FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")
# 下载数据（TCGA中与乳腺癌相关的基因表达原始数据）
query <- GDCquery(project = CancerProject,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

GDCdownload(query = query, directory = DataDirectory)


###表达谱处理###
# 定义样本列表
listSamples <-
  c("TCGA-E9-A1NG-11A-52R-A14M-07", "TCGA-BH-A1FC-11A-32R-A13Q-07", "TCGA-A7-A13G-11A-51R-A13Q-07",
    "TCGA-BH-A0DK-11A-13R-A089-07", "TCGA-E9-A1RH-11A-34R-A169-07", "TCGA-BH-A0AU-01A-11R-A12P-07",
    "TCGA-C8-A1HJ-01A-11R-A13Q-07", "TCGA-A7-A13D-01A-13R-A12P-07", "TCGA-A2-A0CV-01A-31R-A115-07",
    "TCGA-AQ-A0Y5-01A-11R-A14M-07")
# 查询 Illumina HiSeq 平台的数据
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  experimental.strategy = "RNA-Seq",
  platform = "Illumina HiSeq",
  file.type = "results",
  barcode = listSamples,
  legacy = TRUE
)
# 下载 对应样本的信息
GDCdownload(query)
# 预处理
BRCARnaseqSE <- GDCprepare(query)
#处理成行为 geneID  列为样本 (barcode) 的矩阵
BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")

#行名是基因 symbol 与 ID 合并的形式，可以与 rowRanges(BRCARnaseqSE) 获取的基因信息来进行转换
feature <- rowRanges(BRCARnaseqSE)
rownames(BRCAMatrix) <- feature[rownames(BRCAMatrix), "gene_id"]$gene_id

#差异表达分析
#默认使用的 edgeR，要使用limma包时，设置pipeline=“limma”
# 标准化
dataNorm <- TCGAanalyze_Normalization(tabDF = BRCAMatrix, geneInfo =  TCGAbiolinks::geneInfo)

# 使用分位数来过滤基因
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm, method = "quantile", qnt.cut =  0.25
)

# 挑选正常样本：NT
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = c("NT"))
# 挑选肿瘤样本：TP
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = c("TP"))

# 差异表达分析
dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[, samplesNT],
  mat2 = dataFilt[, samplesTP],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT")
# 获取差异表达基因的表达水平
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
  dataDEGs, "Tumor", "Normal",
  dataFilt[, samplesTP], dataFilt[, samplesNT]
)

#HTSeq数据
#获取数据
CancerProject <- "TCGA-BRCA"
query <- GDCquery(
  project = CancerProject,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts")

samplesDown <- getResults(query,cols=c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "TP")
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "NT")
dataSmTP_short <- dataSmTP[1:10]
dataSmNT_short <- dataSmNT[1:10]
#储存路径
DataDirectory <- paste0("./",gsub("-","_",CancerProject))
# 数据名称
FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")
queryDown <- GDCquery(
  project = CancerProject,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts",
  barcode = c(dataSmTP_short, dataSmNT_short)
)
# 下载并保存到指定路径
GDCdownload(query = queryDown, directory = DataDirectory)
# 预处理数据并保存到本地
dataPrep <- GDCprepare(query = queryDown, save = TRUE, directory =  DataDirectory,
                       save.filename = FileNameData)

# 获取 count 数据
dataPrep <- TCGAanalyze_Preprocessing(
  object = dataPrep, cor.cut = 0.6, datatype = "HTSeq - Counts")


# 对数据进行标准化
dataNorm <- TCGAanalyze_Normalization(
  tabDF = dataPrep, geneInfo = TCGAbiolinks::geneInfoHT, method = "gcContent"
)
#比较标准化前后数据分布的差别
boxplot(dataPrep, outline = FALSE, names = FALSE)
boxplot(dataNorm, outline = FALSE, names = FALSE)
#数据过滤和差异表达分析
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm, method = "quantile",  qnt.cut =  0.25
)   
dim(dataFilt)
dim(dataNorm)
dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[, dataSmTP_short],
  mat2 = dataFilt[, dataSmNT_short],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT"
)

###miRNA
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Gene expression",
  data.type = "miRNA gene quantification",
  legacy = TRUE)

head(getResults(query)$file_name)
CancerProject <- "TCGA-BRCA"
DataDirectory <- paste0("./", gsub("-", "_", CancerProject))
FileNameData <-paste0(DataDirectory, "_", "miRNA_gene_quantification", ".rda")
# 查询对应样本
query.miR <- GDCquery(
  project = CancerProject,
  data.category = "Gene expression",
  data.type = "miRNA gene quantification",
  file.type = "mirna",
  legacy = TRUE)
samplesDown.miR <- getResults(query.miR, cols = c("cases"))

dataSmTP.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "TP")
dataSmNT.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR, typesample = "NT")
dataSmTP_short.miR <- dataSmTP.miR[1:10]
dataSmNT_short.miR <- dataSmNT.miR[1:10]
#下载并预处理数据
queryDown.miR <- GDCquery(
  project = CancerProject,
  data.category = "Gene expression",
  data.type = "miRNA gene quantification",
  file.type = "mirna",
  legacy = TRUE,
  barcode = c(dataSmTP_short.miR, dataSmNT_short.miR))

GDCdownload(query = queryDown.miR, directory = DataDirectory)

dataAssy.miR <- GDCprepare(
  query = queryDown.miR,
  save = TRUE,
  save.filename = FileNameData,
  summarizedExperiment = TRUE,
  directory = DataDirectory)
#提取出 read_count 数据，并进行差异表达分析
library(dplyr)
library(tibble)
dataAssy.miR %<>%
  column_to_rownames(var = "miRNA_ID") %>%
  dplyr::select(starts_with("read_count_")) %>%
  rename_all(function(x) gsub("read_count_", "", x))

dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataAssy.miR, method = "quantile",
  qnt.cut =  0.25)


#GO富集分析
Genelist <- rownames(dataDEGsFiltLevel)
ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist)
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP), 
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = Genelist,
                        filename="a.pdf",
                        nBar = 10)

#生存分析
clin.gbm <- GDCquery_clinic("TCGA-BRCA", "clinical")
TCGAanalyze_survival(clin.gbm,"gender",main = "TCGA Set\n BRCA",
                     height = 10, width=10,
                     filename = "./survival.pdf")

#差异甲基化
# 查询 TCGA-BRCA 项目中既包含甲基化又有表达的样本
samples <- matchedMetExp("TCGA-BRCA", n = 10)
# 查询下载
query <- GDCquery(
  project = c("TCGA-BRCA"),
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 450",
  legacy = TRUE,
  barcode = samples,
  
)
GDCdownload(query)
met <- GDCprepare(query, save = FALSE)

# 删除 NA 行
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
#癌症与正常样本之间平均甲基化水平
TCGAvisualize_meanMethylation(met, groupCol = "sample_type", filename = NULL)

met_res <- TCGAanalyze_DMC(
  met,
  # 使用 sample_type 列来分组
  groupCol = "sample_type",
  # 分组名称
  group1 = "Primary Tumor",
  group2 = "Solid Tissue Normal",
  p.cut = 0.05,
  diffmean.cut = 0.15,
  save = FALSE,
  legend = "State",
  plot.filename = "./BRCA.png",
   
)


###########stringdb包###################
BiocManager::install("STRINGdb")
library(STRINGdb)
## 首先选定物种及数据库的版本，人类NCBI分类法标识符9606#######
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=200, input_directory="")
STRINGdb$methods()
STRINGdb$help("plot_network")  
##内置数据集
data(diff_exp_example1)
head(diff_exp_example1)
example1_mapped <- string_db$map( diff_exp_example1,
                                  "gene", removeUnmappedRows = TRUE )   
hits <- example1_mapped$STRING_id[1:200]
string_db$plot_network( hits )
##按p值进行过滤，并添加一个颜色列
example1_mapped_pval05 <- string_db$add_diff_exp_color( subset(example1_mapped, pvalue<0.05), logFcColStr="logFC" )
#将有效载荷信息发布到字符串服务器
payload_id <- string_db$post_payload( example1_mapped_pval05$STRING_id,colors=example1_mapped_pval05$color )
string_db$plot_network( hits, payload_id=payload_id )

#富集
enrichment <- string_db$get_enrichment( hits )
head(enrichment, n=20)
#########设置背景########
backgroundV <- example1_mapped$STRING_id[1:2000]
string_db$set_background(backgroundV)

#获取注释集
annotations <- string_db$get_annotations( hits )
head(annotations, n=20)

#聚类
clustersList <- string_db$get_clusters(example1_mapped$STRING_id[1:600])
# plot first 4 clusters
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList[[i]])
}

######获得蛋白质信息，并可以查询特定蛋白质
string_proteins <- string_db$get_proteins() 
#以tp53和atm为例，先获取其标识符
tp53 = string_db$mp( "tp53" )
atm = string_db$mp( "atm" )
#查找与目标蛋白质之间存在相互作用的蛋白质
string_db$get_neighbors( c(tp53, atm) )
#查找蛋白质之间是否存在互作
string_db$get_interactions( c(tp53, atm) )
#获取同源蛋白
string_db$get_paralogs(tp53)
#获取来自不同物种的蛋白质之间的最佳同系物。
string_db$get_homologs_besthits(tp53)

