setwd('E:/NWPU/third/trypathway/')

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install('GEOquery')

library(GEOquery)
library(Biobase)
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*6)
options( 'download.file.method.GEOquery' = 'libcurl' ) 
gset <- getGEO('GSE28127',destdir = ".",
               AnnotGPL = F,
               getGPL = F)
save(gset,file = 'GSE28127.gset.Rdata')

#取第一个元素
ob=gset[[1]]
#得到其表达矩阵
exprSet=exprs(ob)
##ob
#查看其样本名字
samples=sampleNames(ob)
pdata=pData(ob)
group_list=as.character(pdata[,2])
dim(exprSet)
#查看该表达矩阵的前几行
exprSet[1:5,1:5]

write.table( exprSet, "snp_GSE28127.txt", row.names=TRUE, col.names=TRUE, sep="\t",quote = FALSE)
