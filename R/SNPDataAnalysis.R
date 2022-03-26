setwd('E:/NWPU/third/trypathway/')
snpstarter <- read.table(file="snp_GSE28127.txt", header=TRUE, sep="\t", row.names=1,stringsAsFactors = FALSE)

snpMatrix<-as.matrix(snpstarter)
#释放内存
rm(snpstarter)
gc()

rowname<-row.names(snpMatrix)

locs <-strsplit(rowname, "-")
rowname_snp<-sapply(locs, "[", 1)

tail(rowname_snp)

snpTemp<-chartr("AB","01",snpMatrix)
#释放内存
rm(snpMatrix)
gc()

snpNum<-apply(snpTemp,2,as.numeric)
#释放内存
rm(snpTemp)
gc()

snpNum[snpNum==11]<-2
row.names(snpNum)<-rowname_snp

#install.packages("DMwR2")
library(DMwR2)
#列出所有缺失值小于20%（0.2）的样本号
snpOutNA <- snpNum[-manyNAs(snpNum,0.2),]
#释放内存
rm(snpNum)
gc()

#imputation
for(i in 1:nrow(snpOutNA)) {
  uniqv <- unique(snpOutNA[i,])
  snpOutNA[i, is.na(snpOutNA[i,])] <- uniqv[which.max(tabulate(match(snpOutNA[i,], uniqv)))]
}

## cal MAF
n0 <- apply(snpOutNA==0,1,sum,na.rm=T)
n1 <- apply(snpOutNA==1,1,sum,na.rm=T)
n2 <- apply(snpOutNA==2,1,sum,na.rm=T)

n <- n0 + n1 + n2

## calculate allele frequencies
p <- ((2*n0)+n1)/(2*n)
q <- 1 - p
maf <- pmin(p, q)

#释放内存
gc()

snpoutmaf<-snpOutNA[maf>0.05,]
#释放内存
rm(snpOutNA)
gc()

write.table( snpoutmaf, "SNP_first.txt", row.names=TRUE, col.names=TRUE,sep="\t",quote = FALSE)

