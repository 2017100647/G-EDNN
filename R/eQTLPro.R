setwd('E:/NWPU/third/trypathway/')

#Get eQTL
#只取GTex与snp\gene相关的两列
# file.remove("eQTLdata_liver.txt")
# con <- file("Liver.allpairs.txt", "r")
# line=readLines(con,n=1)
# i=0
# while( length(line) != 0 ) {
#   linetemp<-strsplit(line,split = "\t")
#   linedata<-linetemp[[1]][c(1,2,6,7)]
#   if(linetemp[[1]][6]>0.05 && linetemp[[1]][7]<0.1)
#   {
#    cat(linedata, file = "eQTLdata_liver.txt", sep = "\t", append = TRUE)
#    cat("", file = "eQTLdata_liver.txt", sep = "\n", append = TRUE)
#   }
#   line=readLines(con,n=1)
#   i=i+1
#   print(i)
# }
# close(con)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPlocs.Hsapiens.dbSNP142.GRCh37")
gc()
# 用snp_sample相关的GTex数据
snp_sample_raw <- read.table(file="SNP_first.txt", header=TRUE, sep="\t",stringsAsFactors = FALSE)

library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
snp_ref <- SNPlocs.Hsapiens.dbSNP142.GRCh37
snp_sample_pos<-snpsById(snp_ref, rownames(snp_sample_raw), ifnotfound="drop")
seqlevelsStyle(snp_sample_pos) <- "NCBI"
snp_sample_pos_dataframe<-as.data.frame(snp_sample_pos)
snp_sample_pos_dataframe_sel<-snp_sample_pos_dataframe[c("RefSNP_id","seqnames","pos")]
snp_sample_pos_dataframe_sel$seqnames<-sub("ch","",snp_sample_pos_dataframe_sel$seqnames)
tail(snp_sample_pos_dataframe_sel)
write.table( snp_sample_pos_dataframe_sel, "snp_sample_pos_dataframe_sel.txt", row.names=FALSE, col.names=TRUE, sep="\t",quote = FALSE)
#SNP name chr pos
snp_sample_pos_dataframe_sel <- read.table(file="snp_sample_pos_dataframe_sel.txt", header=TRUE, sep="\t")

#2、gtex_Blood转换GRanges class
gtex_Blood <- read.table(file="eQTLdata_liver.txt", header=TRUE, sep="\t")
locs <-strsplit(gtex_Blood$variant_id, "_")
gtex_Blood$seqnames<-sapply(locs, "[", 1)
gtex_Blood$pos <-sapply(locs, "[", 2)
tail(gtex_Blood)
#限定条件
gtex_Blood1<-gtex_Blood[which(gtex_Blood$pval_nominal<0.05),]
gtex_Blood<-gtex_Blood1[which(gtex_Blood1$maf>0.05),]

snps_sample_in_gtex<- merge(snp_sample_pos_dataframe_sel, gtex_Blood, by=c("seqnames","pos"),all= FALSE)

snps_sample_in_gtex$ENSEMBL<- sub("(ENSG[0-9]+)\\.[0-9]+", "\\1",snps_sample_in_gtex$gene_id)

snps_eQTL<-snps_sample_in_gtex[c("ENSEMBL","RefSNP_id")]

#Get Pos associated

library(biomaRt) 
snp_sample_raw <- read.table(file="SNP_first.txt", header=TRUE, sep="\t",stringsAsFactors = FALSE)

#snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp") 
snp_mart = useEnsembl("ENSEMBL_MART_SNP", dataset="hsapiens_snp",mirror = "asia") 

snp_attributes = c("refsnp_id", "chr_name", "chrom_start", "associated_gene", "ensembl_gene_stable_id", "minor_allele_freq") 

getENSG <- function(rs, mart = snp_mart) { 
  results <- getBM(attributes = snp_attributes,filters = "snp_filter", values = rs, mart = mart) 
  return(results) 
} 

SNP_Pos<-getENSG(rownames(snp_sample_raw))
names(SNP_Pos)[names(SNP_Pos) == "refsnp_id"] <- "RefSNP_id"
names(SNP_Pos)[names(SNP_Pos) == "ensembl_gene_stable_id"] <- "ENSEMBL"
SNP_Pos<-SNP_Pos[c("ENSEMBL","RefSNP_id")]

SNP_ENSG<-rbind(snps_eQTL, SNP_Pos) 

library(org.Hs.eg.db)
library("clusterProfiler")
# 采用bitr()函数进行转换
gene_symbol <- bitr(SNP_ENSG$ENSEMBL, fromType="ENSEMBL", toType= "SYMBOL", OrgDb=org.Hs.eg.db)
# 查看转换的结果
head(gene_symbol)

#删除等于NA的列
gene_symbol =na.omit(gene_symbol)
ENSG_symbol=gene_symbol[!duplicated(gene_symbol$SYMBOL),]#!duplicated去掉重复，从而把V1中所有的元素变成唯一

SNP_gene<- merge(SNP_ENSG, ENSG_symbol, by="ENSEMBL",all = FALSE)

SNP_gene<-na.omit(SNP_gene)
tail(SNP_gene)

write.table( SNP_gene, "SNP_gene.txt", row.names=FALSE, col.names=TRUE, sep="\t",quote = FALSE)


PPI <- read.table(file="PPI.txt", header=FALSE, sep=" ")

PPI<-PPI[which(PPI[,3]>0.3),]

Protein2Symbol<-read.table(file="protein_to_gene.txt", header=TRUE, sep="\t")

Protein2Symbol_pro<-Protein2Symbol[which(Protein2Symbol$Gene%in%SNP_gene$SYMBOL),]
#两列基因与蛋白质都要匹配上
PPIfilter<-PPI[which((PPI[,1]%in%Protein2Symbol_pro[,1]) & (PPI[,2]%in%Protein2Symbol_pro[,1])),]

PPI2Gene<-apply(PPIfilter, 2, function(x) Protein2Symbol_pro$Gene[match(x, Protein2Symbol_pro$Protein)])

#第一列与第二列交换，去重复，这样就可以看出第一列与第二列交换有没有一样的
PPI2Gene<-data.frame(t(apply(PPI2Gene,1,function(x)sort(x))))

#先去除一次重复,看本身有没有重复
PPI2Gene<-PPI2Gene[!duplicated(PPI2Gene),]

Gene2Gene<-PPI2Gene[as.character(PPI2Gene$X1) != as.character(PPI2Gene$X2),]

write.table( Gene2Gene, "Gene2Gene.txt", row.names=FALSE, col.names=TRUE, sep="\t",quote = FALSE)
# 用gene相关的GTex数据
Gene2Gene <- read.table(file="Gene2Gene.txt",header=TRUE, sep="\t",stringsAsFactors = FALSE)

Gene_sel<-union(as.character(Gene2Gene$X1),as.character(Gene2Gene$X2))

gene_rowname<-as.data.frame(Gene_sel)

names(gene_rowname)<- "SYMBOL"
snps_gene_in_gtex<- merge(SNP_gene, gene_rowname, by="SYMBOL",all = FALSE)
write.table( snps_gene_in_gtex, "snps_gene_eQTL.txt", row.names=FALSE, col.names=TRUE, sep="\t",quote = FALSE)

SNPname_Last<-intersect(rownames(snp_sample_raw),snps_gene_in_gtex$RefSNP_id)

SNPpro_last<-snp_sample_raw[rownames(snp_sample_raw) %in% SNPname_Last,]

sampleR = t(SNPpro_last)
write.table( sampleR, "SNP_second.txt", row.names=TRUE, col.names=TRUE, sep="\t",quote = FALSE)
