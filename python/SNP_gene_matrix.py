import pandas as pd
import csv

def list(l1, l2):
    for i in range(0, len(l1)):
        if l1[i] not in l2:
            l2.append(l1[i])
    return l2

def list_with_repeat(l1, l2):
    for i in range(0, len(l1)):
        l2.append(l1[i])
    return l2

def read(fname):
    data = pd.read_csv(fname,sep = "\t",header = 0)
    print("data over")
    SNP = []
    gene = []
    SNP_name = []
    gene_name = []
    SNP_list = list_with_repeat(data["RefSNP_id"], SNP)
    gene_list = list_with_repeat(data["SYMBOL"], gene)
    SNP_name = list(data["RefSNP_id"], SNP_name)
    gene_name = list(data["SYMBOL"], gene_name)
    SNP_name.sort()
    gene_name.sort()
    print("read over")
    return SNP_list, gene_list,SNP_name,gene_name

def toMatrix(SNP_list, gene_list,SNP_name,gene_name):
    matrix = [[0 for _ in range(len(SNP_name))] for _ in range(len(gene_name))]
    print("have matrix")
    for i in range(0,len(SNP_list)):
        x = gene_list[i]
        y = SNP_list[i]
        x_index = gene_name.index(x)
        y_index = SNP_name.index(y)
        if matrix[x_index][y_index]==0:
            matrix[x_index][y_index] += 1
            print(x_index,y_index)
        else:
            pass
    return matrix

def save_matrix_name(SNP_name,gene_name,savenamepath):
    with open(savenamepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(SNP_name)
        writer.writerow(gene_name)

def save_matrix(matrix,savepath):
    with open(savepath, 'w', newline='') as ff:
        writer = csv.writer(ff)
        #write row
        for j in range(0,len(gene_name)):
            writer. writerow(matrix[j])

if __name__ == '__main__':
    openpath = 'E:/NWPU/third/trypath/snps_gene_eQTL.txt'
    SNP_list, gene_list,SNP_name,gene_name = read(openpath)
    matrix = toMatrix(SNP_list, gene_list,SNP_name,gene_name)

    savenamepath = 'E:/NWPU/third/trypath/row_col_name.csv'
    save_matrix_name(SNP_name, gene_name, savenamepath)

    savepath = 'E:/NWPU/third/trypath/snp_gene_matrix.csv'
    save_matrix(matrix,savepath)