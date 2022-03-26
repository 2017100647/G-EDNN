import pandas as pd
import csv
import numpy as np

def list(l1, l2):
    for i in range(0, len(l1)):
        if l1[i] not in l2:
            l2.append(l1[i])
    return l2

def read(fname):
    gene = []
    data = []
    with open(fname, 'r') as csvfile:
        reader = csv.reader(csvfile,delimiter = "\t")
        for row in reader:
            gene = list(row, gene)
    with open(fname, 'r') as csvfile1:
        reader1 = csv.reader(csvfile1,delimiter = "\t")
        data = [row1 for row1 in reader1]

    print("data over")

    gene.sort()

    print("read over")
    return gene,data

def toMatrix(gene,data):
    matrix = [[0 for _ in range(len(gene))] for _ in range(len(gene))]
    print("have matrix")

    for i in range(0,len(data)):
        for j in range(0,len(data[i])):
            x = data[i][j]
            for m in range(0,len(data[i])):
                y = data[i][m]
                x_index = gene.index(x)
                y_index = gene.index(y)
                if matrix[x_index][y_index]==0:
                    matrix[x_index][y_index] += 1
                    print(x_index,y_index)
                else:
                    pass
    print([*zip(*matrix)])
    return matrix

def save_matrix_name(gene,savenamepath):
    with open(savenamepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(gene)

def save_matrix(matrix,savepath):
    with open(savepath, 'w', newline='') as ff:
        writer = csv.writer(ff)
        for j in range(0,len(gene)):
            writer. writerow(matrix[j])


if __name__ == '__main__':
    #openpath = 'E:/NWPU/third/data_GSE95496/Gene2Gene.txt'  E:\NWPU\third\trypathway
    openpath = 'E:/NWPU/third/trypath/Gene2Gene.txt'
    #openpath = 'E:/NWPU/third/trypath/Gene2Gene.txt'
    gene,data = read(openpath)
    matrix = toMatrix(gene,data)

    #savenamepath = 'E:/NWPU/third/data_GSE95496/gene_name.csv'
    #savenamepath = 'E:/NWPU/third/trypath/gene_name.csv'
    savenamepath = 'E:/NWPU/third/trypath/gene_name.csv'
    save_matrix_name(gene, savenamepath)

    #savepath = 'E:/NWPU/third/data_GSE95496/gene_gene_matrix.csv'
    #savepath = 'E:/NWPU/third/trypath/gene_gene_matrix.csv'
    savepath = 'E:/NWPU/third/trypath/gene_gene_matrix.csv'
    matrix_T = np.transpose(matrix)
    matrix_last = matrix_T+matrix
    matrix_last[matrix_last>1]=1
    save_matrix(matrix_last,savepath)
