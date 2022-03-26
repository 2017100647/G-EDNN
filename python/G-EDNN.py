import tensorflow as tf
import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.utils import shuffle
import time

SNP_Path = "E:\\NWPU\\third\\trypathway\\SNP_second.txt"
SNP_gene_Path = "E:\\NWPU\\third\\trypathway\\snp_gene_matrix.csv"
gene_gene_Path =  "E:\\NWPU\\third\\trypathway\\gene_gene_matrix.csv"
sample_Path = "E:\\NWPU\\third\\trypathway\\samples.txt"


out_file = "E:\\NWPU\\third\\trypathway\\out0.1.csv"
out_file1 = "E:\\NWPU\\third\\trypathway\\ytest0.1.csv"
out_file2 = "E:\\NWPU\\third\\trypathway\\ys0.1.csv"

# SNP_Path = "E:\\NWPU\\third\\trypath\\SNP_second.txt"
# SNP_gene_Path = "E:\\NWPU\\third\\trypath\\snp_gene_matrix.csv"
# gene_gene_Path =  "E:\\NWPU\\third\\trypath\\gene_gene_matrix.csv"
# sample_Path = "E:\\NWPU\\third\\trypath\\samples.txt"
#
#
# out_file = "E:\\NWPU\\third\\trypath\\out.csv"

def readcsv(path):
    data = []
    with open(path, 'r', encoding='gbk', errors='ignore') as f:
        for line in f:
            data.append(line.split(','))

    data = np.array(data,dtype='int32')
    return data

def readtxt(path):
    data = []
    with open(path, 'r', encoding='gbk', errors='ignore') as f:
        for line in f:
            data.append(line.split('\t'))

    data = np.array(data)
    return data

## load in data
partition0 = readcsv(SNP_gene_Path)
partition0 = np.transpose(partition0)
partition1 = readcsv(gene_gene_Path)

# #read data
SNPdata = pd.read_table(SNP_Path,sep='\t',header=0)
SNP_raw_data = np.array(SNPdata, dtype='int32')
sample = readtxt(sample_Path)
label_vec = np.array(sample[:,-1], dtype=int)

# # read test data
# SNP_raw_data = readcsv(SNP_Path)
# SNP_data=np.array(SNP_raw_data)
# m,n=np.shape(SNP_data)
# sample = readcsv(sample_Path)
# label_vec = np.array(sample[:,-1], dtype=int)

labels = []
for l in label_vec:
    if l == 1:
        labels.append([0,1])
    else:
        labels.append([1,0])
labels = np.array(labels,dtype=int)

#归一化数据
# Normal_SNP=np.zeros((m,n))
# nummin=np.min(SNP_raw_data)
# nummax=np.max(SNP_raw_data)
# dif=nummax-nummin
# for i in range(m):
#     for j in range(n):
#         Normal_SNP[i][j]=(SNP_data[i][j]-nummin)/dif

## train/test data split,取0.8行
cut = int(0.1*np.shape(SNP_raw_data)[0])
SNP_raw_data, labels = shuffle(SNP_raw_data, labels)
x_train = SNP_raw_data[:cut, :]
x_test = SNP_raw_data[cut:, :]
y_train = labels[:cut, :]
y_test = labels[cut:, :]

start =  time.process_time()

#学习率
learning_rate = 0.0001
#数据送入多少次
training_epochs = 200
batch_size = 8
#精度
error=0.01
#加入一个正则项来防止过拟合
L2 = True
#显示
display_step = 1
#第一层、第二层是否dropout
droph0 = True
droph1 = True
droph2 = True
droph3 = True
droph4 = True

n_hidden_0 = np.shape(partition0)[0]
n_hidden_1 = np.shape(partition1)[0]
n_hidden_2 = 256
n_hidden_3 = 64
n_hidden_4 = 16
n_classes = 2
n_features0 = np.shape(SNP_raw_data)[1]
n_features1 = np.shape(partition1)[1]

## the constant limit for feature selection
gamma_c = 20
gamma_numerator = np.sum(partition1, axis=0)
gamma_denominator = np.sum(partition1, axis=0)
gamma_numerator[np.where(gamma_numerator>gamma_c)] = gamma_c

## initiate training logs
loss_rec = np.zeros([training_epochs, 1])
training_eval = np.zeros([training_epochs, 2])

#初始化权值和阈值
weights = {
    'h0': tf.Variable(tf.random.truncated_normal(shape=[n_features0, n_features1], stddev=0.1)),
    'h1': tf.Variable(tf.random.truncated_normal(shape=[n_features1, n_hidden_1], stddev=0.1)),
    'h2': tf.Variable(tf.random.truncated_normal(shape=[n_hidden_1, n_hidden_2], stddev=0.1)),
    'h3': tf.Variable(tf.random.truncated_normal(shape=[n_hidden_2, n_hidden_3], stddev=0.1)),
    'h4': tf.Variable(tf.random.truncated_normal(shape=[n_hidden_3, n_hidden_4], stddev=0.1)),
    'out': tf.Variable(tf.random.truncated_normal(shape=[n_hidden_4, n_classes], stddev=0.1))
}

biases = {
    'b0': tf.Variable(tf.zeros([n_features1])),
    'b1': tf.Variable(tf.zeros([n_hidden_1])),
    'b2': tf.Variable(tf.zeros([n_hidden_2])),
    'b3': tf.Variable(tf.zeros([n_hidden_3])),
    'b4': tf.Variable(tf.zeros([n_hidden_4])),
    'out': tf.Variable(tf.zeros([n_classes]))
}

#输入占位
x = tf.placeholder(tf.float32, [None, n_features0])
y = tf.placeholder(tf.int32, [None, n_classes])

keep_prob = tf.placeholder(tf.float32)
lr = tf.placeholder(tf.float32)

#构建图：前向传播
def multilayer_perceptron(x, weights, biases, keep_prob):
    #layer_0= tf.add(tf.matmul(x, tf.multiply(weights['h0'], partition0)), biases['b0'])
    layer_0= tf.add(tf.matmul(x, weights['h0']), biases['b0'])
    layer_0 = tf.nn.relu(layer_0)
    if droph0:
        layer_0 = tf.nn.dropout(layer_0, keep_prob=keep_prob)

    layer_1 = tf.add(tf.matmul(layer_0, tf.multiply(weights['h1'], partition1)), biases['b1'])
    #layer_1 = tf.add(tf.matmul(layer_0, weights['h1']), biases['b1'])
    layer_1 = tf.nn.relu(layer_1)
    if droph1:
        layer_1 = tf.nn.dropout(layer_1, keep_prob=keep_prob)

    layer_2 = tf.add(tf.matmul(layer_1, weights['h2']), biases['b2'])
    layer_2 = tf.nn.relu(layer_2)
    if droph2:
        layer_2 = tf.nn.dropout(layer_2, keep_prob)

    layer_3 = tf.add(tf.matmul(layer_2, weights['h3']), biases['b3'])
    layer_3 = tf.nn.relu(layer_3)
    if droph3:
        layer_3 = tf.nn.dropout(layer_3, keep_prob)

    layer_4 = tf.add(tf.matmul(layer_3, weights['h4']), biases['b4'])
    layer_4 = tf.nn.relu(layer_4)
    if droph4:
        layer_4 = tf.nn.dropout(layer_4, keep_prob)

    out_layer = tf.matmul(layer_4, weights['out']) + biases['out']
    return out_layer

# Construct model
pred = multilayer_perceptron(x, weights, biases, keep_prob)

# Define loss and optimizer
cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred, labels=y))
if L2:
    reg = tf.nn.l2_loss(weights['h0']) +tf.nn.l2_loss(weights['h1']) + tf.nn.l2_loss(weights['h2']) + tf.nn.l2_loss(weights['h3']) + tf.nn.l2_loss(weights['out'])
    cost = tf.reduce_mean(cost + 0.0025 * reg)
#Adam算法
optimizer = tf.train.AdamOptimizer(learning_rate=lr).minimize(cost)
#optimizer = tf.train.GradientDescentOptimizer(learning_rate=lr).minimize(cost)#梯度下降法

## Evaluation
correct_prediction = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))
y_score = tf.nn.softmax(logits=pred)

var_left = tf.reduce_sum(tf.abs(tf.multiply(weights['h1'], 256)), 0)#partition1)), 0)
var_right = tf.reduce_sum(tf.abs(weights['h2']), 1)
var_importance = tf.add(tf.multiply(tf.multiply(var_left, gamma_numerator), 1./gamma_denominator), var_right)

config = tf.ConfigProto(allow_soft_placement=True)
config.gpu_options.allocator_type = 'BFC'
#config.gpu_options.per_process_gpu_memory_fraction = 0.95
config.gpu_options.allow_growth = True
session = tf.Session(config= config)

with session as sess:

    sess.run(tf.global_variables_initializer())
    total_batch = int(np.shape(x_train)[0] / batch_size)

    ## Training cycle
    for epoch in range(training_epochs):
        avg_cost = 0.
        x_tmp, y_tmp = shuffle(x_train, y_train)
        # Loop over all batches
        for i in range(total_batch-1):
            batch_x, batch_y = x_tmp[i*batch_size:i*batch_size+batch_size], \
                                y_tmp[i*batch_size:i*batch_size+batch_size]

            _, c= sess.run([optimizer, cost], feed_dict={x: batch_x, y: batch_y, keep_prob: 0.8,
                                                        lr: learning_rate
                                                        })
            # Compute average loss
            avg_cost += c / total_batch

        del x_tmp
        del y_tmp

        ## Display logs per epoch step
        if epoch % display_step == 0:
            loss_rec[epoch] = avg_cost
            acc, y_s = sess.run([accuracy, y_score], feed_dict={x: x_train, y: y_train, keep_prob: 1})
            auc = metrics.roc_auc_score(y_train, y_s)
            training_eval[epoch] = [acc, auc]
            print ("Epoch:", '%d' % (epoch+1), "cost =", "{:.9f}".format(avg_cost),
                    "Training accuracy:", round(acc,3), " Training auc:", round(auc,3))

        if avg_cost <= 0.1:
            print("Early stopping.")
            break

    end =  time.process_time()
    print("运行耗时(s)：",end-start)

    ## Testing cycle
    acc, y_s = sess.run([accuracy, y_score], feed_dict={x: x_test, y: y_test, keep_prob: 1})
    auc = metrics.roc_auc_score(y_test, y_s)
    var_imp = sess.run([var_importance])
    var_imp = np.reshape(var_imp, [n_features1])
    print("*****=====", "Testing accuracy: ", acc, " Testing auc: ", auc, "=====*****")
    #print("*****=====", "y_test: ", y_test, " y_s: ", y_s, "=====*****")

np.savetxt(out_file, var_imp, delimiter=",")
np.savetxt(out_file1, y_test, delimiter=",")
np.savetxt(out_file2, y_s, delimiter=",")