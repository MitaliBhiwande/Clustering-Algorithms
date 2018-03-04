
# coding: utf-8

# In[161]:

import sys
import numpy as np
import random
from scipy.spatial import distance
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


# In[162]:

def fetch_data(filename,eps,minPts):
    with open(filename) as textFile:
        lines = [line.replace("\n","").split("\t") for line in textFile]
        original_data = np.array(lines, dtype='float')
        data=[]
        for i in range(0,original_data.shape[0]):
            data.append((original_data[i:i+1,2:original_data.shape[1]][0]))
        compute_dbscan((data),original_data,eps,minPts)


# In[163]:

def regionQuery(data,gene_data, eps):
    neighbors=[]
    for i in range(len(data)):
        if(distance.euclidean(data[gene_data],data[i]) < eps):
            neighbors.append(i)
    #neighbors.append(gene_data)
    return neighbors


# In[164]:

def expandCluster(data, gene_data, neighborPts, eps,minPts,C, label):
    label[gene_data]=(C)
    i=0
    while i<len(neighborPts):
        Pt=neighborPts[i]
                
        if(label[Pt]==-1):
            label[Pt]=C
        if(label[Pt]=="unvisited"):
            label[Pt]=C
            new_neighborPts=list(regionQuery(data,Pt,eps))
            if(len(new_neighborPts) >= minPts):
                neighborPts = list((neighborPts + new_neighborPts))
        i=i+1


# In[165]:

def plot_PCA(geneexp,labels,plot_title):
    data=geneexp[:,2:]
    pca = PCA(n_components=2)
    data = np.matrix(data).T
    pca.fit(data)
    data_pca = pca.components_
    title =  plot_title
    plt.figure(figsize=(8,6))
    px=data_pca[0,]
    py=data_pca[1,]
    unique = list(set(labels))
    colors = [plt.cm.jet(float(i)/max(unique)) for i in unique]
    for i, u in enumerate(unique):
        xi = [px[j] for j  in range(len(px)) if labels[j] == u]
        yi = [py[j] for j  in range(len(px)) if labels[j] == u]
        plt.scatter(xi, yi, c=colors[i], label=str(u))
    
    plt.legend()
    plt.xlabel("PC 1")
    plt.ylabel("PC 2")
    plt.title(title)
    plt.show()


# In[166]:

def compute_dbscan(data,original_data,eps,minPts):
    C=0
    label=list(data)
    for i in range(len(label)):
        label[i]="unvisited"
    for i in range(len(data)):
        gene_data=i
        if(label[i]=="unvisited" or label[i]==-1): 
            neighborPts = list(regionQuery(data,gene_data, eps))
            if(len(neighborPts) < minPts):
                label[i]=-1
            else:
                C=C+1
                expandCluster(data,gene_data,neighborPts, eps, minPts,C, label)
    print((label))
               
   
    compute_coeffs(original_data,label)
    plot_PCA(original_data,label,"Clustering results for DBSCAN")


# In[167]:

def compute_coeffs(original_data,label):
    m11=0
    m10=0
    m01=0
    m00=0
    for i in range(0,len(original_data)):
        gene1_ground_truth=int(original_data[i:i+1,1:2][0])
        gene1_data=original_data[i:i+1,2:original_data.shape[1]][0]
        gene1_cluster_number=label[i]
        for j in range(1,original_data.shape[0]):
            gene2_ground_truth=int(original_data[j:j+1,1:2][0])
            gene2_data=original_data[j:j+1,2:original_data.shape[1]][0]
            gene2_cluster_number=label[j]
            if(gene1_ground_truth == gene2_ground_truth and gene1_cluster_number == gene2_cluster_number):
                m11=m11+1
            elif(gene1_ground_truth == gene2_ground_truth and not(gene1_cluster_number == gene2_cluster_number)):
                m10=m10+1
            elif(not(gene1_ground_truth == gene2_ground_truth) and gene1_cluster_number == gene2_cluster_number):
                m01=m01+1
            elif(not(gene1_ground_truth == gene2_ground_truth) and not(gene1_cluster_number == gene2_cluster_number)):
                m00=m00+1
                
    jacard_coeff=float(m11)/float((m11+m01+m10))
    print("jaccard_coefficient: ", jacard_coeff)
    rand_index=float((m11+m00))/float((m11+m00+m10+m01))
    print("rand_index: ", rand_index)

    


# In[168]:

#fetch_data("iyer.txt",1.5,3)
fetch_data(sys.argv[1],float(sys.argv[2]),int(sys.argv[3]))


# In[ ]:




# In[ ]:



