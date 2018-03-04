import numpy as np
import sys
from scipy.spatial import distance
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

data_file_path = sys.argv[1]
output_file_path = sys.argv[2]

def fetch_data(data_file_path, output_file_path):
    original_data = []
    cluster_indices = []
    clusterAttr = []
    current_cluster = "";
    cluster = []
    flag = True
    with open(data_file_path) as textFile:
        lines = [line.replace("\n","").split() for line in textFile]
        original_data = np.array(lines, dtype='float')
    clusters = [None]*len(original_data)
    with open(output_file_path) as textFile:
        lines = [line.replace("\n","").split() for line in textFile]
    for line in lines :
        if(line[0] == current_cluster or flag):
            cluster.append(int(line[1]) - 1)
            if(flag):
                current_cluster = line[0]
            flag = False
        else:
            current_cluster = line[0]
            cluster_indices.append(cluster)
            cluster = []
            cluster.append(int(line[1]) - 1)
        clusters[int(line[1]) - 1] = int(line[0])
    cluster_indices.append(cluster)
    return original_data, cluster_indices, clusters

def compute_coeffs(original_data,cluster_indices):
    m11=0
    m10=0
    m01=0
    m00=0
    for i in range(0,len(original_data)):
        gene1_ground_truth=int(original_data[i:i+1,1:2][0])
        gene1_data=original_data[i:i+1,2:original_data.shape[1]][0]
        for cluster in range(len(cluster_indices)):
            if(i in cluster_indices[cluster]):
                gene1_cluster_number=cluster
                break
        for j in range(1,original_data.shape[0]):
            gene2_ground_truth=int(original_data[j:j+1,1:2][0])
            gene2_data=original_data[j:j+1,2:original_data.shape[1]][0]
            for cluster in range(len(cluster_indices)):
                if(j in cluster_indices[cluster]):
                    gene2_cluster_number=cluster
                    break
            if(gene1_ground_truth == gene2_ground_truth and gene1_cluster_number == gene2_cluster_number):
                m11=m11+1
            elif(gene1_ground_truth == gene2_ground_truth and not(gene1_cluster_number == gene2_cluster_number)):
                m10=m10+1
            elif(not(gene1_ground_truth == gene2_ground_truth) and gene1_cluster_number == gene2_cluster_number):
                m01=m01+1
            elif(not(gene1_ground_truth == gene2_ground_truth) and not(gene1_cluster_number == gene2_cluster_number)):
                m00=m00+1
                
    jacard_coeff=float(m11)/float(m11+m01+m10)
    print("jaccard_coefficient: ", jacard_coeff)
    rand_index=float(m11+m00)/float(m11+m00+m10+m01)
    print("rand_index: ", rand_index)  

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

original_data, cluster_indices, clusters = fetch_data(data_file_path, output_file_path)
compute_coeffs(original_data,cluster_indices)
plot_PCA(original_data,clusters,"KmeansHadoop")