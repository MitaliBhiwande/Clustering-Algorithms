
# coding: utf-8

# In[40]:

import numpy as np
import sys
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import math

def hac(file,number_of_clusters):

    global geneCount
    global expCount
    global geneexp
    global numberofclusters
    
    path=file
    numberofclusters=number_of_clusters
    # read file as lines into genes array
    file = open(path , 'r')
    genes=file.readlines();
    geneCount=len(genes)
    expCount=len(genes[0].split('\t'))
    
    #print ('parameters', geneCount, '|', expCount)
    
    #matrix for storing gene exp
    geneexp=np.zeros((geneCount,expCount),dtype=np.float64 )
    ground_t=[]
    index_dict={}
    for i in range (geneCount):
        geneexp[i]=np.array(genes[i].split('\t'))
        #index_dict[i]=i+1
        index_dict[i]=str(i+1)
        ground_t.append(int(genes[i].split('\t')[1]))
    
    if(numberofclusters>geneCount):
        sys.exit("Not enough values")
        

    # euclidean distance
    distmatrix = squareform(pdist(geneexp[0:geneCount,2:expCount], 'euclidean'))
    #print(distmatrix)
    
    counter=-1

    while counter!=(geneCount-int(numberofclusters)-1):
        
        counter=0
        for k in index_dict:
            if index_dict[k]=='empty':
                counter+=1
        #print('counter', geneCount-int(numberofclusters)-1)
        
        minimum=(distmatrix[distmatrix>0]).min()
        min_index=np.where(distmatrix==minimum)[-1]
        #print(min_index)
        x=min_index[-2]
        y=min_index[-1]
        for i in range(0,geneCount):
        # recompute dist matrix after combining clusters
            distmatrix[x][i]=min(distmatrix[x][i],distmatrix[y][i])
            distmatrix[y][i]= 0     # diagonally mirrored hence redundant value
            distmatrix[i][x]=min(distmatrix[i][x],distmatrix[i][y])
            distmatrix[i][y]= 0
    
        #combine clusters in map and delete old entry
        index_dict[x]=index_dict[x]+','+(index_dict[y])
        index_dict[y]="empty"
        
       
    marker=0
    clusters=[]
    #print(index_dict)
    
    for i in range(numberofclusters):
        clusters.append([])
    for i in range(geneCount):
        if index_dict[i]!='empty':
            clusters[marker].append(index_dict[i])
            marker+=1
    #print (len(clusters))

    finalclusters=[]
    for i in range(len(clusters)):
        finalclusters.append([])
    
    for i in range(len(clusters)):    
        cluster_set=clusters[i]
        temp=[]
        for clust in cluster_set:
            for pt in clust.split(','):
                temp.append(int(pt))
        finalclusters[i]=temp  
    #print(finalclusters)
    return finalclusters,ground_t,


# In[41]:

def hac_labels(clusterlist):
    
    #print (len(clusterlist))
    result=[0]*geneCount
    clustername=1
    
    for labels in clusterlist:
        for l in labels:
            result[l-1]=clustername
        clustername+=1
    
    #print(result)
    return result


# In[42]:

def plot_PCA(geneexp,labels,plot_title):
    data=geneexp[:,2:expCount]
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
     


# In[43]:

def compute_coeffs(original_data,cluster_indices):
    m11=0
    m10=0
    m01=0
    m00=0
    
    for i in range(0,len(original_data)):
        gene1_ground_truth=int(original_data[i:i+1,1:2][0])
        gene1_data=original_data[i:i+1,2:original_data.shape[1]][0]
        for cluster in range(len(cluster_indices)):
            if((i+1) in cluster_indices[cluster]):
                gene1_cluster_number=cluster
                break
        for j in range(0,original_data.shape[0]):
            gene2_ground_truth=int(original_data[j:j+1,1:2][0])
            gene2_data=original_data[j:j+1,2:original_data.shape[1]][0]
            for cluster in range(len(cluster_indices)):
                if(j+1 in cluster_indices[cluster]):
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
                
    jacard_coeff=float(m11)/float((m11+m01+m10))
    print("jaccard_coefficient: ", jacard_coeff)
    rand_index=float((m11+m00))/float((m11+m00+m10+m01))
    print("rand_index: ", rand_index)

    


# In[47]:


h_clusters,ground_t=hac(sys.argv[1],int(sys.argv[2]))
i = 1
for c in h_clusters:
    print('cluster ',i)
    print(c)
    i = i + 1
print('')
compute_coeffs(geneexp,h_clusters)
h_labels=hac_labels(h_clusters)
#print(h_labels)
plot_PCA(geneexp,h_labels, 'Hierarchical clustering PCA plot')


# In[48]:

plot_PCA(geneexp,ground_t,'Ground Truth PCA plot')


# In[ ]:



