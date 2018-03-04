
# coding: utf-8
"""
Description: Compute K mean clusters and generate scatter plots using PCA for a given dataset 
Group Members: Mitali Bhiwande | Sumedh Ambokar | Tejasvi Sankhe
"""
# In[207]:
import sys
import numpy as np
import random
from scipy.spatial import distance
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


# In[208]:
# Function to fetch data from input file and compute initial centroids if needed.
def fetch_data(filename, num_of_clusters, centroidsArr, num_of_iterations):
    centroid_list=[]
    global iteration_num
    iteration_num=0
    with open(filename) as textFile:
        lines = [line.replace("\n","").split("\t") for line in textFile]
        original_data = np.array(lines, dtype='float')
        if(centroidsArr == []):
            centroid_array=np.array(get_random_IDs(original_data,num_of_clusters))
        else:
            centroid_array=np.array(centroidsArr)
        for i in range(len(centroid_array)):
            centroid_list.append(original_data[centroid_array[i]-1])
        centroids=np.asarray(centroid_list)
        centroids=centroids[:,2:]
        compute_kmean(original_data,centroids,num_of_iterations,num_of_clusters,iteration_num)


# In[209]:
# This function provides random locations of gene data as initial clusters
def get_random_IDs(original_data,number_of_clusters):
    return (np.random.choice(original_data.shape[0], number_of_clusters, replace=False))


# In[210]:
#Function computes the distances, comapres it and assigns points to clusters as per minimum distance.
def compute_kmean(original_data, centroids,num_of_iterations,num_of_clusters,iteration_num):
    clusters=[[] for _ in range(len(centroids))]
    cluster_indices=[[] for _ in range(len(centroids))]
    for i in range(0,original_data.shape[0]):
        gene_data=original_data[i:i+1,2:original_data.shape[1]]
        dist=[]
        for j in range(0,centroids.shape[0]):
            dist.append(distance.euclidean(gene_data[0],centroids[j]))
        clusters[dist.index(min(dist))].append(gene_data[0])
        cluster_indices[dist.index(min(dist))].append(i+1)
    iteration_num=iteration_num + 1
    if(iteration_num <= num_of_iterations):
        print("Iteration Number : ", iteration_num)
    compute_centroids(centroids,clusters,original_data,cluster_indices,num_of_iterations,num_of_clusters,iteration_num)           


# In[211]:
#Computes the centroids of new clusters formed and checks for convergence of kmeans
def compute_centroids(centroids,clusters,original_data,cluster_indices,num_of_iterations,num_of_clusters,iteration_num):
    average=[]
    for i in range(len(clusters)):
        average.append(list(np.average(clusters[i], axis=0)))
    #average = [[float("%.4f"%item) for item in cluster_ele] for cluster_ele in average]
    z = [tuple(y) for y in centroids]
    x = [tuple(y) for y in average]
    a=set(x) - set(z)
    if(len(a)==0 or iteration_num==num_of_iterations):
        compute_coeffs(original_data,cluster_indices)
        draw_plots(original_data,clusters,num_of_clusters)
    else:
        compute_kmean(original_data,np.array(average),num_of_iterations,num_of_clusters,iteration_num)


# In[212]:
# computes the external index for similarity measure
def compute_coeffs(original_data,cluster_indices):
    for i in range(len(cluster_indices)):
        print("Cluster : ",i+1)
        print(cluster_indices[i])
    m11=0
    m10=0
    m01=0
    m00=0
    for i in range(0,len(original_data)):
        gene1_ground_truth=int(original_data[i:i+1,1:2][0])
        gene1_data=original_data[i:i+1,2:original_data.shape[1]][0]
        for cluster in range(len(cluster_indices)):
            if(i+1 in cluster_indices[cluster]):
                gene1_cluster_number=cluster
                break
        for j in range(1,original_data.shape[0]):
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


# In[213]:
#plots clusters 
def draw_scatter_plot(x,y,clusters,num_of_clusters, plot_name):
    clusters = [int(ele[0]) for ele in clusters]
    categories = np.unique(clusters)
    colors = np.linspace(0, 1, len(categories))
    use_colours = dict(zip(categories, colors))
    unique = list(use_colours.values())
    plot_size = plt.rcParams["figure.figsize"]
    plot_size[0] = 12
    plot_size[1] = 9
    plt.rcParams["figure.figsize"] = plot_size
    colors = [plt.cm.jet(float(i)/max(unique)) for i in unique]
    for u in use_colours.keys():
        xi = [x[j] for j  in range(len(x)) if clusters[j] == u]
        yi = [y[j] for j  in range(len(x)) if clusters[j] == u]
        plt.scatter(xi, yi, c=colors[unique.index(use_colours.get(u))], label=str(u))
    #plt.scatter(centroids[:, 0][i], centroids[:, 1][i],marker='x', s=500, linewidths=5,color='black', zorder=10)
    plt.xlabel("PC 1")
    plt.ylabel("PC 2")
    title = "K-Means Clustering for "+str(num_of_clusters)+" centroids" + " on " + plot_name
    plt.title(title)
    plt.legend()
    plt.show()


# In[214]:
# makes the data ready for plot by applying PCA
def draw_plots(original_data,clusters,num_of_clusters):
    data=[]
    ground_truth=[]
    pca = PCA(n_components=2)
    for i in range(0,original_data.shape[0]):
        ground_truth.append((list(original_data[i:i+1,1:2][0])))
        data.append((original_data[i:i+1,2:original_data.shape[1]][0]))
    pca.fit(data)
    new_data=np.array(pca.transform(data))
    x=[]
    y=[]
    for i in range(0,new_data.shape[0]):
        x.append((new_data[i:i+1,0:1][0]))
        y.append(new_data[i:i+1,1:2][0])
    
    draw_scatter_plot(np.array(x),np.array(y),ground_truth,num_of_clusters, "Original Data")
    
    cluster_data=[]
    cluster_number=[]
    for i in range(len(clusters)):
        for j in range(len(clusters[i])):
            cluster_data.append(clusters[i][j])
            cluster_number.append([i+1])
    pca.fit(cluster_data)
    new_cluster_data=np.array(pca.transform(cluster_data))
    clusterx=[]
    clustery=[]
    for i in range(0,new_cluster_data.shape[0]):
        clusterx.append((new_cluster_data[i:i+1,0:1][0]))
        clustery.append(new_cluster_data[i:i+1,1:2][0])
    draw_scatter_plot(np.array(clusterx),np.array(clustery),cluster_number,num_of_clusters, "Kmeans Output")
    
"""

# In[215]:

fetch_data("cho.txt",5,[5,25,32,100,132],10)
fetch_data("cho.txt",5,[],100)


# In[216]:

fetch_data("iyer.txt",11,[45,2,56,23,1,67,3,11,5,89,180],10)
fetch_data("iyer.txt",11,[],100)


# In[217]:

fetch_data("new_dataset_1.txt",3,[3,5,9],10)
fetch_data("new_dataset_1.txt",3,[],10)


# In[ ]:"""

if len(sys.argv) > 1:

    filePath = sys.argv[1]
    clusters_no = int(sys.argv[2])
    centroid_id = []
    iterations = -1
    index = 0
    if sys.argv[3] == "random":
        iterations = int(sys.argv[4])
    else:    
        for x in range(clusters_no):
            centroid_id.append(int(sys.argv[x + 3]))
            index = x + 3
        iterations = int(sys.argv[index + 1])
    fetch_data(filePath, clusters_no, centroid_id, iterations)
else:
    print("Not sufficient Inputs")



