################################  K-MEANS CLUSTERING ############################################

To run the code just pass the data file path as an input command line parameter
- Below is a sample execution command:

if specifying the cluster centers are specified, the command is:
>> python KmeanClustering.py cho.txt 5 5 10 15 20 25 10

Here number of clusters = 5, gene IDs of centroids: 5 10 15 20 25 , number of iterations= 10

if the initial cluster centroids are not entered:
>> python KmeanClustering.py cho.txt 5 random 10
Here number of clusters = 5, gene IDs of centroids:random , number of iterations= 10

Input:

- Number of clusters, gene IDs of centroids or random and number of iterations.

Output:

- The code will output the clusters and the gene id's within that cluster.
- The external index i.e Jaccard coefficient and Rand Index
- Plots using PCA for the clustering algorithm results and the ground truth clusters.


################################  HIERARCHICAL AGGLOMERATIVE CLUSTERING ############################

- To run the code pass the data file path and the number of clusters as an input using the command line.
- Below is a sample execution command:

>> python HAC.py cho.txt 5
>> python HAC.py iyer.txt 10

Input:

- The txt file that contains the gene data with all the attributes.
- Number of clusters to be formed.

Output:

- The code will output the clusters and the gene id's within that cluster.
- The external index i.e Jaccard coefficient and Rand Index
- Plots using PCA for the clustering algorithm results and the ground truth clusters.


############################################## DBSCAN ################################################

- To run the code pass the data file path, epsilon and the minimum points criteria as an input using command line.
- Below is a sample execution command:

>> python DBSCAN.py cho.txt 1.04 4
>> python DBSCAN.py iyer.txt 1.5 2

Input:

- The txt file that contains the gene data with all the attributes.
- Epsilon value for identifying the neighbors for any point.
- Minimum number of neighbors the point should have to be present in a cluster. Else it will considered as noise.

Output:

- The code will output the clusters and the gene id's within that cluster.
- The external index i.e Jaccard coefficient and Rand Index
- Plots using PCA for the clustering algorithm results.


################################# KMEANS USING HADOOP MAPREDUCE ########################################

Pre-requisite:
* A single folder with following files in it:
1. Kmeans.jar [Hadoop Mapreduce code to perform Kmeans on given dataset] 
2. KmeansHadoop.sh [A script file to execute all the commands in a sequence]
3. Data text file [A text file with complete gene data]
4. Centroid text file [A text file with initial centroids specified manually. First column represents a unique centroid number for each record and all others are the attributes.]
5. KmeansHadoop.py [A python file to calculate jacard coefficient and plot graphs]

Steps:
1. Navigate to the directory of the folder with all specific files present.
2. Start hadoop dfs and yarn using standard commands to activate Datanode, Namenode, ResourceManager and Nodemanager.
3. Execute shell script file with appropriate parameters. Below is the list of parameters need to be passed. All the files specified should be present in the same directory
	i. Data file name.
	ii. Jar file name for hadoop mapreduced code.
	iii. Centroid file name. Centroid file needs to be populate manually. Initial centroids needs to be added to the file.

Example of execution:
./KmeansHadoop.sh new_dataset_1.txt centroids_new.txt Kmeans.jar