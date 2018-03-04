#!/bin/bash

hadoop fs -test -e /Kmeans
if [ $? == 0 ]
    then
        hdfs dfs -rm -r /Kmeans
fi
hdfs dfs -mkdir /Kmeans
hdfs dfs -mkdir /Kmeans/Output
hdfs dfs -put $1 /Kmeans
hdfs dfs -put $2 /Kmeans
hadoop jar $3 /Kmeans/$1 /Kmeans/Output/output /Kmeans/$2
hdfs dfs -get /Kmeans/Output/output_final/part-r-00000 ./
python KmeansHadoop.py $1 part-r-00000
