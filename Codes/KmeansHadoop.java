import java.util.*;
import java.io.*;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.FileSystem;

@SuppressWarnings("deprecation")
public class KmeansHadoop {
	
	private static final String SPLITTER = "\t| ";
	private static final String OUTPUT_FILE_NAME = "/part-r-00000";
	
	public static Boolean compareCentroids(HashMap<Integer, ArrayList<Double>> oldCentroids, HashMap<Integer, ArrayList<Double>> newCentroids) {
		
		for(Integer centroidIndex : oldCentroids.keySet()) {
			if(!oldCentroids.get(centroidIndex).equals(newCentroids.get(centroidIndex))) {
				return false;
			}
		}
		return true;
	}
	

	public static class Map extends Mapper<LongWritable,Text,IntWritable,Text> {
		
		private static HashMap<Integer, ArrayList<Double>> centroids = new HashMap<Integer, ArrayList<Double>>();	

		protected void setup(Context context) throws IOException, InterruptedException {
			
	       Configuration conf = context.getConfiguration();
	       String centroidFilePath = conf.get("centroidFilePath");
	       FileSystem fileSys = FileSystem.get(conf);
	       BufferedReader buffRdr = new BufferedReader(new InputStreamReader(fileSys.open(new Path(centroidFilePath))));
			String line = buffRdr.readLine();
			ArrayList<Double> centroid = new ArrayList<Double>();
			while (line != null) {
				centroid = new ArrayList<Double>();
				String[] strArr = line.split("\t| ");
				for(int i = 1; i < strArr.length; i++) {
					centroid.add(Double.parseDouble(strArr[i]));
				}
				centroids.put(Integer.parseInt(strArr[0]), centroid);
				line = buffRdr.readLine();
			}
			buffRdr.close();
	   }
		
		public void map(LongWritable key, Text value,Context context) throws IOException,InterruptedException{
			
			String[] geneAttr = value.toString().split(SPLITTER);
			ArrayList<Double> gene = new ArrayList<Double>();
			for(int i = 2; i < geneAttr.length; i++) {
				gene.add(Double.parseDouble(geneAttr[i]));
			}
			Integer centroidId = findNearestCentroid(gene);
			context.write(new IntWritable(centroidId), value);
		}
		
		public static Integer findNearestCentroid(ArrayList<Double> gene) {
			
			Integer minCentroidId = null;
			Double minDistance = null, euclideanDist = null;
			for(Integer centroidIndex : centroids.keySet()) {
				euclideanDist = 0.0;
				for(int i = 0; i < gene.size(); i++) {
					euclideanDist += Math.pow(gene.get(i)-centroids.get(centroidIndex).get(i), 2);
				}
				euclideanDist = Math.sqrt(euclideanDist);
				if(minDistance == null || euclideanDist < minDistance) {
					minDistance = euclideanDist;
					minCentroidId = centroidIndex;
				}
			}	
			return minCentroidId;
		}
		
	}
	
	public static class Reduce extends Reducer<IntWritable,Text,IntWritable,Text> {
		
		public void reduce(IntWritable key, Iterable<Text> values, Context context) throws IOException,InterruptedException {
			
			ArrayList<ArrayList<Double>> cluster = new ArrayList<ArrayList<Double>>();
			ArrayList<Double> gene = new ArrayList<Double>();
			String[] geneAttr = null;
			String out = "";
			Configuration conf = context.getConfiguration();
			String isFinal = conf.get("isFinal");
	        System.out.println(isFinal);
	        if(isFinal.equals("false")) {
				for(Text strGene: values){
					gene = new ArrayList<Double>();
					geneAttr = strGene.toString().split("\t| ");
					for(int i = 2; i < geneAttr.length; i++) {
						gene.add(Double.parseDouble(geneAttr[i]));
					}
					cluster.add(gene);
				}
				ArrayList<Double> centroids = calculateCentroids(cluster);
				if(centroids.size() > 0) {
					for(Double attr : centroids) {
						out += String.valueOf(attr) + " ";
					}
					context.write(key, new Text(out.substring(0, out.length() - 1)));
				}
	        }else if(isFinal.equals("true")) {
		       	for(Text strGene: values){
		       		geneAttr = strGene.toString().split("\t| ");
		       		out = geneAttr[0] + " ";
		       		for(int i = 2; i < geneAttr.length; i++) {
		       			out += geneAttr[i] + " ";
						}
		       		context.write(key, new Text(out.substring(0, out.length() - 1)));
				}
	        }
		}
		
		
		public static ArrayList<Double> calculateCentroids(ArrayList<ArrayList<Double>> cluster){
		

			Integer attrLen = cluster.get(0).size();
			ArrayList<Double> centroids = new ArrayList<Double>(Arrays.asList(new Double[attrLen]));
			Double attrSum = null;
			for(int i = 0; i < attrLen; i++) {
				attrSum = 0.0;
				for(ArrayList<Double> gene : cluster) {
					attrSum += gene.get(i);
				}
				centroids.set(i, attrSum/cluster.size());
			}
			return centroids;
		}
	}

	public static void main(String[] args) throws Exception {

		String geneDataPath = args[0];
		String outputPath = args[1];
		String initialCentroidFilePath = args[2];
		Integer iterator = 0;
		Boolean completed = false, kmeansCompleted = false;
		HashMap<Integer, ArrayList<Double>> prevCentroids = null;
		HashMap<Integer, ArrayList<Double>> newCentroids = null;
		while(true) {
			System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Iteration Number " + iterator);
			Configuration conf= new Configuration();
			if(iterator == 0) {
				conf.set("centroidFilePath", initialCentroidFilePath);
			}else {
				conf.set("centroidFilePath", outputPath + (iterator - 1) + OUTPUT_FILE_NAME);
			}
			if(kmeansCompleted) {
				conf.set("isFinal", "true");
			}else {
				conf.set("isFinal", "false");
			}
			Job job = new Job(conf,"KMeans");
			job.setJarByClass(KmeansHadoop.class);
			job.setMapperClass(Map.class);
			job.setReducerClass(Reduce.class);
			job.setOutputKeyClass(IntWritable.class);
			job.setOutputValueClass(Text.class);
			job.setInputFormatClass(TextInputFormat.class);
			job.setOutputFormatClass(TextOutputFormat.class);
			FileInputFormat.addInputPath(job, new Path(geneDataPath));
			if(kmeansCompleted) {
				FileOutputFormat.setOutputPath(job, new Path(outputPath + "_final"));
			}else {
				FileOutputFormat.setOutputPath(job, new Path(outputPath + iterator));
			}
			completed = job.waitForCompletion(true);
			if (completed) {
				if(kmeansCompleted) {
					break;
				}
				prevCentroids = new HashMap<Integer, ArrayList<Double>>();
				newCentroids = new HashMap<Integer, ArrayList<Double>>();
				Path filePath = null;
				ArrayList<Double> centroid = new ArrayList<Double>();
				BufferedReader buffRdr = null;
				String line = "";
				FileSystem fileSys = FileSystem.get(new Configuration());
				
				// Read the previous centroids from reducer output.
				if(iterator == 0) {
					filePath = new Path(initialCentroidFilePath);
				}else {
					filePath = new Path(outputPath + (iterator - 1) + OUTPUT_FILE_NAME);
				}
				buffRdr = new BufferedReader(new InputStreamReader(fileSys.open(filePath)));
				line = buffRdr.readLine();
				while (line != null) {
					centroid = new ArrayList<Double>();
					String[] strArr = line.split("\t| ");
					for(int i = 1; i < strArr.length; i++) {
						centroid.add(Double.parseDouble(strArr[i]));
					}
					prevCentroids.put(Integer.parseInt(strArr[0]), centroid);
					line = buffRdr.readLine();
				}
				buffRdr.close();
				// Read the new computed centroids from reducer output.
				filePath = new Path(outputPath + iterator + OUTPUT_FILE_NAME);
				buffRdr = new BufferedReader(new InputStreamReader(fileSys.open(filePath)));
				line = buffRdr.readLine();
				while (line != null) {
					centroid = new ArrayList<Double>();
					String[] strArr = line.split("\t| ");
					for(int i = 1; i < strArr.length; i++) {
						centroid.add(Double.parseDouble(strArr[i]));
					}
					newCentroids.put(Integer.parseInt(strArr[0]), centroid);
					line = buffRdr.readLine();
				}
				buffRdr.close();
				if(compareCentroids(prevCentroids, newCentroids)) {
					kmeansCompleted = true;
				}
				iterator++;
	       }
		}
		System.exit(kmeansCompleted ? 0 : 1);
	}
}