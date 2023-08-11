package probeshiftr;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;

public class LowComplexityFilter implements  OligoFilter {

	String faFile;
	HashMap<String, Integer> patternHash;
	int threads;
	int kmerLength;
	double mean_entropy;
	double var_entropy;
	int Npatterns; // set in NormFactors
	HashMap<String, Double> kmerProb;
	double cutOff;
	
	LowComplexityFilter(String faFile, int threads, int kmerLength, double complexityCutoff){
		
		this.faFile = faFile;
		this.threads = threads;
		this.kmerLength = kmerLength;
		this.cutOff = complexityCutoff;
		
		patternHash = getOligoHash();
		patternSize();
		setNormFactors();
		
	}
	
	private void patternSize() {
		
		this.Npatterns = 0;
		
		for(String kmer: this.patternHash.keySet()) {
			
			this.Npatterns += this.patternHash.get(kmer);
		}
	}
	
	private void setNormFactors() {
	
		
		this.kmerProb = new HashMap<String, Double>(); 
		this.mean_entropy = 0.0d;
		this.var_entropy = 0.0d;
		
		int Nkmers = patternHash.keySet().size();
		
		for(String kmer: patternHash.keySet()) {
			
			kmerProb.put(kmer, (double) patternHash.get(kmer) / (double) Npatterns);
			
			mean_entropy += kmerProb.get(kmer);
		}
		
		mean_entropy /= (double) Nkmers;
		
		
		for(String kmer: patternHash.keySet()) {
			
			var_entropy += Math.pow(kmerProb.get(kmer) - mean_entropy, 2);
		}
		
		var_entropy /= (double) Nkmers;
		
	}
	
	private HashMap<String, Integer> getOligoHash(){
		
		HashMap<String, Sequence> seqs;
		HashMap<String, Integer> oligoHash = null;
		
		try{

			System.out.println("Read genome file");
			
			seqs = SeqIO.readFastaGenom(faFile);
			
			System.out.println("Read " + seqs.size() + " sequences");
			
			String[] idArray = new String[seqs.keySet().size()];
			seqs.keySet().toArray(idArray);
			
			ExecutorService threadPool = Executors.newFixedThreadPool(threads);
			 
			// submit jobs to be executing by the pool
			ArrayList<Future<HashMap<String, Integer>>> submitList = new ArrayList<Future<HashMap<String, Integer>>>();
			 
			// store results from each thread in the list
			ArrayList<HashMap<String, Integer>> kmerCountList = new ArrayList<HashMap<String, Integer>>();
			
			try {
				
				for (int i = 0; i < idArray.length; i++) {
					
					System.out.println("Create oligoHash of " + i + " chromosome");
					
					Sequence backgroundSequence = seqs.get(idArray[i]);
					
					Callable<HashMap<String, Integer>> worker = new kmerHashCallable(backgroundSequence, kmerLength);
					Future<HashMap<String, Integer>> submit = threadPool.submit(worker);
					submitList.add(submit);
				}
				
				int nbJobs = submitList.size();
				int queueSize = nbJobs;
				 
				do{
					queueSize = ((ThreadPoolExecutor) threadPool).getQueue().size();	
					System.out.println(queueSize);
					
				}while(queueSize > 0);
				 
				 for(int i=0; i<submitList.size(); i++){
				 
					kmerCountList.add(submitList.get(i).get());
				 }
				  
				threadPool.shutdown();
				
				oligoHash = combineCounters(kmerCountList);

			} catch (Exception e) {
				  
				System.out.println(e.getMessage());
			}
			
		}catch(IOException e){
			
			System.out.println(e.getMessage());
		}
		
		
		
		// Combine kmerCount ArrayList into one HashMap
		
		return oligoHash;
		
	}
	
	private HashMap<String, Integer> combineCounters(ArrayList<HashMap<String, Integer>> counterList){
	
		HashMap<String, Integer> combinedHash = new HashMap<String, Integer>();
		
		if(counterList.size() >= 1) {
			
			combinedHash = counterList.remove(0);
		
		}
		
		if(counterList.size() >= 1) {
			
			for(HashMap<String, Integer> singleHash: counterList) {
				
				for(String key: singleHash.keySet()) {

					if(combinedHash.containsKey(key)) {
						
						combinedHash.put(key, combinedHash.get(key) + singleHash.get(key));
					}
				}
			}
		
		}
		
		return combinedHash;
	}
	
	private double informationContent(String oligo) {
		
		double infContent = 0.0d;
		
		HashMap<String, Integer> oligoHash = Oligo.kmerHash(new Sequence(null, oligo), kmerLength);
		
		for(int l = 0; l < oligo.length() - kmerLength + 1; l++) {
		
			System.out.println(oligo.substring(l, l + kmerLength));
			
			infContent += kmerEntropy(oligoHash, oligo.substring(l, l + kmerLength));

		}

		infContent = normalizeInfContent(oligo.length(), infContent);
		
		return infContent;
		
	}
	
	private double normalizeInfContent(int oligoLength, double infContent) {
	
		double scaling_factor = (double) (oligoLength - kmerLength + 1);
		
		return (infContent - scaling_factor * mean_entropy)/(scaling_factor * Math.sqrt(var_entropy));

	}
	
	private double kmerEntropy(HashMap<String, Integer> oligoHash, String kmer) {
		
		double p = kmerProb.get(kmer);
		
		double kmerEnt = p * Math.log(p)/Math.log(2);
		
		return kmerEnt;
	}
	
	
	@Override
	public boolean filterOligo(String oligo) {
		
		boolean aboveCutoff = false;
		
		double infContent = informationContent(oligo);
		
		aboveCutoff = infContent > this.cutOff;
		
		// TODO Auto-generated method stub
		return aboveCutoff;
	}
	
	public HashMap<String, Double> getkmerProb(){
		
		return kmerProb;
	}
}
