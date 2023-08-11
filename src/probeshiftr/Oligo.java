package probeshiftr;

import java.util.HashMap;


public class Oligo {

	Oligo(){
		
	}
	
	public static HashMap<String, Integer> kmerHash(Sequence seq, int k){
		
		
		HashMap<String, Integer> kmerTable = new HashMap<String, Integer>();
		
		for(int i = 0; i < seq.getSeq().length() - k; i++) {
			
			String kmer = seq.getSeq().substring(i, i + k);
			Integer count = 1;
			
			if(kmerTable.containsKey(kmer)) {
				
				count = kmerTable.get(kmer) + 1;
				
			}
			kmerTable.put(kmer, count);
		}
		
		return(kmerTable);
	}	
}
