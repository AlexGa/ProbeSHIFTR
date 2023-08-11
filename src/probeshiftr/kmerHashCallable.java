package probeshiftr;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.HashMap;
import java.util.concurrent.Callable;

class kmerHashCallable implements Callable<HashMap<String,Integer>> {

	private Sequence seq;
	private int kmer;
	private HashMap<String,Integer> oligoHashMap;
	
	
	kmerHashCallable(Sequence seq, int kmer){
		
		super();
		this.seq = seq;
		this.kmer = kmer;
	}
	
	@Override
	public HashMap<String,Integer> call() throws Exception{
		
		try{
			
			oligoHashMap = Oligo.kmerHash(seq, kmer);	
			System.out.println("Done");
			
		}catch(Exception e){	

			System.out.println(e.getStackTrace());
			System.out.println(e.getMessage());
			
			StringWriter sw = new StringWriter();
			PrintWriter pw = new PrintWriter(sw);
			e.printStackTrace(pw);
			String sStackTrace = sw.toString(); // stack trace as a string
			System.out.println(sStackTrace);
		}
		return oligoHashMap;
	}
}