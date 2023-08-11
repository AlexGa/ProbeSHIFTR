package probeshiftr;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;


public class Gene extends Feature{

	private ArrayList<MetaFeature> transcripts = null;

	public Gene() {
		super();
		// TODO Auto-generated constructor stub
	}

	public Gene(Feature feature) {
		super(feature);
		// TODO Auto-generated constructor stub
	}

	public Gene(String id, String chromosome, String source, String strand, int start,
			int end, String score, String offset, String description) {
		super(id, chromosome, source, strand, start, end, score, offset, description);
		// TODO Auto-generated constructor stub
	}

	public Gene(String id, String chromosome, String source, String strand, int start, int end) {
		super(id, chromosome, source, strand, start, end);
		// TODO Auto-generated constructor stub
	}

	public Gene(String id, String chromosome, String source,
			String strand, int start, int end, String score,
			String offset, String description, ArrayList<MetaFeature> transcripts) {
		super(id, chromosome, source, strand, start, end, score, offset, description);
		this.transcripts = transcripts;
	}
	
	public ArrayList<MetaFeature> getTranscripts() {
		return transcripts;
	}

	public void setTranscripts(ArrayList<MetaFeature> transcripts) {
		this.transcripts = transcripts;
	}
	
	/* Compare exon structure of each transcript and return index of transcript that matches
	 * 
	 * return -1 --> no similar transcript is found
	 * */
	public int isTranscriptOf(MetaFeature transcript){
		
		int tx_index = -1;
		
		for(int i=0; i<this.transcripts.size(); i++){
			boolean is_tx = this.transcripts.get(i).equals(transcript);
			if(is_tx){
				tx_index = i;
				break;
			}
		}
		/*Iterator<MetaFeature> txIter = this.transcripts.iterator();
		while(txIter.hasNext()){
			MetaFeature known_tx = txIter.next();
			boolean is_tx = known_tx.equals(transcript);
			if(is_tx != -1){
				break;
			}
		}*/
		
		return tx_index;
	}
	
	public GenomicLocation[] overlappingExonIntervals(int position) {
		
		GenomicLocation[] locArr = new GenomicLocation[transcripts.size()];
		
		int i = 0;
		for(MetaFeature transcript: this.getTranscripts()) {
			locArr[i] = transcript.location(position);
			i++;
		}
		
		return locArr;
	}
	
	public void addTranscript(MetaFeature transcript){
		
		if(this.transcripts == null) {
			
			this.transcripts = new ArrayList<MetaFeature>();
		
		}
		
		this.transcripts.add(transcript);
	}
	
	public void addTranscript(MetaFeature transcript, String source, String info){
		
		HashMap<String, String> geneDescHash = this.description2Hash();
		HashMap<String, String> txDescHash = transcript.description2Hash();
		
		for(String gkey: geneDescHash.keySet()){
			for(String tkey: txDescHash.keySet()){
				
				if(gkey.equals(tkey)){
					txDescHash.put(tkey, geneDescHash.get(tkey));
				}
			}
		}
		txDescHash.put("info", info);
		txDescHash.put("transcript_source", source);
		
		transcript.Hash2Description(txDescHash);
		transcript.setSource(source);
		
		if(transcript.getStart() < this.getStart()){
			this.setStart(transcript.getStart());
		}
		
		if(transcript.getEnd() > this.getEnd()){
			this.setEnd(transcript.getEnd());
		}
		
		transcript.updateSubFeatureDescriptions("gene_id", geneDescHash.get("gene_id"));
		this.transcripts.add(transcript);
	}
	
	public String GTFoutput(){
		
		String rows = this.getChromosome()+"\t"+this.getSource()+"\t"+"gene"+"\t"+this.getStart()+"\t"+this.getEnd()+"\t"+this.getScore()+"\t"+this.getStrand()+"\t"+this.getFrame()+"\t"+this.getDescription()+"\n";
		
		if(this.transcripts == null){
			
			System.err.println(this.getId().length() + "has no transcripts!!! " +this.getChromosome()+"\t"+ this.getStart()+"\t"+this.getEnd());
		
		}else{
		
			Iterator<MetaFeature> txIter = this.transcripts.iterator();
			
			while(txIter.hasNext()){
				
				MetaFeature tx = txIter.next();
				tx.setCdsFrames();
				rows += tx.printGTFRows();
				
			} 
		}
		
		
		return rows;
	}
	
	/** Check if gene is just as long as minimum start and maximum end of transcripts*/
	public void CheckInterval() {
	
		Interval2d longestInterval = this.getLongestInterval();
		if(!longestInterval.equals(this.getInterval())) {
			this.setInterval(longestInterval);
		}
		
	}
	
	private Interval2d getLongestInterval() {
		
		int minStart = Integer.MAX_VALUE;
		int maxEnd = Integer.MIN_VALUE;
		
		try {
		for(MetaFeature transcript: this.getTranscripts()) {
			
			transcript.CheckInterval();
			if(minStart > transcript.getStart()) {
				minStart = transcript.getStart();
			}
			
			if(maxEnd < transcript.getEnd()) {
				maxEnd = transcript.getEnd();
			}
		}
		}catch(NullPointerException e) {
			
			System.err.println(this.getId() + " failed interval check!");
			System.exit(1);
		}
		
		return new Interval2d(minStart, maxEnd);
	}
	
 	public void updateTranscript(int index, ArrayList<Feature> exons, String source){
		
		MetaFeature tx = this.transcripts.get(index);
		tx.updateExons(exons);
		
		if(tx.getStart() < this.getStart()){
			this.setStart(tx.getStart());
		}
		
		if(tx.getEnd() > this.getEnd()){
			this.setEnd(tx.getEnd());
		}

		HashMap<String, String> geneDescHash = this.description2Hash();
		HashMap<String, String> txDescHash = tx.description2Hash();
		
		for(String gkey: geneDescHash.keySet()){
			for(String tkey: txDescHash.keySet()){
				
				if(gkey.equals(tkey)){
					txDescHash.put(tkey, geneDescHash.get(tkey));
				}
				
			}
		}
		
		txDescHash.put("transcript_source", source);
		tx.updateSubFeatureDescriptions("gene_id", geneDescHash.get("gene_id"));
		//tx.setSource(source);
		
		this.transcripts.remove(index);
		this.transcripts.add(index, tx);
		
	}
 	
 	public boolean containsTranscript(String id){
 		
 		boolean contains = false;
 		
 		for(MetaFeature tx: this.getTranscripts()){
 			
 			if(tx.getId().equals(id)){
 				contains = true;
 				break;
 			}
 		}
 		
 		return contains;
 	}
 	
 	public MetaFeature getTranscriptByName(String id){
 		
 		
 		for(MetaFeature tx: this.getTranscripts()){
 			
 			if(tx.getId().equals(id)){
 				return tx;
 			}
 		}
 		
 		return null;
 		
 	}
 	
}
