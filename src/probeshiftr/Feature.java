package probeshiftr;

import java.util.HashMap;


public class Feature implements Comparable<Feature>{

	private String id;
	private String chromosome;
	private String source = null;
	private String strand;
	/*private int start;
	private int end;*/
	private Interval2d interval;
	private String score;
	private String frame;
	private String description; 
	private int coverage;
	
	public Feature() {
		super();
	}
	

	public Feature(String id, String chromosome, String source, String strand,
			int start, int end, String score, String frame, String description, int coverage) {
	
		this.id = id;
		this.chromosome = chromosome;
		this.source = source;
		this.strand = strand;
		this.interval = new Interval2d(start, end);
		this.score = score;
		this.frame = frame;
		this.description = description;
		this.coverage = coverage;
		
	}
	
	public Feature(Feature feature) {
		
		this(feature.getId(),  feature.getChromosome(), feature.getSource(), 
			 feature.getStrand(), feature.getInterval().getStart(), feature.getInterval().getEnd(),
			 feature.getScore(), feature.getFrame(),feature.getDescription(),feature.coverage);

	}
	public Feature(String id, String chromosome, String source, String strand,
			int start, int end) {
		super();
		this.id = id;
		this.chromosome = chromosome;
		this.source = source;
		this.strand = strand;
		/*this.start = start;
		this.end = end;*/
		this.interval = new Interval2d(start, end);
		this.score = ".";
		this.frame = ".";
		this.coverage = 0;
	}

	public Feature(String id, String chromosome, String source, String strand,
			int start, int end, String score, String frame, String description) {
		super();
		this.id = id;
		this.chromosome = chromosome;
		this.source = source;
		this.strand = strand;
		/*this.start = start;
		this.end = end;*/
		this.interval = new Interval2d(start, end);
		this.score = score;
		this.frame = frame;
		this.description = description;
		this.coverage = 0;
	}
	
	public Feature(String id, String chromosome, String strand,
			int start, int end, String score, String offset, String description) {
		this.id = id;
		this.chromosome = chromosome;
		this.strand = strand;
		/*this.start = start;
		this.end = end;*/
		this.interval = new Interval2d(start, end);
		this.score = score;
		this.description = description;
	}

	public Feature(String id, String chromosome, String strand,
			int start, int end) {
		this.id = id;
		this.chromosome = chromosome;
		this.strand = strand;
		/*this.start = start;
		this.end = end;*/
		this.interval = new Interval2d(start, end);
		this.score = ".";
		this.frame = ".";
	}
	
	public Feature(String id, String chromosome, String strand,
			int start, int end, int coverage) {
		this.id = id;
		this.chromosome = chromosome;
		this.strand = strand;
		/*this.start = start;
		this.end = end;*/
		this.interval = new Interval2d(start, end);
		this.score = ".";
		this.frame = ".";
		this.coverage = coverage;
		this.description = "";
	}
	
	public String getChromosome() {
		return chromosome;
	}
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	public String getSource() {
		return source;
	}
	public void setSource(String source) {
		this.source = source;
	}
	public String getStrand() {
		return strand;
	}
	public void setStrand(String strand) {
		this.strand = strand;
	}
	public int getStart() {
		return this.interval.getStart();
	}
	public int getEnd() {
		return this.interval.getEnd();
	}
	public void setStart(int start) {
		this.interval.setStart(start);
	}
	public void setEnd(int end) {
		this.interval.setEnd(end);
	}
	public String getScore() {
		return score;
	}
	public void setScore(String score) {
		this.score = score;
	}
	public String getFrame() {
		return frame;
	}
	public void setFrame(String frame) {
		this.frame = frame;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getDescription() {
		return description;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	public int getCoverage() {
		return coverage;
	}

	public void setCoverage(int coverage) {
		this.coverage = coverage;
	}

	public Interval2d getInterval() {
		return interval;
	}

	public void setInterval(Interval2d interval) {
		this.interval = interval;
	}

	@Override
	public String toString() {
		return "Feature [id=" + id + ", chromosome=" + chromosome + ", strand="
				+ strand + ", start=" + interval.getStart() + ", end=" + interval.getEnd() + ", score="
				+ score + ", frame=" + frame + ", description=" + description
				+ "]";
	}
	
	@Override
	public int compareTo(Feature f) {
		
		int comp = this.chromosome.compareTo(f.chromosome);
		
		if(comp == 0){
			comp = this.strand.compareTo(f.strand);
			if(comp == 0){
				comp = Integer.compare(this.interval.getStart(), f.interval.getStart());
			}
			if(comp == 0){
				comp = Integer.compare(this.interval.getEnd(), f.interval.getEnd());
			}
		}
		return comp;
	}
	
	public int compareTo(int position) {
		
		int comp =  this.interval.compareTo(position);
		
		return comp;
	}
	
	public int overlapLength(Feature f){
		
		int overlap_length = 0;
		
		int comp = this.chromosome.compareTo(f.chromosome);
		
		if(comp == 0){
			comp = this.strand.compareTo(f.strand);
			if(comp == 0){
				
				overlap_length = this.interval.overlap(f.getInterval());
//				if((f.start > this.end) || (this.start > f.end)){
//					overlap_length = 0;
//				}else{
//					int ov_start = Math.max(this.start, f.start);
//					int ov_end = Math.min(this.end, f.end);
//					overlap_length = ov_end - ov_start;
//				}
			}
		}
		return overlap_length;
		
	}
	
	public boolean overlap(Feature f) {
		
		boolean overlap = false;

		overlap = this.overlapLength(f) > 0;
		
		return overlap;
	}
	
	public boolean equals(Feature f){
		
        boolean isEqual = false;

        if (f != null && f instanceof Feature)
        {
        	isEqual = this.compareTo(f) == 0;
        }

        return isEqual;
    }

	protected HashMap<String, String> description2Hash(){
		
		HashMap<String, String> dHash = new HashMap<String, String>();
		
		String[] desc =  this.description.split("; ");
		
		for(String entry: desc){
			
			String[] keyVal = entry.split(" \"");
			dHash.put(keyVal[0], keyVal[1].replace("\"", "").replace(";", ""));
		}
		
		return dHash;
	}
	
	protected void Hash2Description(HashMap<String, String> dHash){
		
		StringBuilder desc = new StringBuilder();
		
		for(String key: dHash.keySet()){
			
			desc.append(key + " \"" + dHash.get(key) + "\"; ");
		}
		
		this.setDescription(desc.toString());
	}
	public void updateDescription(String key, String value){
		
		HashMap<String, String> descHash = this.description2Hash();
		descHash.put(key, value);
		this.Hash2Description(descHash);
		
	}
	
	public int getWidth(){
		
		return this.interval.length();//this.end - this.start + 1;
	}

}
