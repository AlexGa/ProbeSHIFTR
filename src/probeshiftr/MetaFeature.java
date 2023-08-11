package probeshiftr;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.IntStream;


public class MetaFeature extends Feature{
	
	ArrayList<Feature> exons = null;
	ArrayList<Feature> cds = null;
	ArrayList<Feature> introns = null;
	
	public MetaFeature() {
		super();
		// TODO Auto-generated constructor stub
	}
	public MetaFeature(String id, String chromosome, String source,
			String strand, int start, int end, String score,
			String frame, String description) {
		super(id, chromosome, source, strand, start, end, score, frame, description);
		// TODO Auto-generated constructor stub
	}
	public MetaFeature(String id, String chromosome, String source,
			String strand, int start, int end) {
		super(id, chromosome, source, strand, start, end);
		// TODO Auto-generated constructor stub
	}
	public MetaFeature(String id, String chromosome, String source,
			String strand, int start, int end, String score,
			String frame, String description, ArrayList<Feature> exons, ArrayList<Feature> cds, ArrayList<Feature> introns) {
		super(id, chromosome, source, strand, start, end, score, frame, description);
		this.exons = exons;
		this.cds = cds;
		this.introns = introns;
	}
	public MetaFeature(String id, String chromosome, String source,
			String strand, int start, int end, String score,
			String frame, String description, ArrayList<Feature> exons, ArrayList<Feature> introns) {
		super(id, chromosome, source, strand, start, end, score, frame, description);
		this.exons = exons;
		this.introns = introns;
	}
	public MetaFeature(String id, String chromosome, String source,
			String strand, int start, int end, String score,
			String frame, String description, ArrayList<Feature> exons) {
		super(id, chromosome, source, strand, start, end, score, frame, description);
		this.exons = exons;
	}
	
	public MetaFeature(MetaFeature metafeature) {
		super(metafeature);
		this.exons = metafeature.getExons();
		this.cds =  metafeature.getCds();
		this.introns =  metafeature.getIntrons();
	}
	
	public void addExon(Feature exon) {
	
		if(this.exons == null) {
			
			this.exons = new ArrayList<Feature>();
		}
		
		this.exons.add(exon);
	}
	
	public void addIntron(Feature intron) {
		
		if(this.introns == null) {
			
			this.introns = new ArrayList<Feature>();
		}
		
		this.introns.add(intron);
	}
	
	public void addCds(Feature cds) {
		
		if(this.cds == null) {
			
			this.cds = new ArrayList<Feature>();
		}
		
		this.cds.add(cds);
	}
	
	public void calculateIntrons(){
		try{
			
			if(this.getExons().size() == 1) {
				
				this.introns = null;
				
			}else {
				
				this.introns = new ArrayList<Feature>();
				Feature[] exonArray = new Feature[this.exons.size()];
				this.exons.toArray(exonArray);
				//this.sort(exonArray);
				Arrays.sort(exonArray);
				
				for(int i=0; i<exonArray.length - 1; i++){
					Feature intron = new Feature(exonArray[i]);
					
					/*intron.setStart(exonArray[i].getEnd()+1);
					intron.setEnd(exonArray[i+1].getStart());*/
					
					Interval2d intronInterval = new Interval2d(exonArray[i].getInterval().getEnd()+1, exonArray[i+1].getInterval().getStart()-1);
					intron.setInterval(intronInterval);
					this.introns.add(i, intron);
					
					intron.setDescription(intron.getDescription().replace("exon", "intron"));
					HashMap<String, String> descHash= intron.description2Hash();
					
					
					if(descHash.containsKey("intron_number")) {
						
						descHash.put("intron_number", String.valueOf(Integer.parseInt(descHash.get("intron_number")) - 1));
						
					}
					

					if(descHash.containsKey("intron_id")) {
						
						String[] oldId = descHash.get("intron_id").split("intron");
						
						if(oldId.length == 2) {
							if(!oldId[1].matches("[0-9]+")) {
								
								Pattern p = Pattern.compile("\\d+");
								Matcher m = p.matcher(oldId[1]);
								
								if(m.find()) {
								
									int number = Integer.parseInt(m.group()) - 1;
									oldId[1].replace(m.group(), String.valueOf(number));
							    
								}else {
							    	
									oldId[1] += "(-1)";
							    
								}
								
							}else {
							
								oldId[1] = String.valueOf(Integer.parseInt(oldId[1]) - 1);
							}
							
							descHash.put("intron_id", oldId[0] + "intron" + oldId[1]);	
							
						}else {

							descHash.put("intron_id", oldId[0] + ".intron" + descHash.get("intron_number"));
						}
						
						
					}
					
					intron.Hash2Description(descHash);
				}
			
			}
			
		}catch(Exception e){

			System.err.println("Error at intron calculation: " + this.getId());
			System.err.println(this.toString());
			e.printStackTrace(new java.io.PrintStream(System.err));
			System.err.println();
			
		}
	}
	/*public void sort(Feature[] array){
		// Using mergesort to put features in the right order based on their start coordinates
		
		if(array.length > 1){
			int mid = array.length/2;
			Feature[] leftArr = new Feature[mid];
			for(int i=0; i<mid; i++){
				leftArr[i] = array[i];
			}
			Feature[] rightArr = new Feature[array.length-mid];
			for(int i=mid, j=0; j<(array.length-mid); i++, j++){
				rightArr[j] = array[i];
			}
			sort(leftArr);
			sort(rightArr);
			
			int i = 0;
			int j = 0;
			int k = 0;
			
			while((i < leftArr.length) & (j < rightArr.length)){
				if(leftArr[i].getStart() < rightArr[j].getStart()){
					array[k] = leftArr[i];
					i++;
				}else{
					array[k] = rightArr[j];
					j++;
				}
				k++;
				while(i < leftArr.length){
					array[k] = leftArr[i];
					i++;
					k++;
				}
				while(j < rightArr.length){
					array[k] = rightArr[j];
					j++;
					k++;
				}
			}
		}
		
	}*/
	public ArrayList<Feature> getExons() {
		return exons;
	}
	public void setExons(ArrayList<Feature> exons) {
		
		this.exons = exons;
	}
	public ArrayList<Feature> getCds() {
		return cds;
	}
	public void setCds(ArrayList<Feature> cds) {
		this.cds = cds;
	}
	
	public void setCdsFrames() {
		
		if(this.cds != null) {
			
			Collections.sort(this.cds);
			int frame = 0;
			
			if(this.getStrand().equals("-")) {
			
				Feature cds_item = this.cds.get(this.cds.size() - 1);
				cds_item.setFrame(String.valueOf(frame));
				
				this.cds.set(this.cds.size() - 1, cds_item);
				
				for(int i = this.cds.size()-1; i >= 1; i--) {
					
					frame = (3 - this.cds.get(i).getWidth() % 3) % 3;
					cds_item = this.cds.get(i-1);
					cds_item.setFrame(String.valueOf(frame));
					this.cds.set(i-1, cds_item);
				}
				
			}else if(this.getStrand().equals("+")){
				
				Feature cds_item = this.cds.get(0);
				cds_item.setFrame(String.valueOf(frame));
				
				this.cds.set(0, cds_item);
				

				for(int i = 0; i < this.cds.size() - 1; i++) {
					
					frame = (3 - (this.cds.get(i).getWidth() - frame) % 3) % 3;
					cds_item = this.cds.get(i+1);
					cds_item.setFrame(String.valueOf(frame));
					this.cds.set(i+1, cds_item);
				}
				
			}
		}
	}
	
	public ArrayList<Feature> getIntrons() {
		return introns;
	}
	public void setIntrons(ArrayList<Feature> introns) {
		this.introns = introns;
	}	
	
	public String printGTFRows(){
		
		String rows = null;
		rows = this.getChromosome()+"\t"+this.getSource()+"\t"+"transcript"+"\t"+this.getStart()+"\t"+this.getEnd()+"\t"+this.getScore()+"\t"+this.getStrand()+"\t"+this.getFrame()+"\t"+this.getDescription()+"\n";
		
		
		if(this.getExons() != null){
			ArrayList<Feature> exons = this.getExons();
			Collections.sort(exons);
			Iterator<Feature> exonIter = exons.iterator();
			while(exonIter.hasNext()){
				Feature exon = exonIter.next();
				rows += exon.getChromosome()+"\t"+this.getSource()+"\t"+"exon"+"\t"+exon.getStart()+"\t"+exon.getEnd()+"\t"+exon.getScore()+"\t"+exon.getStrand()+"\t"+exon.getFrame()+"\t"+exon.getDescription()+"\n";
			}
		}else{
			System.err.println("Transcript "+ this.getId() + " does not have exons!");
		}
		
		if(this.getCds() != null){
			ArrayList<Feature> cdses = this.getCds();
			Collections.sort(cdses);
			Iterator<Feature> cdsIter = cdses.iterator();
			while(cdsIter.hasNext()){
				Feature cds = cdsIter.next();
				rows += cds.getChromosome()+"\t"+this.getSource()+"\t"+"CDS"+"\t"+cds.getStart()+"\t"+cds.getEnd()+"\t"+cds.getScore()+"\t"+cds.getStrand()+"\t"+cds.getFrame()+"\t"+cds.getDescription()+"\n";
			}
		}
		return rows;
	}
	
	public String printSpliceAnnotation() {
	
		String line = null;
		HashMap<String, String> descHash = this.description2Hash();
		
		line = String.join("\t", descHash.get("gene_id"), this.getId(), this.getChromosome(), this.getStrand(), String.valueOf(this.getStart()),  String.valueOf(this.getEnd()));
		
		if(this.getCds() != null & this.get_cds_length() > 0) {
			
			line += "\t" + String.join("\t", String.valueOf(this.getCds().get(0).getStart()), String.valueOf(this.getCds().get(this.getCds().size()-1).getEnd()));
		
		}else {
			
			line += "\t" + String.join("\t", ".", ".");
			
		}
		
		line += "\t" + String.valueOf(this.getExons().size());
		
		
		String exonStarts = String.valueOf(this.getExons().get(0).getStart()), 
	  		   exonEnds = String.valueOf(this.getExons().get(0).getEnd());
			
		if(this.getExons().size() == 1) {
				
			line += "\t" + String.join("\t", exonStarts, exonEnds, ".", ".", ".");
				
		}else {
				
			for(int i = 1; i < this.getExons().size(); i++) {
					
				exonStarts += "," + String.valueOf(this.getExons().get(i).getStart());
				exonEnds += "," + String.valueOf(this.getExons().get(i).getEnd());
			}
				
			String intronCov = String.valueOf(this.getIntrons().get(0).getCoverage());
				
			for(int i = 1; i < this.getIntrons().size(); i++){
					
				intronCov += "," + this.getIntrons().get(i).getCoverage();
			}
			
			line += "\t" + String.join("\t", exonStarts, exonEnds, intronCov);
		}
		
 		return line;
	}
	
	/** Compare exon structures
	 * 
	 * ! Equality is also reached if the 5' or 3' UTR containing exon (first and last exon) is shorter compared to reference
	 * 
	 * */
	public boolean equals2(MetaFeature mf){
		
		boolean isEqual = false;
		
		if((this.exons != null) && (mf.getExons().size() == this.exons.size())){
			
			if(this.exons.size() == 1){
				/* 
				 * if different mono exonic transcripts exist they might only differ in their 3' or 5' ends
				 * thus they are in this context equal because they share the same CDS sequence and different length in start or end
				 * can be due to uncertainties of the transcriptional machinery
				 * */
				isEqual = true; //this.exons.get(0).equals(mf.exons.get(0));
				
			}else{
				
				Collections.sort(this.exons);
				Collections.sort(mf.exons);

				isEqual = this.exons.get(0).getEnd() == mf.exons.get(0).getEnd();
				
				if(this.exons.size() == 2){
					isEqual &= this.exons.get(1).getStart() == mf.exons.get(1).getStart();
				}else{
					for(int i = 1; i < this.exons.size() - 1; i++){
						//isEqual &= this.exons.get(i).equals(mf.exons.get(i));
						isEqual &= this.exons.get(i).getInterval().equals(mf.exons.get(i).getInterval());
						/*isEqual &= this.exons.get(i).getStart() == mf.exons.get(i).getStart();
						isEqual &= this.exons.get(i).getEnd() == mf.exons.get(i).getEnd();*/
					}
					
					isEqual &= this.exons.get(this.exons.size() - 1).getStart() == mf.exons.get(this.exons.size() - 1).getStart() ;
				}
			}
		}
		return isEqual;
	}
	
	public boolean equals(MetaFeature mf){
		
		boolean isEqual = false;
		
		if((this.exons != null) && (mf.getExons().size() == this.exons.size())){
			
			if(this.exons.size() == 1){
				/* 
				 * if different mono exonic transcripts exist they might only differ in their 3' or 5' ends
				 * thus they are in this context equal because they share the same CDS sequence and different length in start or end
				 * can be due to uncertainties of the transcriptional machinery
				 * */
				isEqual = true; //this.exons.get(0).equals(mf.exons.get(0));
				
			}else{
			
				this.calculateIntrons();
				mf.calculateIntrons();
				
				ArrayList<Feature> introns = mf.getIntrons();
				
				Collections.sort(this.introns);
				Collections.sort(introns);
				
				
				isEqual = this.introns.get(0).equals(introns.get(0));
				
				for(int i=1; i<introns.size(); i++){
					
					isEqual &= this.introns.get(i).equals(introns.get(i));
				}
				
//				if(this.getId().equals("AT5G38100.2")){
//					
//					System.out.println(isEqual);
//					for(int i=0; i<introns.size(); i++){
//						
//						System.out.println("Known: " + this.getIntrons().get(i).toString());
//						System.out.println("New: " + introns.get(i).toString());
//					}
//				}
				
			}
		}
		return isEqual;
	}

	public ArrayList<Feature> circExonsFromInterval(Interval2d circInt){
		
		Collections.sort(this.exons);
		ArrayList<Feature> new_exons = null;
		
		
		Iterator<Feature> exonIter = this.exons.iterator();
		Feature curExon = null;
		
		while(exonIter.hasNext()){
			
			curExon = exonIter.next();
			if(new_exons == null){
				
				if(curExon.getEnd() > circInt.getStart()){
					
					if(curExon.getStart() > circInt.getEnd()) {
						break;
					}
					
					new_exons = new ArrayList<Feature>();
					Feature update_exon = new Feature(curExon);

					/*
					 * circRNA starts in an intergenic region -> first exon is host transcript exon
					 * othwise -> exon start shifted to circRNA start
					*/
					if(circInt.getStart() >= curExon.getStart())
						update_exon.setStart(circInt.getStart());
			
					
					if(curExon.getEnd() >= circInt.getEnd()){
						
						update_exon.setEnd(circInt.getEnd());
						new_exons.add(update_exon);
						break;
						
					}else{
						new_exons.add(update_exon);
					}
					
				}
				
			}else{
				
				Feature update_exon = new Feature(curExon);
				
				if(curExon.getEnd() >= circInt.getEnd()){
					
					// circRNA ends in intron region
					if(curExon.getStart() > circInt.getEnd()) {
						break;
					}else {
						update_exon.setEnd(circInt.getEnd());
						new_exons.add(update_exon);
						break;
					}
					
					
					
				}else{
					
					new_exons.add(update_exon);
				}
			}	
		}
	
		return new_exons;
	}
	public ArrayList<Feature> interval2dToCds(Interval2d orfInt){
		
		Collections.sort(this.exons);
		ArrayList<Feature> new_cds = new ArrayList<Feature>();
		
		//this.setCds(null);
		
		Iterator<Feature> exonIter = this.exons.iterator();
		Feature curExon = null;
		
		while(exonIter.hasNext()){
			
			curExon = exonIter.next();
			if(new_cds.isEmpty()){
				
				if(curExon.getEnd() > orfInt.getStart()){
					
					Feature cds_feat = new Feature(curExon);
					cds_feat.setStart(orfInt.getStart());
					
					if(curExon.getEnd() >= orfInt.getEnd()){
						
						cds_feat.setEnd(orfInt.getEnd());
						new_cds.add(cds_feat);
						break;
						
					}else{
						new_cds.add(cds_feat);
					}
					
				}
				
			}else{
				
				Feature cds_feat = new Feature(curExon);
				
				if(curExon.getEnd() >= orfInt.getEnd()){
					
					cds_feat.setEnd(orfInt.getEnd());
					new_cds.add(cds_feat);
					break;
					
				}else{
					
					new_cds.add(cds_feat);
				}
			}	
		}
		
		this.setCds(new_cds);
		this.setCdsFrames();
		return cds;
	}	
	
	public int get_cds_length(){
		
		int cds_length = 0;
		
		if(this.cds != null){
			if(this.getStrand().equals("+")){
				for(int i=0; i<this.cds.size(); i++)
					cds_length +=  this.cds.get(i).getInterval().length();
//					cds_length += this.cds.get(i).getEnd() - this.cds.get(i).getStart() + 1;
			}else{
				for(int i=(this.cds.size()-1); i>=0; i--)
					cds_length +=  this.cds.get(i).getInterval().length();
//					cds_length += this.cds.get(i).getEnd() - this.cds.get(i).getStart() + 1;
			}
			
		}
		
		return cds_length;
		
	}
	
	public int get_transcript_length(){
		
		int transcript_length = 0;

		for(int i=0; i<this.exons.size(); i++)
			transcript_length += this.exons.get(i).getInterval().length();//this.exons.get(i).getEnd() - this.exons.get(i).getStart() + 1;
	
		return transcript_length;
		
	}
	
	public void updateExons(ArrayList<Feature> exons){
		
		/* Check and update transcript boundaries*/
		
		Collections.sort(this.exons);
		Collections.sort(exons);
		
		String addInfo = "";
		if(this.exons.size() > 1){
			
			if(exons.get(0).getStart() < this.exons.get(0).getStart()){
			
				this.setStart(exons.get(0).getStart());
			
				Feature startExon = this.exons.get(0);
				startExon.setStart(exons.get(0).getStart());
			
				this.exons.remove(0);
				this.exons.add(0, startExon);
			
				if(this.getStrand().equals("+")){
					addInfo = "alt_5utr";
				}
				if(this.getStrand().equals("-")){
					addInfo = "alt_3utr";
				}
			}
		
			if(this.exons.get(this.exons.size()-1).getEnd() < exons.get(exons.size() - 1).getEnd()){

				this.setEnd(exons.get(exons.size()-1).getEnd());
				
				Feature endExon = this.exons.get(this.exons.size()-1);
				endExon.setEnd(exons.get(exons.size()-1).getEnd());
			
				this.exons.remove(this.exons.size() -1);
				this.exons.add(this.exons.size(), endExon);
				
				if(this.getStrand().equals("+")){
					
					addInfo += addInfo.equals("") ? "alt_3utr" : "_alt_3utr";
				}
				if(this.getStrand().equals("-")){
					
					addInfo = addInfo.equals("") ? "alt_5utr" : "alt_5utr_"+addInfo;
				}
				
			}

		}else{
			if((exons.get(0).getStart() < this.exons.get(0).getStart()) & (this.exons.get(0).getEnd() < exons.get(0).getEnd())){
				
				this.setStart(exons.get(0).getStart());
				this.setEnd(exons.get(0).getEnd());
				
				Feature exon = this.exons.get(0);
				
				exon.setStart(exons.get(0).getStart());
				exon.setEnd(exons.get(0).getEnd());
				
				this.exons.remove(0);
				this.exons.add(0, exon);
				
				addInfo = "alt_5utr_alt_3utr";
				
			}else if(exons.get(0).getStart() < this.exons.get(0).getStart()){
				
				this.setStart(exons.get(0).getStart());
				
				Feature exon = this.exons.get(0);
				
				exon.setStart(exons.get(0).getStart());
				
				this.exons.remove(0);
				this.exons.add(0, exon);
				
				if(this.getStrand().equals("+")){
					addInfo = "alt_5utr";
				}
				if(this.getStrand().equals("-")){
					addInfo = "alt_3utr";
				}
				
			}else if(this.exons.get(0).getEnd() < exons.get(0).getEnd()){
				
				this.setEnd(exons.get(0).getEnd());
				
				Feature exon = this.exons.get(0);
				
				exon.setEnd(exons.get(0).getEnd());
				
				this.exons.remove(0);
				this.exons.add(0, exon);
				
				if(this.getStrand().equals("+")){
					addInfo = "alt_3utr";
				}
				if(this.getStrand().equals("-")){
					addInfo = "alt_5utr";
				}
			}
		}
		
		if(addInfo != ""){
			HashMap<String, String> dHash = this.description2Hash();
			dHash.put("info", addInfo);
			this.Hash2Description(dHash);
		}
	}
	
	public void updateSubFeatureDescriptions(String key, String value){
		
		if(this.exons != null){
			
			ArrayList<Feature> exList = new ArrayList<Feature>();
			
			for(Feature ex: this.getExons()){
				
				ex.updateDescription(key, value);
				exList.add(ex);
			}
			this.setExons(exList);
		}
		
		if(this.cds != null){
			
			ArrayList<Feature> cdsList = new ArrayList<Feature>();
			
			for(Feature cds: this.getCds()){
				
				cds.updateDescription(key, value);
				cdsList.add(cds);
			}
			this.setCds(cdsList);
		}

	}
	
	public int whichExon(Interval2d exonInterval) {
		
		int[] idx = IntStream.range(0, this.exons.size()).toArray();
		int index = Arrays.stream(idx).boxed().filter(i -> this.exons.get(i).getInterval().equals(exonInterval)).findFirst().orElse(-1);
		
		return index;
	}
	
	public int whichIntron(Interval2d intronInterval) {
		
		int[] idx = IntStream.range(0, this.introns.size()).toArray();
		int index = Arrays.stream(idx).boxed().filter(i -> this.introns.get(i).getInterval().equals(intronInterval)).findFirst().orElse(-1);
		
		return index;
	}

	public GenomicLocation location(int position) {
		
		Interval2d[] intervals = new Interval2d[2 * this.getExons().size() - 1];
		
		int j = 0;
		for(int i = 0; i < intervals.length; i++) {
			if( i % 2 == 0)
			
				intervals[i] = this.getExons().get(j).getInterval();
			
			else {
				
				if(this.getIntrons() == null) {
				
					this.calculateIntrons();
				
				}
				intervals[i] = this.getIntrons().get(j).getInterval();
				j++;
			}	
		}
		
		int ind = Interval2d.binarySearchOverlap(intervals, position);
		GenomicLocation posLocation = new GenomicLocation(ind >= 0, (ind % 2 == 0), ind/2);
		
		return posLocation;
	}
	
	public void CheckInterval() {
		
		Interval2d longestInterval = this.getLongestInterval();
		
		if(!longestInterval.equals(this.getInterval())) {
			this.setInterval(longestInterval);
		}
		
	}
	
	private Interval2d getLongestInterval() {
		
		int minStart = Integer.MAX_VALUE;
		int maxEnd = Integer.MIN_VALUE;
		for(Feature exon: this.getExons()) {
			
			if(minStart > exon.getStart()) {
				minStart = exon.getStart();
			}
			
			if(maxEnd < exon.getEnd()) {
				maxEnd = exon.getEnd();
			}
		}
		
		
		return new Interval2d(minStart, maxEnd);
	}
}
