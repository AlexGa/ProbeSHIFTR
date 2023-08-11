package probeshiftr;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class SeqUtils {

	SeqUtils(){
		
	}
	
	public HashMap<String,String> extractTranscripts(HashMap<String,String> genome, HashMap<String,ArrayList<MetaFeature>> coords){
		
		HashMap<String,String> transcripts = new HashMap<String,String>();
		
		/*not the best/fastest solution -> run over the hashMap and cut the sequences out*/
		
		Iterator<String> transIter = coords.keySet().iterator();
		
		while(transIter.hasNext()){
			
			String transcriptID = transIter.next();
			
			ArrayList<MetaFeature> transCoords = coords.get(transcriptID);
			String transSeq = "";
			for(int i=0; i<transCoords.size(); i++){
				
				String chr = transCoords.get(i).getChromosome();
				int[] begEnd = {transCoords.get(i).getStart(), transCoords.get(i).getEnd()};
			
				String tmpSeq = "";
				
				if(!genome.containsKey(chr)){
					continue;
				}
				
				tmpSeq = genome.get(chr).substring(begEnd[0]-1, begEnd[1]);
				
				if(transCoords.get(i).getStrand().equals("-")){
	
					tmpSeq = tmpSeq.replaceAll("A", "t");
					tmpSeq = tmpSeq.replaceAll("C", "g");
					tmpSeq = tmpSeq.replaceAll("G", "c");
					tmpSeq = tmpSeq.replaceAll("T", "a");
					tmpSeq = new StringBuffer(tmpSeq.toUpperCase()).reverse().toString();
					
				}
				transSeq += tmpSeq;
			}
			
			transcripts.put(transcriptID, transSeq);
		}
		return transcripts;
		
	}
	
	public static Sequence getConcatFeatureSequence(Sequence seq, ArrayList<Feature> f_blocks){
		
		Sequence final_seq = null;
		
		if(f_blocks != null && f_blocks.size() >= 1){

			Collections.sort(f_blocks);
			StringBuilder seq_data = new StringBuilder();
			
			for(int i=0; i<f_blocks.size(); i++){
				
				try{
					
					String shortSeq = seq.getSeq().substring(f_blocks.get(i).getStart()-1, f_blocks.get(i).getEnd());
					seq_data.append(shortSeq);
					
				}catch(StringIndexOutOfBoundsException e){
					
					e.printStackTrace();
					System.err.println(f_blocks.get(i).toString());
					System.err.println(f_blocks.get(i).getStart()-1 + " - " + f_blocks.get(i).getEnd());
					System.exit(0);
				
				}
				
				/*
				if(f_blocks.get(i).getId().equals("NewID_00017987")){
					System.out.println(shortSeq.length() + " (" + shortSeq.length()%3+ ")");
				}
				*/
			}
			
			final_seq = new Sequence(f_blocks.get(0).getId(), seq_data.toString());
			
			if(f_blocks.get(0).getStrand().equals("-")){
				
				final_seq = final_seq.getReverseComplement();
				
			}
			
			return final_seq;
		}else{
			throw new IllegalArgumentException("Please check your sequence and/or feature coordinates!");
		}
		
	}
	
	private HashMap<Integer, ArrayList<Integer>> getMatches(String regexp, String seq_string){
		
		HashMap<Integer, ArrayList<Integer>> matchIndices = new HashMap<Integer, ArrayList<Integer>>();		
        Matcher m = Pattern.compile(regexp, Pattern.CASE_INSENSITIVE).matcher(seq_string);

        while(m.find()) {
        	
        	int startIndex = m.start();
        	int classIndex = startIndex % 3;
        	
        	if(matchIndices.containsKey(classIndex)){
        		
        		matchIndices.get(classIndex).add(startIndex);
        		
        	}else{
        		
        		ArrayList<Integer> classList = new ArrayList<Integer>();
        		classList.add(startIndex);
        		matchIndices.put(classIndex, classList);
        	
        	}
        	
        }

		return matchIndices;
		
	}
	
	private Interval2d convertORFcoordinates(Interval2d orfInt, ArrayList<Feature> exons){
		
//		if(exons.get(0).getId().equals("NewID_00017987")){//"NewID_00153029")){//NewID_00165038")){
//			System.out.println("before convert: "+ orfInt.toString());
//			/*for(Feature exon: exons){
//				System.out.println(exon.toString());
//				System.out.println(exon.getEnd() - exon.getStart() + 1);
//			}*/
//		}
//		
		int start = orfInt.getStart();
		int end = orfInt.getEnd();
		
		Collections.sort(exons);
		
		for(Feature exon: exons){
			
			int exon_length = exon.getEnd() - exon.getStart() + 1;
			
			if(exon_length <= start){
				
				start -= exon_length;
			
			}else{
			
				start += exon.getStart();
				
				if(start == exon.getEnd()){
				
					start = 2;
					end += 2;
					continue;
				
				}
				break;
			
			}
		}
		
		for(Feature exon: exons){
			
			int exon_length = exon.getEnd() - exon.getStart() + 1;

			if(exon_length < end){
				end -= exon_length;
				
			}else{
				end = exon.getStart() + end - 1;

				if(end > exon.getEnd()){
					System.err.println("Wrong end! (too big)");
				}

				break;
			}
		}

		return new Interval2d(start, end);
	}
	
	public HashMap<Integer, ArrayList<Interval2d>> predictPossibleOrfCoordinates(Sequence seq, ArrayList<Feature> f_blocks){
		
		Collections.sort(f_blocks);
		
		Sequence tx_seq = getConcatFeatureSequence(seq, f_blocks);
		
		String start_exp = "atg";
		String stop_exp = "taa|tag|tga";
		
		HashMap<Integer, ArrayList<Integer>> startIndices = getMatches(start_exp, tx_seq.getSeq());
		HashMap<Integer, ArrayList<Integer>> stopIndices = getMatches(stop_exp, tx_seq.getSeq());
		
		Iterator<Integer> outIter = startIndices.keySet().iterator();
		HashMap<Integer, ArrayList<Interval2d>> frameOrfPairs = new HashMap<>();
		
		while(outIter.hasNext()){
			
			int frame = outIter.next();
			
			if(!startIndices.containsKey(frame)){
				
				continue;
			}
			
			if(!stopIndices.containsKey(frame)){

				continue;
			}
			
			int newFrame = -1;
			
			HashMap<Integer, Integer> uniqueStartStop  = new HashMap<>();
			ArrayList<Integer> startList = startIndices.get(frame);
			ArrayList<Integer> stopList = stopIndices.get(frame);
			
			// Sort start and stop codons
			Collections.sort(startList);
			Collections.sort(stopList);

			int i=0,j=0; 
			int lastStop = 0;
			
			/* ORFs are defined as non overlapping intervals */
			
			while(i<startList.size() & j<stopList.size()){
				
				if(startIndices.get(frame).get(i) < lastStop) {
					i++;
				}else if((startIndices.get(frame).get(i) > stopIndices.get(frame).get(j))){
					j++;
				}else{
					uniqueStartStop.put(startIndices.get(frame).get(i), stopIndices.get(frame).get(j));
					lastStop = stopIndices.get(frame).get(j);
					i++;
					j++;
				}
			}
			
			Iterator<Integer> startIter = uniqueStartStop.keySet().iterator();
			
			while(startIter.hasNext()){
			
					int startCodon = startIter.next();
					int stopCodon = uniqueStartStop.get(startCodon);
					
					Interval2d orfInterval = new Interval2d();
					
					int extStop = stopCodon + 3;
					
					if(f_blocks.get(0).getStrand().equals("-")){
						
						orfInterval.setStart(tx_seq.getSeq().length() - extStop);
						orfInterval.setEnd(tx_seq.getSeq().length() - startCodon);
						
					}else{

						orfInterval.setStart(startCodon);
						orfInterval.setEnd(extStop);
					}
					
					orfInterval = convertORFcoordinates(orfInterval, f_blocks);
					
					if(f_blocks.get(0).getStrand().equals("-")){
					
						newFrame = orfInterval.getEnd() % 3;
					
					}else{
					
						newFrame = orfInterval.getStart() % 3;
					
					}
					if(frameOrfPairs.containsKey(newFrame)){
						
						ArrayList<Interval2d> orfPairList = frameOrfPairs.get(newFrame);
						orfPairList.add(orfInterval);
						frameOrfPairs.put(newFrame, orfPairList);
						
					}else{
						
						ArrayList<Interval2d> orfPairList = new ArrayList<Interval2d>();
						orfPairList.add(orfInterval);
						frameOrfPairs.put(newFrame, orfPairList);
					}
					
			}
		}

		return frameOrfPairs;	
	}
	
	public double txOverlap(Sequence refCds, Sequence newCds) {
		

		Alignment ali = new Alignment(refCds.getSeq(), newCds.getSeq());
		ali.globalAlignment();
		
		return Math.round(ali.identity()/(double) refCds.getSeq().length() * 10000.00d) / 100.00d;
		
	}
	
	public boolean validateIsoform(Gene gene, MetaFeature transcript, Sequence seq, int min_orf_length, double minOverlap){
		
		boolean belongsToGene = false;

		String strand = transcript.getStrand();

		if(gene.getStrand().equals(strand)){
			
			// first predict ORFs of transcript
			HashMap<Integer, ArrayList<Interval2d>> tx_orf_hash = predictPossibleOrfCoordinates(seq, transcript.getExons());
			
			int number_orfs = 0;
			ArrayList<Interval2d> orfList = new ArrayList<Interval2d>();
			
			for(ArrayList<Interval2d> orfs: tx_orf_hash.values()) {

				number_orfs += orfs.size();
				for(Interval2d orfInterval: orfs)
					orfList.add(orfInterval);
			}
			
			ArrayList<MetaFeature> isoformList = gene.getTranscripts();
			
			/*  
			 * 	- scoreMatrix contains for each combination of orf and isoform a score between -1, 1, 2
			 *  - score of -1 implies that a combination of isoforom and orf is not valid
			 */
			
			int[][] scoreMatrix = new int[isoformList.size()][number_orfs];
			double[][] overlapMatrix = new double[isoformList.size()][number_orfs];
			int[] orfLengths = new int[number_orfs]; 
			int[] isoLengths = new int[isoformList.size()]; 
			
			for(int i = 0; i < scoreMatrix.length; i++)
				Arrays.fill(scoreMatrix[i], 0);
			
			for(int i = 0; i < overlapMatrix.length; i++)
				Arrays.fill(overlapMatrix[i], 0.0d);
			
			int maxScore = -1;
			double maxOverlap = 0.0d;
			Integer[] overlapPair = new Integer[2];
			
			ArrayList<Integer[]> maxPairs = new ArrayList<Integer[]>();
			
			/* get cds lengths for each orf interval */
			
			for(int orfIndex = 0; orfIndex < number_orfs; orfIndex++) {
				
				Interval2d curOrf = orfList.get(orfIndex);
				transcript.interval2dToCds(curOrf);	
				orfLengths[orfIndex] = transcript.get_cds_length();
			}
			
			for(int isoIndex = 0; isoIndex < isoformList.size(); isoIndex++) {
				
				MetaFeature isoform = isoformList.get(isoIndex);
				
				// if reference isoform has no annotated CDS it will not be considered -> score = -1
				
				if(isoform.getCds() == null) {
				
					Arrays.fill(scoreMatrix[isoIndex], -1);
					
				}else{
					
					ArrayList<Feature> cdsList = isoform.getCds();
					isoLengths[isoIndex] = isoform.get_cds_length();
					
					// start and end of cds from reference isoform for following comparisons
					
					int isoStart = cdsList.get(0).getStart(); 
					int isoEnd = cdsList.get(cdsList.size() - 1).getEnd(); 
					
					Sequence refIsoformSequence = getConcatFeatureSequence(seq, cdsList);
					
					for(int orfIndex = 0; orfIndex < number_orfs; orfIndex++) {
					
						Interval2d curOrf = orfList.get(orfIndex);
						transcript.interval2dToCds(curOrf);		
						
						/* If cds length of ORF is less than minOverlap (default: 50%) than ORF of isoform the orf is considered not to be valid */
						
						double tx_iso_ratio = (double)transcript.get_cds_length()/(double) isoform.get_cds_length() * 100.0d;
						
						
						if(tx_iso_ratio < minOverlap) {
							
							scoreMatrix[isoIndex][orfIndex] = -2;
							
						}else {

							/* 1.) 
							 * 
							 * 	- Calculate overlap between cds of ref isoform and cds from predicted orfs of novel transcript 
							 *  - Check if overlap is greater than given cdsOverlap (default: 50%)
							 *  - if requirement is not fulfilled combination between isoform and orf is not valid -> score = -1
							 * */
							
						
							Sequence newIsoformSequence = getConcatFeatureSequence(seq, transcript.getCds());
							
							double overlap = txOverlap(refIsoformSequence, newIsoformSequence);
							
							overlapMatrix[isoIndex][orfIndex] = overlap;
							
							if(overlap >= minOverlap) {
								
								scoreMatrix[isoIndex][orfIndex] += 1;
								
								// store information of best combinations
								if(overlapMatrix[isoIndex][orfIndex] > maxOverlap) {
										
										maxOverlap = overlapMatrix[isoIndex][orfIndex];
										overlapPair[0] = isoIndex;
										overlapPair[1] = orfIndex;
								}
								
							}else {
								
								scoreMatrix[isoIndex][orfIndex] = -2;
							}
							
							/* 2.) 
							 * 
							 * 	- Check if the CDS from reference and the predicted CDS share the same start position
							 * 
							 * */
							
							boolean equalStart = (gene.getStrand().equals("+") && curOrf.getStart() == isoStart) || (gene.getStrand().equals("-") && curOrf.getEnd() == isoEnd);
							scoreMatrix[isoIndex][orfIndex] += equalStart ? 1 : -1;
							
							/* 3.) 
							 * 
							 * 	- Check if end from reference CDS is short 
							 *  - if the reference CDS is shorter than it could be an indication that the two cds are based on a frame shift
							 * 
							 * 
							
							boolean shorterEnd = (gene.getStrand().equals("+") && curOrf.getEnd() < isoEnd) || (gene.getStrand().equals("-") && curOrf.getStart() > isoStart);
							scoreMatrix[isoIndex][orfIndex] += shorterEnd ? -1 : 1;
							
							 */
						}

						//store information of best combinations
						if(scoreMatrix[isoIndex][orfIndex] >= maxScore) {
							
							if(scoreMatrix[isoIndex][orfIndex] > maxScore) {
								
								maxScore = scoreMatrix[isoIndex][orfIndex];
								maxPairs = new ArrayList<Integer[]>();
							}
							
							Integer[] pair = new Integer[2];
							pair[0] = isoIndex;
							pair[1] = orfIndex;
							maxPairs.add(pair);
						}
					}
				}
			}
			
			// Check if there is at least one combination that could fit
			if(maxScore > 0) {

				belongsToGene = true;
				
				if(maxPairs.size() == 1) {
					
					Integer[] pair = maxPairs.get(0);
					//System.out.println(orfList.get(pair[1]).toString());
					transcript.interval2dToCds(orfList.get(pair[1]));
					
				}else {
					
					transcript.interval2dToCds(orfList.get(overlapPair[1]));
				}
			}
			
			/*if(belongsToGene) {
				
				System.out.println("==== Isoform - Lengths ====");
				for(int i = 0; i < isoLengths.length; i++) {
					System.out.print(i + ".\t");
				}
				System.out.println();
				
				for(int i = 0; i < isoLengths.length; i++) 
					System.out.print(isoLengths[i] + "\t");
				System.out.println();
				
				System.out.println("==== ORF - Lengths ====");
				for(int i = 0; i < orfLengths.length; i++) {
					System.out.print(i + ".\t");
				}
				System.out.println();
				
				for(int i = 0; i < orfLengths.length; i++) 
					System.out.print(orfLengths[i] + "\t");
				System.out.println();
				
				
				System.out.println("==== Overlap - Matrix ====");
				for(int i = 0; i < overlapMatrix[0].length; i++) {
					System.out.print(i + ".\t");
				}
				System.out.println();
				
				for(int i = 0; i < overlapMatrix.length; i++) {
					for(int j = 0; j < overlapMatrix[i].length; j++)
						System.out.print(overlapMatrix[i][j] + "\t");
					System.out.println();
				}
				System.out.println();

				
				int k = 0;
				for(Interval2d i2d: orfList) {
					System.out.println((k) + ": " + i2d.toString() + " --> " + orfLengths[k]);
					k++;
				}
				
				System.out.println("==== Score - Matrix ====");
				for(int i = 0; i < scoreMatrix.length; i++) {
					for(int j = 0; j < scoreMatrix[i].length; j++)
						System.out.print(scoreMatrix[i][j] + "\t");
					System.out.println();
				}
				System.out.println();
			}*/
			
		}
		return belongsToGene;
		
	}
	
	public boolean validateIsoforms(Gene gene, MetaFeature transcript, Sequence seq, int min_orf_length){
		
		boolean belongsToGene = false;
		ArrayList<Feature> longestCds = null;
		
		// first predict ORFs of transcript
		HashMap<Integer, ArrayList<Interval2d>> tx_orf_hash = predictPossibleOrfCoordinates(seq, transcript.getExons());
		
		String strand = transcript.getStrand();

		System.out.println(transcript.getId());
		if(gene.getStrand().equals(strand)){
			
			// Check if and which isoform of the gene shares then same start codon
			ArrayList<MetaFeature> isoformList = gene.getTranscripts();
			
			HashMap<Integer, Integer[]> scoreMap = new HashMap<Integer, Integer[]>();

			Iterator<MetaFeature> isoIter = isoformList.iterator();
			int isoIndex = 0;
			
			int minScore = Integer.MAX_VALUE;
			
			while(isoIter.hasNext()){
				
				MetaFeature isoform = isoIter.next();
				
				if(isoform.getCds() != null){

					ArrayList<Feature> cdsList = isoform.getCds();
					Collections.sort(cdsList);
					
					int isoform_frame = 0;
					
					if(gene.getStrand().equals("-")){
					
						isoform_frame = cdsList.get(cdsList.size() - 1).getEnd() % 3;
					
					}else{
					
						isoform_frame = cdsList.get(0).getStart() % 3;
					
					}
					
					Set<Integer> orfKey_set = tx_orf_hash.keySet();
					
					double max_overlap = 0.0d;
					Sequence refIsoformSequence = getConcatFeatureSequence(seq, cdsList);
					
					for(Integer orfKey: orfKey_set) {
						
						for(Interval2d curOrf: tx_orf_hash.get(orfKey)) {
							
							transcript.interval2dToCds(curOrf);						
							Sequence newIsoformSequence = getConcatFeatureSequence(seq, transcript.getCds());
							
							max_overlap = Math.max(max_overlap, txOverlap(refIsoformSequence, newIsoformSequence));
						}
					}
					
					System.out.println(isoform.getCds().get(0).getStart() + " -- " + isoform.getCds().get(isoform.getCds().size()-1).getEnd());
					System.out.println(max_overlap);
					
					if(tx_orf_hash.containsKey(isoform_frame)){
						
						ArrayList<Interval2d> startEndList = tx_orf_hash.get(isoform_frame);
						
						Iterator<Interval2d> orfIter = startEndList.iterator();

						Integer[] scoresOrfs = new Integer[startEndList.size()];		
							
						for(int i=0; i<scoresOrfs.length; i++)
							scoresOrfs[i] = Integer.MAX_VALUE;
						
						int isoStart = cdsList.get(0).getStart(); 
						int isoEnd = cdsList.get(cdsList.size() - 1).getEnd(); 

						int orfIndex = 0;
						while(orfIter.hasNext()){
						
							Interval2d startEnd = orfIter.next();
							
							/***************************************************\ 
							| Check if start codon is shared with known isoform |
							\***************************************************/
								
							boolean equalStart = (gene.getStrand().equals("+") && startEnd.getStart() == isoStart) || (gene.getStrand().equals("-") && startEnd.getEnd() == isoEnd);
							boolean equalEnd = (gene.getStrand().equals("+") && startEnd.getEnd() == isoEnd) || (gene.getStrand().equals("-") && startEnd.getStart() == isoStart);
						
							if(equalStart){

								if(equalEnd){
									// best score 1 if start and end overlap with known isoform
									scoresOrfs[orfIndex] = 1;
									
								}else{
									// just the start overlaps with a known isoform
									scoresOrfs[orfIndex] = 2;
								}	
								
								if(scoresOrfs[orfIndex] < minScore){
									
									minScore = scoresOrfs[orfIndex];
								
								}
							}else{
								// orf is in the same reading frame as known isoform
								if(equalEnd){
//									 best score 1 if start and end overlap with known isoform
									scoresOrfs[orfIndex] = 3;
								}else{
									scoresOrfs[orfIndex] = 4;
								}
							}

							orfIndex++;
						}
						
						scoreMap.put(isoIndex, scoresOrfs);
					}
				}
				
				isoIndex++;
			}
			
			// Check the scoreList if there are some ORFs that are equal to one of the known isoforms
			if(scoreMap.size() == 0){
				/*
				Iterator<Integer> orfIter = tx_orf_hash.keySet().iterator();
				int length_longestORF = 0;
				Interval2d longest_orf_int = new Interval2d();
				
				longestCds = new ArrayList<Feature>();
				
				while(orfIter.hasNext()){
					
					int frame = orfIter.next();
					ArrayList<Interval2d> orfInts = tx_orf_hash.get(frame);
					
					for(int i=0; i<orfInts.size(); i++){
						
						Interval2d curOrf = orfInts.get(i);
						
						transcript.interval2dToCds(curOrf);
						
						int cds_length = transcript.get_cds_length();
							
						if((length_longestORF < cds_length) & (cds_length > min_orf_length)){
							
							length_longestORF = cds_length;
							longestCds = transcript.getCds();
							longest_orf_int = curOrf;
						}
					}
				}
				
				/** 
				 * Check if gene locus overlaps longest ORF by at least 50% 
				 
				if(longestCds != null){
					
				
					Interval2d geneInterval = new Interval2d(gene.getStart(), gene.getEnd());
					int overlap = geneInterval.overlap(longest_orf_int);
				
					// nestedness -> how much does the cds/transcript overlap with the corresponding gene locus
					double nestedness = (double) overlap / (double) longest_orf_int.length();
				
				
					if((nestedness > 0.5d)){
						transcript.setCds(longestCds);
						belongsToGene = true;
					}else{
						transcript.setCds(null);
						belongsToGene = false;
					}
				}else{
					transcript.setCds(null);
					belongsToGene = false;
				}*/
			}else{
				
				ArrayList<String> isoOrf_index = new ArrayList<String>();
				
				// Find out which combination had the minimal score
				Iterator<Integer> scoreMapIter = scoreMap.keySet().iterator();
				while(scoreMapIter.hasNext()){
					
					int isoformIndex = scoreMapIter.next();
					for(int j=0; j<scoreMap.get(isoformIndex).length; j++){
						
						if(scoreMap.get(isoformIndex)[j] == minScore){
							
							isoOrf_index.add(isoformIndex + "-" + j);
						}
					
					}
				}
				
				// Check if minScore is smaller than three and if there are more than one ORF-Isoform-combination take the longest
					
				Iterator<String> isoOrf_iter = isoOrf_index.iterator();
				int length_longestORF = 0;
				Interval2d longest_orf_int = new Interval2d();
				
				longestCds = new ArrayList<Feature>();
							
				while(isoOrf_iter.hasNext()){
						
					String[] iso_orf = isoOrf_iter.next().split("-");
					
					int iso = Integer.parseInt(iso_orf[0]);
					int orf = Integer.parseInt(iso_orf[1]);
					
					MetaFeature curIsoform = isoformList.get(iso);
					int isoform_frame = 0;
					
					if(gene.getStrand().equals("-")){
					
						isoform_frame = curIsoform.getCds().get(curIsoform.getCds().size() - 1).getEnd() % 3;
					
					}else{
					
						isoform_frame = curIsoform.getCds().get(0).getStart() % 3;
					
					}
					Interval2d curOrf = tx_orf_hash.get(isoform_frame).get(orf);
						
					// construct CDS coordinates for exon structure of transcript
						
					transcript.interval2dToCds(curOrf);
						
					int cds_length = transcript.get_cds_length();
					
					if(cds_length <= 0){
						
						System.err.println(transcript.getId() + " has ORF with " + cds_length);
						System.err.println(curOrf.toString());
						
						for(int i=0; i<transcript.getExons().size(); i++)
							System.err.println(transcript.getExons().get(i).toString());
						
						System.err.println(transcript.getCds().toString());
						
						throw new IllegalArgumentException("ORF length has to be positive. See above!");
					}
					
					if((length_longestORF < cds_length) & (min_orf_length < cds_length)){
						
						length_longestORF = cds_length;
						longestCds = transcript.getCds();
						longest_orf_int = curOrf;
					}
				}
					
				if(minScore <= 2){
					
					transcript.setCds(longestCds);
					
					if(transcript.get_cds_length() > min_orf_length){

						belongsToGene = true;
						
					}else{

						belongsToGene = false;
						transcript.setCds(null);
					}
				}else{
					// If minScore is three check if transcript in same overlapping at least 50% of the knwon gene
					Interval2d geneInterval = new Interval2d(gene.getStart(), gene.getEnd());
					int overlap = geneInterval.overlap(longest_orf_int);
					
					// nestedness -> how much does the cds/transcript overlap with the corresponding gene locus
					double nestedness = (double) overlap / (double) longest_orf_int.length();
					transcript.setCds(longestCds);
					
					if((nestedness > 0.5d) & (transcript.get_cds_length() > min_orf_length)){
						belongsToGene = true;
					}else{
						transcript.setCds(null);
						belongsToGene = false;
					}
				}
			}
		}

		return belongsToGene;
	}
		
	public static void writeFASTA(HashMap<String,String> seqs, String filename, int offset) throws IOException{
		
		BufferedWriter output = new BufferedWriter(new FileWriter(filename));
		Iterator<String> iter = seqs.keySet().iterator();
		
		while(iter.hasNext()){
			
			String id = iter.next();
			String seq = seqs.get(id);
			
			int start = 0;
			output.write(">"+id+" ("+seq.length()+")\n");
			for(int i=0; i<seq.length()/offset; i++){
				output.write(seq.substring(start, start+offset)+"\n");
				start += offset;
			}
			if(start<seq.length()){
				output.write(seq.substring(start, seq.length())+"\n");
			}
			
		}
		output.close();
	}
	
	public static void writeFASTASeqHash(HashMap<String,Sequence> seqs, String filename, int offset) throws IOException{
		
		BufferedWriter output = new BufferedWriter(new FileWriter(filename));
		Iterator<String> iter = seqs.keySet().iterator();
		
		while(iter.hasNext()){
			
			String id = iter.next();
			String seq = seqs.get(id).getSeq();
			
			int start = 0;
			output.write(">"+id+" ("+seq.length()+")\n");
			for(int i=0; i<seq.length()/offset; i++){
				output.write(seq.substring(start, start+offset)+"\n");
				start += offset;
			}
			if(start<seq.length()){
				output.write(seq.substring(start, seq.length())+"\n");
			}
			
		}
		output.close();
	}

	public static void writeFASTA(Sequence seq, String filename, int offset) throws IOException{
		
		BufferedWriter output = new BufferedWriter(new FileWriter(filename));
			
		String seqData = seq.getSeq();
		
		int start = 0;
		output.write(">"+seq.getId()+" ("+seqData.length()+")\n");
		
		for(int i=0; i<seqData.length()/offset; i++){
		
			output.write(seqData.substring(start, start+offset)+"\n");
			start += offset;
		}
		
		if(start < seqData.length()){
		
			output.write(seqData.substring(start, seqData.length())+"\n");
		}
			
		output.close();
	}

	public static String getFASTAentry(Sequence seq, int offset){
		
		StringBuilder output = new StringBuilder();

		int start = 0;
		output.append(">"+seq.getId()+" ("+seq.getSeq().length()+")\n");
		for(int i=0; i<seq.getSeq().length()/offset; i++){
			output.append(seq.getSeq().substring(start, start+offset).toUpperCase()+"\n");
			start += offset;
		}
		if(start<seq.getSeq().length()){
			output.append(seq.getSeq().substring(start, seq.getSeq().length()).toUpperCase()+"\n");
		}
		
		return output.toString();
	}
	
	public static String getFASTAentry(Sequence seq, int offset, String description){
		
		StringBuilder output = new StringBuilder();

		int start = 0;
		output.append(">"+seq.getId()+" ("+seq.getSeq().length()+") | "+ description +"\n");
		for(int i=0; i<seq.getSeq().length()/offset; i++){
			output.append(seq.getSeq().substring(start, start+offset).toUpperCase()+"\n");
			start += offset;
		}
		if(start<seq.getSeq().length()){
			output.append(seq.getSeq().substring(start, seq.getSeq().length()).toUpperCase()+"\n");
		}
		
		return output.toString();
	}
	
	public void writeSeq(BufferedWriter output, String seq, String id, int offset) throws IOException{
		
		int start = 0;
		
		output.write(">"+id+"("+seq.length()+")\n");
		for(int i=0; i<seq.length()/offset; i++){
			output.write(seq.substring(start, start+offset)+"\n");
			start += offset;
		}
		if(start<seq.length()){
			output.write(seq.substring(start, seq.length())+"\n");
		}
	}
	
	public boolean checkExonBorders(Sequence seq, ArrayList<Feature> exons){
	
		boolean correct = true;
		
		if(exons != null){
			
			if(exons.size() > 1){
				
				Collections.sort(exons);
				
				
				
				for(int i=0; i<exons.size()-1; i++){
					
					
//					if(exons.get(i).getId().equals("atcg00180.1")){
//						System.out.println("Seq: " + seq.toString());
//					}
					// + ": " + exons.get(i).getEnd() + " - " + exons.get(i).getEnd()+2);
					
					String start_gt = seq.getSeq().substring(exons.get(i).getEnd(), exons.get(i).getEnd()+2);
					String end_ag = seq.getSeq().substring(exons.get(i+1).getStart()-3, exons.get(i+1).getStart()-1);
					
					// 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
					correct &= 	(start_gt.equals("GT") & end_ag.equals("AG")) | 
								(start_gt.equals("CT") & end_ag.equals("AC")) |
								(start_gt.equals("GC") & end_ag.equals("AG")) |
								(start_gt.equals("CT") & end_ag.equals("GC")) |
								(start_gt.equals("AT") & end_ag.equals("AC")) |
								(start_gt.equals("GT") & end_ag.equals("AT"));
//					System.out.println(start_gt + " - " + end_ag);
					
					/*if(exons.get(i).getStrand().equals("-")){
						
						correct &= start_gt.equals("CT") & end_ag.equals("AC");
//						System.out.println(start_gt + " - " + end_ag);
					}else if(exons.get(i).getStrand().equals("+")){

						correct &= start_gt.equals("GT") & end_ag.equals("AG");
						
					}else{
						
						correct = false;
						break;
						
					}*/
				}
					
			}
			
		}
		
		return correct;
		
	}
	
	
	
}


