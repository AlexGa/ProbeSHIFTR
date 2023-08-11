package probeshiftr;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.lang3.StringUtils;

/** 
 * 
 * Code from guttmanlab
 * github.com/GuttmanLab/guttmanlab-capture/blob/master/src/capture/filter/LowComplexityFilter.java
 * 
 * 
 **/

public class ComplexityFilter implements OligoFilter{

	public String name = "LowComplexity";

	final char chars[] = new char[] {'A','C','T','G'};
	// Setting K = 6 will check for sequencing containing repeats of up to k = 6 bases (e.g., (GGCACA GGCACA ...))
	private int maxK = 6;
	
	// Sets the length of the repeat that must be present to call a match
	private int[] minLength = new int[] { 10, 10, 12, 12, 18, 18 };
	
	private Map<Integer, HashSet<String>> hashes = null;
	
	public Map<Integer, HashSet<String>> LowComplexityFilter() {
		
		return getLowComplexityHashes(maxK, minLength);
	}
	
	public ComplexityFilter() {
		hashes = getLowComplexityHashes(maxK, minLength);
	}
	
	/**
	 * @param k
	 * @param minLengths
	 * @return
	 * Set up hash sets containing all strings of minLengths (n) bases from k-mer repeats
	 */
	private Map<Integer, HashSet<String>> getLowComplexityHashes(int k, int[] minLengths) {
		
		Map<Integer, HashSet<String>> hashes = new TreeMap<Integer, HashSet<String>>();
		
		for (int repeatK = 0; repeatK < k; repeatK++) {
		
			HashSet<String> currHash = hashes.containsKey(minLengths[repeatK]) ?
				currHash = hashes.get(minLengths[repeatK]) : new HashSet<String>();
	
			List<String> allRepeats = new ArrayList<String>();
			addAllKmerRepeats(allRepeats, "", repeatK+1);
			
			for (String repeatElement : allRepeats) {
			
				addRepeatToHash(currHash, repeatElement, minLengths[repeatK]);
			
			}
			
			hashes.put(minLengths[repeatK], currHash);
		}
		return hashes;
	}
	
	/**
	 * @param allRepeats
	 * @param base
	 * @param k
	 * Recursive function to generate all possible k-mer sequences
	 */
	private void addAllKmerRepeats(List<String> allRepeats, String base, int k) {
		if (base.length() >= k) {
			//logger.info("Adding " + base + " to list.");
			allRepeats.add(base);
		} else {
			for (int i = 0; i < chars.length; i++)
				addAllKmerRepeats(allRepeats, base + chars[i], k);
		}
	}
	
	
	/**
	 * @param hash
	 * @param repeatElement
	 * @param hashLength
	 * Take a repeat element (e.g., GGA) and add all strings of hashLength bases to the hash
	 */
	private void addRepeatToHash(HashSet<String> hash, String repeatElement, int hashLength) {
		repeatElement = StringUtils.repeat(repeatElement, hashLength*2);
		//logger.info("Adding " + repeatElement);
		for (int i = 0; i < hashLength; i++) {
			hash.add(repeatElement.substring(i,hashLength));
		}
	}

	private boolean sequenceContainsNmerInHash(String seq, int n) {
		for (int charIndex = 0; charIndex < seq.length() - n + 1; charIndex++) {
			String substr = seq.substring(charIndex, charIndex + n);
			if (hashes.get(n).contains(substr)) return true;
		}
		return false;
	}
	
	@Override
	public boolean filterOligo(String oligo) {
		for (int period = 0; period < maxK; period++) {
			if (sequenceContainsNmerInHash(oligo, minLength[period])) return true;
		}
		return false;
	}
}
