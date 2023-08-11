package probeshiftr;

public class PolybaseFilter {

	/*
	 * Polybase Filter 
	 * 
	 * */
	
	private String basesToFilter; // ACGT
	private int cutoff;  // rel_cutoff * repeat_length e.g. 0.008 * 15
	private int repeatLength;
	
	public PolybaseFilter() {}

	public PolybaseFilter(String basesToFilter, double rel_cutoff, int repeatLength) {
		this.basesToFilter = basesToFilter;
		this.cutoff = (int) (rel_cutoff * (double) repeatLength);
		
		this.repeatLength = repeatLength;
	}
	
	public PolybaseFilter(String basesToFilter, int cutoff, int repeatLength) {
		this.basesToFilter = basesToFilter;
		this.cutoff = cutoff;
		
		this.repeatLength = repeatLength;
	}
	
	public boolean filterOligo(String oligo) {
		
		char[] seq = oligo.toUpperCase().toCharArray();
		
		for (char c : basesToFilter.toCharArray()) {
			
			int charMatchesInWindow = 0;
			for (int i = 0; i < seq.length; i++) {
				if (seq[i] == c) {
					charMatchesInWindow++;
				}
				
				if (i >= repeatLength) {
					if (seq[i-repeatLength] == c) {
						charMatchesInWindow--;
					}
				}
				
				if (i >= repeatLength - 1) {
					if (charMatchesInWindow >= cutoff) {
						return true;
					}
				}
			}
		}
		return false;
	}
	
}
