package probeshiftr;

public class RepeatFilter implements OligoFilter{

	/** 
	 * 
	 * Masked repeats and N in sequence 
	 * 
	 * **/
	
	private boolean includesLower;
	private boolean includesN;
	private double maxPct;
	
	public RepeatFilter() {}
	
	/**
	 * @param maxRepeatPct Maximum allowable percentage of repeat masked bases (0 to 1)
	 * @param includeLowerCase Lower case bases count as repeats
	 * @param includeN Ns count as repeats
	 */
	public RepeatFilter(double maxRepeatPct, boolean includeLowerCase, boolean includeN) {
		
		if(maxRepeatPct < 0 || maxRepeatPct > 1) {
		
			throw new IllegalArgumentException("Max repeat percentage must be between 0 and 1");
		
		}
		
		if(!includeLowerCase && !includeN) {
		
			throw new IllegalArgumentException("Must include at least one: lowercase or N");
		
		}
		
		maxPct = maxRepeatPct;
		includesLower = includeLowerCase;
		includesN = includeN;
	}
	
	@Override
	public boolean filterOligo(String oligo) {
		
		int size = oligo.length();
		int numRepeats = 0;
		
		for(int i=0; i<size; i++) {
			
			char c = oligo.charAt(i);
			
			if(includesLower && !Character.isUpperCase(c)) {
				numRepeats++;
				continue;
			}
			
			if(includesN && (c == 'n' || c == 'N')) {
				numRepeats++;
				continue;
			}
		}
		double pct = (double) numRepeats / (double) size;
		
		return pct > maxPct;
	}

	
}
