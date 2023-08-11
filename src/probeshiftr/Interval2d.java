package probeshiftr;

public class Interval2d implements Comparable<Interval2d>{
	
	private int start;
	private int end;
	
	public Interval2d(){
		
	}
	
	public Interval2d(int start, int end) {
		super();
		this.start = start;
		this.end = end;
	}
	
	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	@Override
	public int compareTo(Interval2d i) {
		
		int comp = Integer.compare(this.start, i.start);
		
		if(comp == 0){
			comp = -1 * Integer.compare(this.end, i.end);
		}
		
		return comp;
	}

	@Override
	public String toString() {
		
		return "Interval2d [start=" + start + ", end=" + end + "]";
	}
	
	@Override
	public boolean equals(Object obj) {
		
		if (this == obj)
            return true;
        
		if (obj == null)
            return false;
        
		if (getClass() != obj.getClass())
            return false;
        
		Interval2d toTest = (Interval2d) obj;
        
		if (this.start == toTest.getStart() & this.end == toTest.getEnd())
            return true;
        
		return false;
		
	}
	
	public int compareTo(int position) {
		
		
//		System.out.println("\t -> COMPARISON: " + this.toString() + " for " + position);
		
		int comp =  position - this.getStart();
		
//		System.out.println("\t -> COMPARISON: \t Before IF: " + comp);
		
		if(comp >= 0) {
			if(this.end >= position) {		
				comp = 0;
			}else {
				comp = position - this.getEnd();
			}
		}
		
//		System.out.println("\t -> COMPARISON: \t After IF: " + comp);
		
		return comp;
	}

	
	public int overlap(Interval2d i){
		
		int overlap_length = 0;
		
		if((i.start > this.end) || (this.start > i.end)){
			overlap_length = 0;
		}else{
			int ov_start = Math.max(this.start, i.start);
			int ov_end = Math.min(this.end, i.end);
			overlap_length = ov_end - ov_start;
		}
		
		return overlap_length;
	}
	
	public int length(){
		
		return this.end - this.start + 1;
	}
	
	public static int binarySearchOverlap(Interval2d[] intArray, int position){
		
		int l = 0;
		int r = intArray.length -1;
		
		int nothing = -1; 

		while(l <= r){
			
			int mid = l + ((r - l) / 2);

//			System.out.println("BINARY SEARCH: " + "\tl = " + l + "\tr = " + r + "\tmid = " + mid);
			int comp = intArray[mid].compareTo(position);

//			System.out.println("BINARY SEARCH: " + "\tcomp = " + comp);

			if(comp == 0){
//				System.out.println("BINARY SEARCH: " + mid + "\t" + intArray[mid].toString());
				return mid;
			}else{
				if(comp > 0){
					l = mid + 1;
				}else{
					r = mid - 1;
				}
			}
		}
//		System.out.println("BINARY SEARCH: Final interval -> :" + "\tl = " + l + "\tr = " + r);
//		System.out.println("BINARY SEARCH: NO RESULT!");
		
		return nothing;
	}
}
