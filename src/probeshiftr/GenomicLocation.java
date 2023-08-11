package probeshiftr;

public class GenomicLocation {

	private boolean genic;
	private boolean exon;
	private int index; // negativ if genic is false
	
	public GenomicLocation(boolean genic, boolean exon, int index) {
		super();
		this.genic = genic;
		this.exon = exon;
		this.index = index;
	}

	public boolean isGenic() {
		return genic;
	}

	public boolean isExon() {
		return exon;
	}

	public int getIndex() {
		return index;
	}

	@Override
	public String toString() {
		return "GenomicLocation [genic=" + genic + ", exon=" + exon + ", index=" + index + "]";
	}
	
	
}
