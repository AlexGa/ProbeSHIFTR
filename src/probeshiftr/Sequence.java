package probeshiftr;

import java.util.HashMap;


public class Sequence {

	private String id;
	private String seq;
	
	public Sequence() {
		super();
	}
	
	public Sequence(String id, String seq) {
		super();
		this.id = id;
		this.seq = seq;
	}

	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getSeq() {
		return seq;
	}
	public void setSeq(String seq) {
		this.seq = seq;
	}
	
	public Sequence reverse() throws NullPointerException{

		if(this.seq != null){
			
			try{
			
				return new Sequence(this.id, new StringBuffer(this.seq).reverse().toString());
		
			}catch( Exception e ){
			
				RuntimeException runException = new RuntimeException( e.getMessage() );
				runException.setStackTrace( e.getStackTrace() );
				throw runException;
				
			}
			
		}else{
			
			throw new NullPointerException("Attriubute seq from Sequence object is null!");
		
		}
	}
	
	public Sequence complement() throws NullPointerException{
		

		if(this.seq != null){
			
			try{
				
				String compSeq = this.seq.toUpperCase(); 
				
				compSeq = compSeq.replaceAll("A", "t");
				compSeq = compSeq.replaceAll("C", "g");
				compSeq = compSeq.replaceAll("G", "c");
				compSeq = compSeq.replaceAll("T", "a");
				
				compSeq = compSeq.replaceAll("Y", "r");
				compSeq = compSeq.replaceAll("R", "y");
				compSeq = compSeq.replaceAll("W", "w");
				compSeq = compSeq.replaceAll("S", "s");
				compSeq = compSeq.replaceAll("K", "m");
				compSeq = compSeq.replaceAll("M", "k");
				
				compSeq = compSeq.replaceAll("D", "h");
				compSeq = compSeq.replaceAll("H", "d");
				compSeq = compSeq.replaceAll("B", "v");
				compSeq = compSeq.replaceAll("V", "b");
				compSeq = compSeq.replaceAll("X", "x");
				compSeq = compSeq.replaceAll("N", "n");
				
				
				compSeq = compSeq.toUpperCase();
				
				return new Sequence(this.id, compSeq);
				
			}catch( Exception e ){
			
				RuntimeException runException = new RuntimeException( e.getMessage() );
				runException.setStackTrace( e.getStackTrace() );
				throw runException;
				
			}
			
		}else{
			
			throw new NullPointerException("Attriubute seq from Sequence object is null!");
		
		}
	}
	
	public Sequence getReverseComplement(){
		
		return this.reverse().complement();
		
	}

	@Override
	public String toString() {
		return "Sequence [id=" + id + ", seq=" + seq + "]";
	}
	
	public Sequence translate(){
		
		if(this.seq.length() % 3 != 0){
//			System.err.println(this.getSeq());
//			throw new IllegalArgumentException("Sequence \"" + this.id + "\" with length "+this.seq.length() +" cannot be divided by 3!");
			System.err.println("Sequence \"" + this.id + "\" with length "+this.seq.length() +" cannot be divided by 3!");
		}
		
		StringBuilder proteinSeq = new StringBuilder();
		
		HashMap<String, String> codon2aa = new HashMap<String, String>();
	    
		codon2aa.put("GGA","G");
	    codon2aa.put("GGC","G");
	    codon2aa.put("GGG","G");
	    codon2aa.put("GGT","G");
	    
	    codon2aa.put("GGR","G");
	    codon2aa.put("GGY","G");
	    codon2aa.put("GGM","G");
	    codon2aa.put("GGK","G");
	    codon2aa.put("GGS","G");
	    codon2aa.put("GGW","G");
	    codon2aa.put("GGH","G");
	    codon2aa.put("GGB","G");
	    codon2aa.put("GGV","G");
	    codon2aa.put("GGD","G");
	    codon2aa.put("GGN","G");
	    
	    codon2aa.put("GCA","A");
	    codon2aa.put("GCC","A");
	    codon2aa.put("GCG","A");
	    codon2aa.put("GCT","A");
	    
	    codon2aa.put("GCR","A");
	    codon2aa.put("GCY","A");
	    codon2aa.put("GCM","A");
	    codon2aa.put("GCK","A");
	    codon2aa.put("GCS","A");
	    codon2aa.put("GCW","A");
	    codon2aa.put("GCH","A");
	    codon2aa.put("GCB","A");
	    codon2aa.put("GCV","A");
	    codon2aa.put("GCD","A");
	    codon2aa.put("GCN","A");

	    codon2aa.put("GTA","V");
	    codon2aa.put("GTC","V");
	    codon2aa.put("GTG","V");
	    codon2aa.put("GTT","V");
	    
	    codon2aa.put("GTR","V");
	    codon2aa.put("GTY","V");
	    codon2aa.put("GTM","V");
	    codon2aa.put("GTK","V");
	    codon2aa.put("GTS","V");
	    codon2aa.put("GTW","V");
	    codon2aa.put("GTH","V");
	    codon2aa.put("GTB","V");
	    codon2aa.put("GTV","V");
	    codon2aa.put("GTD","V");
	    codon2aa.put("GTN","V");

	    codon2aa.put("CTA","L");
	    codon2aa.put("CTC","L");
	    codon2aa.put("CTG","L");
	    codon2aa.put("CTT","L");
	    
	    codon2aa.put("CTR","L");
	    codon2aa.put("CTY","L");
	    codon2aa.put("CTM","L");
	    codon2aa.put("CTK","L");
	    codon2aa.put("CTS","L");
	    codon2aa.put("CTW","L");
	    codon2aa.put("CTH","L");
	    codon2aa.put("CTB","L");
	    codon2aa.put("CTV","L");
	    codon2aa.put("CTD","L");
	    codon2aa.put("CTN","L");
	    
	    codon2aa.put("TTA","L");
	    codon2aa.put("TTG","L");
	    codon2aa.put("TTR","L");
	    
	    codon2aa.put("ATA","I");
	    codon2aa.put("ATC","I");
	    codon2aa.put("ATT","I");
	    codon2aa.put("ATY","I");
	    codon2aa.put("ATM","I");
	    codon2aa.put("ATW","I");
	    codon2aa.put("ATH","I");

	    codon2aa.put("CCA","P");
	    codon2aa.put("CCC","P");
	    codon2aa.put("CCG","P");
	    codon2aa.put("CCT","P");
	    
	    codon2aa.put("CCR","P");
	    codon2aa.put("CCY","P");
	    codon2aa.put("CCM","P");
	    codon2aa.put("CCK","P");
	    codon2aa.put("CCS","P");
	    codon2aa.put("CCW","P");
	    codon2aa.put("CCH","P");
	    codon2aa.put("CCB","P");
	    codon2aa.put("CCV","P");
	    codon2aa.put("CCD","P");
	    codon2aa.put("CCN","P");

	    codon2aa.put("TTC","F");
	    codon2aa.put("TTT","F");
	    codon2aa.put("TTY","F");

	    codon2aa.put("TAC","Y");
	    codon2aa.put("TAT","Y");
	    codon2aa.put("TAY","Y");
	    
	    codon2aa.put("TGC","C");
	    codon2aa.put("TGT","C");
	    codon2aa.put("TGY","C");

	    codon2aa.put("ATG","M");

	    codon2aa.put("CAC","H");
	    codon2aa.put("CAT","H");
	    codon2aa.put("CAY","H");

	    codon2aa.put("AAA","K");
	    codon2aa.put("AAG","K");
	    codon2aa.put("AAR","K");
	    
	    codon2aa.put("AGA","R");
	    codon2aa.put("AGG","R");
	    codon2aa.put("AGR","R");
	    
	    codon2aa.put("CGA","R");
	    codon2aa.put("CGC","R");
	    codon2aa.put("CGG","R");
	    codon2aa.put("CGT","R");
	    
	    codon2aa.put("CGR","R");
	    codon2aa.put("CGY","R");
	    codon2aa.put("CGM","R");
	    codon2aa.put("CGK","R");
	    codon2aa.put("CGS","R");
	    codon2aa.put("CGW","R");
	    codon2aa.put("CGH","R");
	    codon2aa.put("CGB","R");
	    codon2aa.put("CGV","R");
	    codon2aa.put("CGD","R");
	    codon2aa.put("CGN","R");

	    codon2aa.put("TGG","W");

	    codon2aa.put("AGC","S");
	    codon2aa.put("AGT","S");
	    codon2aa.put("AGY","S");
	    
	    codon2aa.put("TCA","S");
	    codon2aa.put("TCC","S");
	    codon2aa.put("TCG","S");
	    codon2aa.put("TCT","S");
	    
	    codon2aa.put("TCR","S");
	    codon2aa.put("TCY","S");
	    codon2aa.put("TCM","S");
	    codon2aa.put("TCK","S");
	    codon2aa.put("TCS","S");
	    codon2aa.put("TCW","S");
	    codon2aa.put("TCH","S");
	    codon2aa.put("TCB","S");
	    codon2aa.put("TCV","S");
	    codon2aa.put("TCD","S");
	    codon2aa.put("TCN","S");

	    codon2aa.put("ACA","T");
	    codon2aa.put("ACC","T");
	    codon2aa.put("ACG","T");
	    codon2aa.put("ACT","T");
	    
	    codon2aa.put("ACR","T");
	    codon2aa.put("ACY","T");
	    codon2aa.put("ACM","T");
	    codon2aa.put("ACK","T");
	    codon2aa.put("ACS","T");
	    codon2aa.put("ACW","T");
	    codon2aa.put("ACH","T");
	    codon2aa.put("ACB","T");
	    codon2aa.put("ACV","T");
	    codon2aa.put("ACD","T");
	    codon2aa.put("ACN","T");
	    
	    codon2aa.put("GAC","D");
	    codon2aa.put("GAT","D");
	    codon2aa.put("GAY","D");
	    
	    codon2aa.put("GAA","E");
	    codon2aa.put("GAG","E");
	    codon2aa.put("GAR","E");
	    
	    codon2aa.put("AAC","N");
	    codon2aa.put("AAT","N");
	    codon2aa.put("AAY","N");
	    
	    codon2aa.put("CAA","Q");
	    codon2aa.put("CAG","Q");
	    codon2aa.put("CAR","Q");
	    
	    codon2aa.put("RAY","B");

	    codon2aa.put("SAR","Z");
	    
	    codon2aa.put("TGA","");
	    codon2aa.put("TAA","");
	    codon2aa.put("TAG","");
	    
	    for(int i=0; i<this.seq.length()-2; i+=3){
	    	String codon = this.seq.substring(i, i+3).toUpperCase();
	    	if(codon2aa.containsKey(codon)){
	    		proteinSeq.append(codon2aa.get(codon));
	    	}else{
	    		proteinSeq.append("X");
	    		System.err.println("Codon "+ codon + " could not be found in translation table! Set amino acid X for unknown codon at position " + (i+1) + " in sequence " + this.id);
//	    		throw new IllegalArgumentException("Codon "+ codon + " could not be found in translation table! Set X as amino acid");
	    	}
	    	
	    }
	    
	    return new Sequence(this.id, proteinSeq.toString());

	}
}
