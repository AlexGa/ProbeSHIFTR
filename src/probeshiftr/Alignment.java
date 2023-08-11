package probeshiftr;

import java.util.HashMap;

public class Alignment {

	private int[][] d;
	String s1;
	String s2;
	String[] alignment;

	static int k = -1;
	HashMap<Integer, String[]> allAls;

	public Alignment(String s1, String s2) {

		this.s1 = s1;
		this.s2 = s2;
		alignment = new String[2];
	}

	public String[] getAlignment() {
		return alignment;
	}
	
	private void fillDpMatrix() {

		allAls = new HashMap<>();
		d = new int[s1.length() + 1][s2.length() + 1];

		/* Initialize DP-Matrix*/
		for (int i = 0; i < d.length; i++) {
			d[i][0] = i;
		}

		for (int j = 0; j < d[0].length; j++) {
			d[0][j] = j;
		}

		for (int i = 1; i < d.length; i++) {
			for (int j = 1; j < d[0].length; j++) {
				d[i][j] = Math.min(d[i][j - 1] + 1, Math.min(d[i - 1][j] + 1, d[i - 1][j - 1] + dist(s1.charAt(i - 1), s2.charAt(j - 1))));
			}
		}

	}

	public double identity() {
		
		double identity = 0.0d;
		
		if(alignment != null) {
			for(int i = 0; i < alignment[0].length(); i++) 
				identity += alignment[0].charAt(i) == alignment[1].charAt(i) ? 1.0d : 0.0d;
			
			//identity = identity/(double) alignment[0].length();
		}
		
		return identity;
	}
	public void globalAlignment(){

		fillDpMatrix();

		String sa1 = "", sa2 = "";
		
		int i = s1.length(), j = s2.length();
		
		while(i>0 & j>0){
			
			if(d[i][j] == d[i-1][j-1] + dist(s1.charAt(i-1), s2.charAt(j-1))){
				
				sa1 = s1.charAt(i - 1) + sa1;
				sa2 = s2.charAt(j - 1) + sa2;
				i--;
				j--;
			}else if(d[i][j] == d[i][j-1] + 1){
				sa1 = "-" + sa1;
				sa2 = s2.charAt(j - 1) + sa2;
				j--;
			}else{
				sa1 = s1.charAt(i - 1) + sa1;
				sa2 = "-" + sa2;
				i--;
			}
		}
		
		while(i>0){
			sa1 = s1.charAt(i-1) + sa1;
			sa2 = "-" + sa2;
			i--;
		}
		
		while(j>0){
			sa1 = "-" + sa1;
			sa2 = s2.charAt(j-1) + sa2;
			j--;
		}

		alignment[0] = sa1;
		alignment[1] = sa2;

	}
	
	
	private int dist(char i, char j){
		
		return i == j ? 0 : 1;
	}

}
