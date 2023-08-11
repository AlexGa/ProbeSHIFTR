package probeshiftr;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class PerformBLAT {
	
	
	
	private String[] parameterList;
	
	
	public PerformBLAT(String blatPath, int minMatch, int minScore, String dbSeq, String queryFa, String outputFile) {
		
		this.parameterList = new String[6];
		
		this.parameterList[0] =  blatPath;
		this.parameterList[1] =  "-minMatch=" + String.valueOf(minMatch);
		this.parameterList[2] =  "-minScore=" + String.valueOf(minScore);
		this.parameterList[3] =  dbSeq;
		this.parameterList[4] =  queryFa;
		this.parameterList[5] =  outputFile;
		
	}
	
	public int runBLAT() throws IOException, InterruptedException {
		
		ProcessBuilder pb = new ProcessBuilder(this.parameterList);
		
		Process blatProc = pb.start();
		
		BufferedReader blatStdOut = new BufferedReader(new InputStreamReader(blatProc.getInputStream()));

		String line = null;
		
		
		while((line = blatStdOut.readLine()) != null) {
			
			System.out.println(line);
		}

		int exitCode = blatProc.waitFor();
		
		return(exitCode);
		
	}	
}
