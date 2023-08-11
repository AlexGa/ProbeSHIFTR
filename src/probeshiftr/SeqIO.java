package probeshiftr;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

public class SeqIO {

	public SeqIO(){
		
	}
	
	public static HashMap<String, Sequence> readFastaGenom(String filename) throws IOException {

		HashMap<String, Sequence> seqMap = new HashMap<String,Sequence>();
		BufferedReader input = null;
		InputStream gzstream = null;
		
		try{
			
			if(filename.contains(".gz")){
				
				InputStream filestream = new FileInputStream(filename);
				gzstream = new GZIPInputStream(filestream);
				Reader decode = new InputStreamReader(gzstream);
				input = new BufferedReader(decode);
				
			}else{
				
				input = new BufferedReader(new FileReader(filename));
			}

			String line = null;
			String id = "";
			Sequence seq = null;

			StringBuffer strFile = new StringBuffer();

			while ((line = input.readLine()) != null) {

				if (line.length() <= 0) {
					
					continue;
					
				}
				if (line.charAt(0) == '>') {
					
					if(id != ""){
					
						seq = new Sequence(id, strFile.toString());
						seqMap.put(id, seq);
						id = "";
						
					}
					
					id = line.split(" ")[0].replace(">", "");//.toLowerCase();
					strFile = new StringBuffer();
			
				}else{
					strFile.append(line.replaceAll("\n", ""));//.toUpperCase());
				}
				
			}

			seq = new Sequence(id, strFile.toString());
			seqMap.put(id, seq);
			
			input.close();
			return seqMap;
			
		} catch ( Exception e){
			
			RuntimeException runException = new RuntimeException( e.getMessage() );
			runException.setStackTrace( e.getStackTrace() );
			throw runException;
			
		}finally{

			if(gzstream != null){
				
				gzstream.close();
			}
			
			if(input != null){

				input.close();
				
			}
		}
	}
}
