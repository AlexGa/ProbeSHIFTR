package probeshiftr;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

public class Rinterface {
	
	/**
	anno_file <- args[1] # gtf annotation compressed or uncompressed
	blat_dir <- args[2] # directory containing the BLAT output created from java framework
	outDir <- args[3] # output directory containing the final oligo sets
	plot_dir <- args[4] # directory containing the analysis plots for different thresholds

	oligo_fasta_file <- args[5] # fasta file containing oligo sequences from java program
	oligo_length <- args[6] # length of oligos -> infer from sequence length!
	target_fasta_file <- args[7] # fasta file containing the RNA target sequences
	*/
	
	private String[] parameterList;
	
	
	public Rinterface(String R_bin_path, String Rscript_path, String blat_dir, 
					   String outDir, String plot_dir, String oilgo_fasta_file,
					   int oligo_length, String target_fasta_file, String... anno_file) {
		
		
		if(anno_file[0] != null) {
			this.parameterList = new String[9];
			this.parameterList[8] =  anno_file[0];
		}else {
			this.parameterList = new String[8];
		}
		
		this.parameterList[0] =  R_bin_path;
		this.parameterList[1] =  Rscript_path;
		this.parameterList[2] =  blat_dir;
		this.parameterList[3] =  outDir;
		this.parameterList[4] =  plot_dir;
		this.parameterList[5] =  oilgo_fasta_file;
		this.parameterList[6] =  String.valueOf(oligo_length);
		this.parameterList[7] =  target_fasta_file;
		
	}
	
	public int runR() throws IOException, InterruptedException {

		InputStream inputStream = Rinterface.class.getResourceAsStream(this.parameterList[1]);
		
        if (inputStream == null) {
        	return(1);
        	
        }
        
        // Read the R script content from the resource
		ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
		byte[] buffer = new byte[4096];

		int bytesRead;

		while ((bytesRead = inputStream.read(buffer)) != -1) {
			outputStream.write(buffer, 0, bytesRead);
		}

		String scriptContent = new String(outputStream.toByteArray(), StandardCharsets.UTF_8);

		// Write the R script content to a temporary file
		File tempScriptFile = File.createTempFile("designOligos", ".R");

		try (PrintWriter writer = new PrintWriter(new FileWriter(tempScriptFile))) {
			writer.write(scriptContent);
		}

		this.parameterList[1] = tempScriptFile.getAbsolutePath();
		
		System.out.println("Path to script: " + this.parameterList[1]);
		
		ProcessBuilder pb = new ProcessBuilder(this.parameterList);

		System.out.println(Arrays.toString(this.parameterList));

		Process Rproc = pb.start();

		BufferedReader rStdOut = new BufferedReader(new InputStreamReader(Rproc.getInputStream()));

		String line = null;

		while ((line = rStdOut.readLine()) != null) {

			System.out.println(line);
		}

		int exitCode = Rproc.waitFor();

		// Clean up the temporary file
//		tempScriptFile.delete();

		return (exitCode);
	}	
}
