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
					   String outDir, String plot_dir, String oligo_fasta_file,
					   int oligo_length, String target_fasta_file, boolean check_within, 
					   int fragment_size, String... anno_files) {
		
		
		if(anno_files[0] != null) {
			
			this.parameterList = new String[11];
			this.parameterList[10] = anno_files[0];
			
		}else {
			
			this.parameterList = new String[10];
		
		}
		
		this.parameterList[0] =  R_bin_path;
		this.parameterList[1] =  Rscript_path;
		this.parameterList[2] =  blat_dir;
		this.parameterList[3] =  outDir;
		this.parameterList[4] =  plot_dir;
		this.parameterList[5] =  oligo_fasta_file;
		this.parameterList[6] =  String.valueOf(oligo_length);
		this.parameterList[7] =  target_fasta_file;
		this.parameterList[8] =  String.valueOf(check_within);
		this.parameterList[9] =  String.valueOf(fragment_size);
		
	}

	public int runR(boolean verbose) throws IOException, InterruptedException {
		
		InputStream inputStream = Rinterface.class.getResourceAsStream(this.parameterList[1]);
		
        if (inputStream == null) {
        	System.err.println("An error has occured during the execution of R...");
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

		int exitCode = 0;
		 try {
	            // Start the R process
			 	
	            ProcessBuilder processBuilder = new ProcessBuilder(this.parameterList);
	            Process Rproc = processBuilder.start();

	            if(verbose) {
	            	
	            	// Create readers for the standard output and error streams
		            BufferedReader rStdOut = new BufferedReader(new InputStreamReader(Rproc.getInputStream()));
		            BufferedReader rStdErr = new BufferedReader(new InputStreamReader(Rproc.getErrorStream()));

		         // Create threads to handle the output and error streams
		            Thread outputThread = new Thread(new Runnable() {
		                @Override
		                public void run() {
		                    readStream(rStdOut, "OUTPUT (Rscript)");
		                }
		            });
		            Thread errorThread = new Thread(new Runnable() {
		                @Override
		                public void run() {
		                    readStream(rStdErr, "ERROR (Rscript)");
		                }
		            });

		            // Start the threads
		            outputThread.start();
		            errorThread.start();
		            
		            // Wait for the R process and threads to finish
		            exitCode = Rproc.waitFor();
		            outputThread.join();
		            errorThread.join();
		            System.out.println("R process exited with code: " + exitCode);
		            
	            }else {
	            	

		            // Wait for the R process to complete
		            exitCode = Rproc.waitFor();
	            	
	            }
	            
	        

	        } catch (IOException | InterruptedException e) {
	            e.printStackTrace();
	        }

		// Clean up the temporary file
		tempScriptFile.delete();

		return (exitCode);
	}	
	
	   private static void readStream(BufferedReader reader, String streamType) {
	        String line;
	        try {
	            while ((line = reader.readLine()) != null) {
	                System.out.println(streamType + ": " + line);
	            }
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
	    }
}
