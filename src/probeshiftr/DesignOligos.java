package probeshiftr;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.io.FilenameUtils;

public class DesignOligos {

	public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {

		Options options = new Options();

		options.addOption("t", "target-fasta", true, "fasta containing target sequence(s)");
		options.addOption("o", "output-dir", true, "directory to store oligo designs and temporary files (e.g. BLAT results)");
		options.addOption("d", "database-fasta", true, "database sequence file for BLAT searches (e.g. genome.fasta)");
		options.addOption("D", "database-dir", true, "directory containing fasta files for BLAT searches");
		options.addOption("l", "oligo-length", true, "length of antisense oligos");
		options.addOption("s", "get-sense", false, "create sense oligos (Default: false)");
		
		options.addOption("match", "min-match", true, "minimal matches for BLAT searches (Default: 1)");
		options.addOption("score", "min-score", true, "mininmal score for BLAT searches (Default: 10)");
//		options.addOption("i", "identity-cutoff", true, "filter oligos based on matches in BLAT search");
		
		options.addOption("bf", "bases2filter", true, "bases to filter for polybases in oligos (Default: ACGT)");
		options.addOption("rmf", "repeat-masking-format", true, "are repeats masked wth lower cases (Default: true)");
		options.addOption("nm", "n-repeats", true, "relative freuqnecies of N repeats in oligos (Default: 0.001)");
		options.addOption("r", "max-repeats", true, "maximal percetage of repeats (Default: 0.07)");
		options.addOption("n", "include-n", true, "include N in repeat filtering (Default: true)");
		
		options.addOption("pbl", "polybase-length", true, "relative length of polybases within oligo (Default: 0.8)");
		options.addOption("rbl", "repeat-length-polybases", true, "length of polybase repreats within sequences (Default: 15)");
		options.addOption("ncpf", "no-complexity-filter", false, "set of complexity filter (Default: set on)");
		options.addOption("cw", "check-within-target", false, "check for homologies within target sequence (Default: false)");
		options.addOption("irk", "include-repetitive-kmers", false, "include repetitive kmers (Default: false)");
		options.addOption("f", "fragment-size", true, "if gaps between oligos are longer than the fragment size then oligos with too many off-targets are used to fill these gaps (Default: 200)");
		
		options.addOption("g", "gtf-files", true, "gtf/gff annotation files to consider protein-coding exonic off-target regions (separated by semicolon without spaces e.g. anno1.gtf;anno2.gtf)");
		options.addOption("p", "processes", true, "number of parallel threads/processes to run BLAT comparisons");
		options.addOption("blat", "blat-path", true, "path to BLAT executable (Default: assumed to be in the environmental variable PATH)");
		options.addOption("rscript", "rscript-path", true, "path to Rscript executable (Default: assumed to be in the environmental variable PATH)");
		options.addOption("log", "log-file", false, "write log file to evaluate sequence complexity filters (Default: false)");
		options.addOption("v", "verbose", false, "Print debugging output (Default: false)");
		options.addOption("h", "help", false, "print this message");

		CommandLineParser parser = new DefaultParser();
		HelpFormatter formatter = new HelpFormatter();
		
		String toolName = "java -jar ProbeSHIFTR.jar";

		// oligo length and parallelization
		int oligo_length = 0;
		int threads = 1;
		boolean do_sense = false;
		boolean check_within = false;
		boolean include_rep_kmers = false;
		boolean write_log = false;
		boolean verbose = false;
		
		// BLAT prarameters
		int minMatch = 1; 
		int minScore = 10; 
		
		// Parameter repeat filtering
		boolean repeatLowerMask = true;
		boolean repeatIncludeN = true;
		double maxRepeatPercentage = 0.07;
		double maxRepeatNPercentage = 0.001;
//		double complexityCutoff = 0.08;
		
		// Parameter complexity filtering
		boolean no_complexity_filter = false;
		
		// Parameter poly base filtering
		String bases2Filter = "ACTG";
		double relativePolyBaseLength = 0.8;
		int lengthPolyBaseRepeats = 15; 
		int fragment_size = 200;
		
		String target_fasta_file = null; 
		String genome_database_fasta_file = null; 
		
		String blatPath = null; 
		String super_output_dir = null; 
		String genome_database_fasta_dir = null; 

		
		String path_to_Rscript = "/callableRscript.R";
		
		// only necessary if filtering on transcripts not on subsequences in the genome
		String transcript_annotation_files = null; 
		
		String plot_dir = null; 
		String fasta_out_dir = null; 
		String blat_dir = null; 
		
		String R_bin_path = null;
		
		ProgramFinder prgFinder = new ProgramFinder();
		
		try {
			CommandLine line = parser.parse(options, args);

			if (line.hasOption("h")) {
				formatter.printHelp(toolName, options, true);
				System.exit(0);
			}

			if (line.hasOption("t")) {
				
				target_fasta_file = line.getOptionValue("t");

				File file = new File(target_fasta_file);
				
				if (!file.exists()) {
					System.err.println("Could not find input target fasta file "+ target_fasta_file +"!");
					formatter.printHelp(toolName, options, true);
					System.exit(0);
				}
			}else {
				
				System.err.println("Please define a fasta file containing RNA target sequences.");
				formatter.printHelp(toolName, options, true);
				System.exit(0);
				
			}
			
			if(line.hasOption("l")) {
				
				oligo_length = Integer.parseInt(line.getOptionValue("l"));
			
			}else {
				
				System.err.println("Please define the length of the desired oligo sequences.");
				formatter.printHelp(toolName, options, true);
				System.exit(0);
			}
			
			if (line.hasOption("blat")) {
				
				blatPath = line.getOptionValue("blat");

				File file = new File(blatPath);
				
				if (!file.exists()) {
					System.err.println("Could not find BLAT application at: \n"+ blatPath +"\n Please check your BLAT path.");
					formatter.printHelp(toolName, options, true);
					System.exit(0);
				}
			}else {
				
				blatPath = prgFinder.findPath("blat");
				
				System.out.println("Using default BLAT path: " + blatPath);
				
				File file = new File(blatPath);
				
				if (!file.exists()) {
					System.err.println("Could not find BLAT application at: \n"+ blatPath +"\n Please check your BLAT path.");
					formatter.printHelp(toolName, options, true);
					System.exit(0);
				}
				
			}
			
			if (line.hasOption("rscript")) {
				
				R_bin_path = line.getOptionValue("rscript");

				File file = new File(path_to_Rscript);
				
				if (!file.exists()) {
					System.err.println("Could not find Rscript application at: \n"+ R_bin_path +"\n Please check your Rscript path.");
					formatter.printHelp(toolName, options, true);
					System.exit(0);
				}
			}else {
				
				R_bin_path = prgFinder.findPath("Rscript"); //"/usr/local/bin/Rscript";
				System.out.println("Using default Rscript path: " + R_bin_path);
				
				File file = new File(R_bin_path);
				
				if (!file.exists()) {
					System.err.println("Could not find Rscript application at: \n"+ R_bin_path +"\n Please check your Rscript path.");
					formatter.printHelp(toolName, options, true);
					System.exit(0);
				}
				
			}
			
			
	
			if (line.hasOption("d") & line.hasOption("D")) {
				
				System.err.println("Please define either a fasta file containing the fasta file (-d) for BLAT comparisons or define a directory that contains all fasta files (-D).");
				formatter.printHelp(toolName, options, true);
				System.exit(0);
				
			}
			
			if (line.hasOption("d") || line.hasOption("D")) {
				
				if (line.hasOption("d")) {
					genome_database_fasta_file = line.getOptionValue("d");
					genome_database_fasta_dir = null;
				}
				
				if (line.hasOption("D")) {
					genome_database_fasta_dir = line.getOptionValue("D");
					genome_database_fasta_file = null;
				}

			}else {
				
				System.err.println("Please define either a fasta file containing the fasta file (-d) for BLAT comparisons or define a directory that contains all fasta files (-D).");
				formatter.printHelp(toolName, options, true);
				System.exit(0);
				
			}

			if (line.hasOption("o")) {
				
				super_output_dir = line.getOptionValue("o");
				
				File outputDir = new File(super_output_dir);
				
				
				plot_dir = super_output_dir + File.separatorChar + "plots_final_oligo_designs";
				fasta_out_dir = super_output_dir + File.separatorChar + "final_oligo_designs";
				blat_dir = super_output_dir + File.separatorChar + "BLAT_results";

				if (!outputDir.exists()){
					outputDir.mkdirs();
				}
				
				File file_blat_dir = new File(blat_dir);
				
				if (!file_blat_dir.exists()){
					file_blat_dir.mkdirs();
				}
				
			}else {
				
				System.err.println("Please define a name for the output directory!");
				formatter.printHelp(toolName, options, true);
				System.exit(0);
				
			}
			
			if (line.hasOption("score")) {
				
				minScore = Integer.parseInt(line.getOptionValue("score"));
			}
			
			if (line.hasOption("match")) {
				
				minMatch = Integer.parseInt(line.getOptionValue("match"));
			}
			
			
			if (line.hasOption("bf")) {
				
				bases2Filter = line.getOptionValue("bf");
			}
			
			if (line.hasOption("rmf")) {
				
				repeatLowerMask = line.getOptionValue("rmf") == "lower" ? true : false;
				
			}
			
			if (line.hasOption("nm")) {
				
				maxRepeatNPercentage = Double.parseDouble(line.getOptionValue("nm"));
				
			}
			
			if (line.hasOption("r")) {
				
				maxRepeatPercentage = Double.parseDouble(line.getOptionValue("r"));
				
			}
			
			if (line.hasOption("f")) {
				
				fragment_size = Integer.parseInt(line.getOptionValue("f"));
				
			}
			
			if (line.hasOption("pbl")) {
				
				lengthPolyBaseRepeats = Integer.parseInt(line.getOptionValue("pbl"));
				
			}
			
			if (line.hasOption("rbl")) {
				
				relativePolyBaseLength = Double.parseDouble(line.getOptionValue("rbl"));
				
			}
			
			if (line.hasOption("ncpf")) {
				
				no_complexity_filter = true;
				
			}
			
			if (line.hasOption("n")) {
				
				repeatIncludeN = Boolean.getBoolean(line.getOptionValue("n"));
				
			}
			

			if (line.hasOption("s")) {
				
				do_sense = true;
				
			}


			if (line.hasOption("log")) {
				
				write_log = true;
				
			}

			if (line.hasOption("cw")) {
				
				check_within = true;
				
			}
			
			if (line.hasOption("irk")) {
				
				include_rep_kmers = true;
				
			}

			if (line.hasOption("v")) {
				
				verbose = true;
			}
			
			if (line.hasOption("g")) {
				
				transcript_annotation_files = line.getOptionValue("g");
				
				String[] tx_anno_file_arr = transcript_annotation_files.split(";");
				
				for(String single_anno_filename: tx_anno_file_arr) {
					
					File single_anno_file = new File(single_anno_filename);
					
					if(!single_anno_file.exists()) {
						
						System.err.println("Could not find " + single_anno_filename +
								 		   ". Please check the path of your annotation files.");
						formatter.printHelp(toolName, options, true);
						System.exit(0);
					}
				}
				
			}
			
			if(line.hasOption("p")) {
				
				threads = Integer.parseInt(line.getOptionValue("p"));
			}
			
		}catch (ParseException exp) {

			System.err.println("Parsing failed.  Reason: " + exp.getMessage());
			formatter.printHelp(toolName, options);
		}
		
		if(repeatLowerMask) {
			
			System.out.println("Repeats are masked by lower cases.");
		
		}else {

			System.out.println("Repeats are masked by upper cases.");
		
		}

		String oligo_fasta_file = super_output_dir + File.separatorChar + "Unfiltered_oligos_antisense_" + oligo_length + "nt.fa";
		String log_file = super_output_dir + File.separatorChar + "All_kmer_oligos_antisense_" + oligo_length + "nt.log";
				
		HashMap<String, Sequence> seqHash = SeqIO.readFastaGenom(target_fasta_file);
	
		RepeatFilter rf = new RepeatFilter(maxRepeatPercentage, repeatLowerMask, repeatIncludeN);
		RepeatFilter rfn = new RepeatFilter(maxRepeatNPercentage, false, repeatIncludeN);
		ComplexityFilter cf = new ComplexityFilter();
		PolybaseFilter pf = new PolybaseFilter(bases2Filter, relativePolyBaseLength, lengthPolyBaseRepeats);
//		LowComplexityFilter lf = new LowComplexityFilter(genomeFileName, threads, kmerLength, complexityCutoff);
	
		int[] filtered = new int[seqHash.keySet().size()];
		
		ArrayList<HashMap<Integer, Sequence>> probeSet = new ArrayList<HashMap<Integer, Sequence>>();
		
		BufferedWriter out_log = null;
		
		if(write_log) {
			
			out_log = new BufferedWriter(new FileWriter(log_file));
			out_log.write("seq\tkmer\tstart\tend\tRepeatFilter\tRepeatFilterN\tComplexFilter\tPolyBaseFilter\n");
		}
		
		
		int set_i = 0;
		for(String seq: seqHash.keySet()) {
			
			filtered[set_i] = 0;
			String faSeq = seqHash.get(seq).getSeq();
			String id = seqHash.get(seq).getId();
			
			System.out.println("Processing target: " + id);
			
			HashMap<String, Integer> kmerHash = Oligo.kmerHash(seqHash.get(seq), oligo_length);
			HashMap<Integer, Sequence> kmerSet = new HashMap<Integer, Sequence>();

			for(String kmer: kmerHash.keySet()) {
				
				boolean repeat = rf.filterOligo(kmer);
				boolean repeatN = rfn.filterOligo(kmer);
				boolean poly = pf.filterOligo(kmer);
//				boolean lowComplex = lf.filterOligo(kmer);
				boolean complex = false;
				boolean filterOligo = false;
				
				// Use complexity filter
				if(! no_complexity_filter) {

					complex = cf.filterOligo(kmer);
				}
				
				filterOligo = filterOligo || repeat || repeatN || complex || poly;
				
				if(write_log) {
					
					int startIndex = -1;
					String startIndices = null;
					String endIndices = null;
					
					ArrayList<Integer> oligoOccurencesArray = new ArrayList<Integer>();
					
					while((startIndex = faSeq.indexOf(kmer, startIndex +1)) != -1) {
						oligoOccurencesArray.add(startIndex);
					}
					
					final int oligo_length2add = oligo_length;

					startIndices = oligoOccurencesArray.stream().map(String::valueOf).collect(Collectors.joining(";"));
					endIndices = oligoOccurencesArray.stream().map(n -> String.valueOf(n + oligo_length2add - 1)).collect(Collectors.joining(";"));
						
					out_log.write(seq + "\t" + kmer + "\t" + startIndices + "\t" + endIndices + "\t"+ repeat + "\t" + repeatN + "\t" + complex + "\t" + poly + "\n");
				}
				
				if(filterOligo) {

					filtered[set_i]++;
					
				}else {
					
					int startIndex = -1 ;
					
					if(include_rep_kmers) {
						
						while((startIndex = faSeq.indexOf(kmer, startIndex +1)) != -1) {

							Sequence kmerSeq = new Sequence(id, kmer);

							String final_oligo = null; 
							String kmerId = null;

							if(do_sense) {

								final_oligo = kmerSeq.getSeq();
								kmerId = "Probe_Sense_" + id + "_" + startIndex + "_" + (startIndex + oligo_length - 1);

								
							}else {
								
								final_oligo = kmerSeq.getReverseComplement().getSeq();
								kmerId = "Probe_Antisense_" + id + "_" + startIndex + "_" + (startIndex + oligo_length - 1);

							}
							
							/*if(write_log) {
								out_log.write(seq + "\t" + kmer + "\t" + startIndex + "\t" + (startIndex + oligo_length - 1) + "\t"+ repeat + "\t" + repeatN + "\t" + complex + "\t" + poly + "\n");
							}*/
							
							Sequence kmerProbe = new Sequence(kmerId, final_oligo);
							
							kmerSet.put(startIndex, kmerProbe);
						}
					}else {
						
						startIndex = faSeq.indexOf(kmer, 0);

						Sequence kmerSeq = new Sequence(id, kmer);

						String final_oligo = null; 
						String kmerId = null;

						if(do_sense) {

							final_oligo = kmerSeq.getSeq();
							kmerId = "Probe_Sense_" + id + "_" + startIndex + "_" + (startIndex + oligo_length - 1);
								
						}else {
								
							final_oligo = kmerSeq.getReverseComplement().getSeq();
							kmerId = "Probe_Antisense_" + id + "_" + startIndex + "_" + (startIndex + oligo_length - 1);

						}
						
						/*
						 if(write_log) {
							out_log.write(seq + "\t" + kmer + "\t" + startIndex + "\t" + (startIndex + oligo_length - 1) + "\t"+ repeat + "\t" + repeatN + "\t" + complex + "\t" + poly + "\n");
						}
						*/
						
						Sequence kmerProbe = new Sequence(kmerId, final_oligo);
							
						kmerSet.put(startIndex, kmerProbe);
					}
					
					
				}
				
			}
			probeSet.add(kmerSet);
			
//			System.out.println(filtered[set_i] + ": Oligos filtered for " + id);
			set_i++;
		}
		
		if(write_log) {
			out_log.close();
		}
		/**
		 *  Write probe sequences into fasta file 
		 * */
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(oligo_fasta_file));
		
		for(int i = 0; i<probeSet.size(); i++) {
			
			HashMap<Integer, Sequence> kmerSet = probeSet.get(i);
			
			Integer[] startPositions = new Integer[kmerSet.size()];
					
			kmerSet.keySet().toArray(startPositions);
			
			Arrays.sort(startPositions);
			
			
			for(Integer start: startPositions) {
				
				out.write(">" + kmerSet.get(start).getId()+ "\n");
				out.write(kmerSet.get(start).getSeq()+ "\n");
			}

		}
		out.close();
		
		System.out.println("Starting BLAT alignment using " + threads + " parallel processes...");
		
		// create a pool of threads
		ExecutorService threadPool = Executors.newFixedThreadPool(threads);
		 
		// submit jobs to be executing by the pool
		ArrayList<Future<Void>> submitList = new ArrayList<Future<Void>>();
		 
		if(genome_database_fasta_file != null) {
			
			HashMap<String, Sequence> seqMap = SeqIO.readFastaGenom(genome_database_fasta_file);
			
			String[] seqIDs = new String[seqMap.size()];
			
			seqMap.keySet().toArray(seqIDs);
			
			String genomePath = FilenameUtils.getFullPath(genome_database_fasta_file);
			
			for (int i = 0; i < seqIDs.length; i++) {
				
				System.out.print("Process " + (i+1) + " of " + seqIDs.length + " sequences\r");
				
				FilenameUtils.getFullPath(genome_database_fasta_file);
				
				String genomeTmpFile = genomePath + File.separator + seqIDs[i] + ".fa";
				
				SeqUtils.writeFASTA(seqMap.get(seqIDs[i]), genomeTmpFile, 80);
				
				String outputFile = blat_dir + File.separator + seqIDs[i] + "_BLAT_oligo.pls";
				
				PerformBLAT blProcess = new PerformBLAT(blatPath, minMatch, minScore, genomeTmpFile, oligo_fasta_file, outputFile);
			
				Callable<Void> worker = new BLATCallable(blProcess);
				Future<Void> submit = threadPool.submit(worker);
				submitList.add(submit);
			}

		}else {
			
			File folder = new File(genome_database_fasta_dir);
			File[] files = folder.listFiles();
			
			for (int i = 0; i < files.length; i++) {
				
				System.out.print("Read "+(i+1)+" of "+files.length+" files\r");
				
				String chrFile = files[i].getPath();
				
				FilenameUtils.getBaseName(chrFile);
				
				String outputFile = blat_dir + File.separator + FilenameUtils.getBaseName(chrFile) + "_BLAT_oligo.pls";
				
				PerformBLAT blProcess = new PerformBLAT(blatPath, minMatch, minScore, chrFile, oligo_fasta_file, outputFile);
			
				Callable<Void> worker = new BLATCallable(blProcess);
				Future<Void> submit = threadPool.submit(worker);
				submitList.add(submit);
				
				//System.out.println(outputFile);
			}
			
		}
		
		try {
			 System.out.println("Perform BLAT searches...");
			 
			 int nbJobs = submitList.size();
			 int queueSize = nbJobs;
			 
			 do{
				 queueSize = ((ThreadPoolExecutor) threadPool).getQueue().size();
				 System.out.print((nbJobs-queueSize)+" of "+nbJobs+"\r");
				 Thread.sleep(10000);
						 
			 }while(queueSize > 0);
			 
			  threadPool.shutdown();

			 try {
		       // Wait for all tasks to complete or timeout after a specified duration
				 while (!threadPool.awaitTermination(10, TimeUnit.SECONDS)) {
		                System.out.println("Waiting for all tasks to finish...");
		            }
		     } catch (InterruptedException e) {
		    	 // Re-interrupt the current thread if interrupted during awaitTermination
		    	 threadPool.shutdownNow();
		         Thread.currentThread().interrupt();
		     }
			  
			 Thread.sleep(10000);
			 System.out.println("Finished!");
			 System.out.println("BLAT results stored at: " + blat_dir);

			  
			} catch (InterruptedException e) {
			  System.out.println(e.getMessage());
			} 
		
		
		
		
		
		/** Run Rscript to calculate finale oligo sets and create analysis plots */
		
		Rinterface blProcess = new Rinterface(R_bin_path, path_to_Rscript, blat_dir, fasta_out_dir, plot_dir, 
											  oligo_fasta_file, oligo_length, target_fasta_file, check_within, 
											  fragment_size, transcript_annotation_files);
		blProcess.runR(verbose);
	}

}

class BLATCallable implements Callable<Void> {

	private PerformBLAT process;
	
	BLATCallable(PerformBLAT process){
		super();
		this.process = process;
	}
	
	@Override
	public Void call() throws Exception{
		
		try{
			process.runBLAT();			
			
		}catch(Exception e){	

			System.out.println(e.getStackTrace());
			System.out.println(e.getMessage());
			
			StringWriter sw = new StringWriter();
			PrintWriter pw = new PrintWriter(sw);
			e.printStackTrace(pw);
			String sStackTrace = sw.toString(); // stack trace as a string
			System.out.println(sStackTrace);
		}
		return null;
	}
}
