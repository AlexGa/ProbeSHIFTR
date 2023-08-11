library(dplyr)
library(readr)
library(data.table)

gaps_path <- function(x, x_end){
  
  x_e <- c(-1, x[,2])
  x_s <- c(x[,1], x_end)
  
  return(x_s - x_e - 1)
}

genesFromOligos <- function(x) {
  unlist(lapply(strsplit(x,"_"),function(l) paste(l[3:(length(l)-2)],collapse="_")))
}

args <- commandArgs(trailingOnly = TRUE)

# args[1] <- "/Users/alg22/LRIB/ProbeDesigner/data/hg38/Homo_sapiens.GRCh38.106.gtf.gz"
# args[2] <- "/Users/alg22/LRIB/Jens/Oligos_Cdifficile/OligoDesign/TEST_targets_antisense_50/BLAT_results"
# args[3] <- "/Users/alg22/LRIB/Jens/Oligos_Cdifficile/OligoDesign/TEST_targets_antisense_50/final_oligo_designs"
# args[4] <- "/Users/alg22/LRIB/Jens/Oligos_Cdifficile/OligoDesign/TEST_targets_antisense_50/plots_final_oligo_designs"
# args[5] <- "/Users/alg22/LRIB/Jens/Oligos_Cdifficile/OligoDesign/TEST_targets_antisense_50/Unfiltered_oligos_antisense_50nt.fa"
# args[6] <- 50
# args[7] <- "/Users/alg22/LRIB/Jens/Oligos_Cdifficile/RaiA_long_nc70.fa"

anno_file <- args[1] # gtf annotation compressed or uncompressed
blat_dir <- args[2] # directory containing the BLAT output created from java framework
outDir <- args[3] # output directory containing the final oligo sets
plot_dir <- args[4] # directory containing the analysis plots for different thresholds

oligo_fasta_file <- args[5] # fasta file containing oligo sequences from java program
oligo_length <- as.numeric(args[6]) # length of oligos -> infer from sequence length!
target_fasta_file <- args[7] # fasta file containing the RNA target sequences

# Not yet implemented
ignored_homologous_regions <- c()

if(length(args) > 7){
  
  ignored_homologous_regions <- args[8] # vector -> put into list format, each element represents one target sequence
  
}

# Read input data 
#----------------

# GTF annotation file
anno <- NULL
if(!file.exists(anno_file)){
  
  anno <- rtracklayer::import(anno_file)
  print("Using transcript annotation file and filter based on mature spliced transcripts...")
  
}
# Read target RNA sequences from fasta 
target_seqs <- seqinr::read.fasta(target_fasta_file)
target_seqs <- data.frame(target_name = attr(target_seqs, "name"), 
                          target_length = unlist(lapply(target_seqs, length)))


# Read oligo sequences 
oligo_seqs <- seqinr::read.fasta(oligo_fasta_file)

# Files fro BLAT comparisons 
files <- dir(blat_dir)

if(length(files) > 0){
  blat_df_big <- do.call("rbind", sapply(files, function(f){
    readr::read_delim(file.path(blat_dir, f), skip = 6, col_names = F, show_col_types = FALSE)
  }, simplify = FALSE))
  
  colnames(blat_df_big) <- c("matches","mismatch", "repmatch", "Ns", 
                             paste0(rep(c("query.", "target."), each = 2), rep(c("gap.count", "gap.bases"), 2)), "strand",
                             paste0(rep(c("query.", "target."), each = 4), rep(c("name", "size", "start", "end"), 2)), 
                             "block.count", "block.sizes", "query.starts", "target.starts")
}else{
  
  blat_df_big <- NULL
}


if(!file.exists(outDir)){
  dir.create(outDir)
}

if(!file.exists(plot_dir)){
  dir.create(plot_dir)
}

blat_dt <- blat_df_big %>% as.data.table()
blat_dt <- blat_dt[,list(target.name, target.start, target.end, strand, query.name), ]

probe_anno <- GenomicRanges::makeGRangesFromDataFrame(
  with(blat_df_big %>% dplyr::select(target.name, target.start, target.end, query.name, matches, strand),
       data.frame(chr = target.name,
                  start = target.start,
                  end = target.end,
                  strand = strand,
                  oligo_name = query.name,
                  matches = matches)),
  keep.extra.columns = TRUE)

blat_df <- c()

if(!is.null(anno)){
  # concentrate only on mature spliced sequences
  anno_target <- anno[anno$type %in% "exon"]
  
  ixn <- IRanges::findOverlaps(query = probe_anno, 
                               subject = anno_target)
  
  blat_df <- probe_anno[unique(ixn@from)] %>% as.data.frame() %>% unique()
  
}else{
  
  blat_df <- probe_anno %>% as.data.frame() %>% unique()
}

blat_df <- data.frame(id = 1:nrow(blat_df), blat_df)

ignored_blat_results <- c()

for(target_gene in names(ignored_homologous_regions)){
  
  homo_genes <- paste0(c(target_gene, ignored_homologous_regions[[target_gene]]), collapse = "|")
  
  ignored_regions <- anno %>% as.data.frame() %>% filter(grepl(pattern = homo_genes, 
                                                               x = gene_name, ignore.case = T), 
                                                         type %in% "exon") %>% 
    GenomicRanges::makeGRangesFromDataFrame()
  
  blat_df_sub_query <- blat_df %>% filter(grepl(pattern = target_gene, x = oligo_name))
  
  blat_df_sub_query_g_ranges <- GenomicRanges::makeGRangesFromDataFrame(with(blat_df_sub_query, 
                                                                             data.frame(chr = seqnames, 
                                                                                        strand = strand,
                                                                                        start = start, 
                                                                                        end = end,
                                                                                        id = id)),keep.extra.columns = T)
  
  ixn <- IRanges::findOverlaps(query = blat_df_sub_query_g_ranges, 
                               subject = ignored_regions, ignore.strand = TRUE)
  
  # remove unspecific within ignored (multiple in the region)
  
  homo_specific_ignored <- blat_df_sub_query_g_ranges[unique(ixn@from)]
  # homo_specific_ignored <- homo_specific_ignored %>% group_by(oligo_name) %>% 
  #                          mutate(N = n()) %>% 
  #                         dplyr::select(oligo_name, N, gene) %>% filter(N <= 1)
  
  ignored_blat_results <- c(ignored_blat_results, homo_specific_ignored$id)
}

blat_df <- blat_df %>% filter(!id %in% ignored_blat_results)

measures_tbl <- data.frame(matrix(ncol = 11, nr = 0))
colnames(measures_tbl) <- c("gene_name", "probe_id", "mismatches", "filtered", "coverage", "#gaps", "cum_gaps", "min_gaps", "max_gaps", "mean_gaps", "sd_gaps")
row_count <- 1

for(match_threshold in c(seq(0, oligo_length, 5), oligo_length)){
  
  print(paste0("Caculate oligo design based on threshold: ", match_threshold))
  
  max_match <- oligo_length # query length
  super_plot_dir <- paste0(outDir, "/match_threshold_",match_threshold)
  super_fasta_dir <- paste0(outDir, "/match_threshold_",match_threshold)
  
  if(!file.exists(super_fasta_dir)){
    dir.create(super_fasta_dir, recursive = T)
  }
  
  if(!file.exists(super_plot_dir)){
    dir.create(super_plot_dir, recursive = T)
  }
  
  blat_df_sub <- blat_df %>% dplyr::filter(matches >= match_threshold) %>%
    dplyr::mutate(gene = genesFromOligos(oligo_name))
  
  oligo_names <- names(oligo_seqs)
  
  n_seqs_per_target <- oligo_names %>% as_tibble() %>% tidyr::extract(col = value, regex = "Probe_Antisense_(.*)_.*?_.*", into = "gene_name", remove = F) %>% group_by(gene_name) %>% count()
  
  if(dim(blat_df_sub)[1] != 0){
    
    blat_removed <- blat_df_sub %>% group_by(oligo_name, gene) %>% 
      summarise(N = n()) %>% filter(N > 1)
    
    oligo_names <- oligo_names[!oligo_names %in% blat_removed$oligo_name]
  }
  
  n_filtered_per_target <- oligo_names %>% as_tibble() %>% tidyr::extract(col = value, regex = "Probe_Antisense_(.*)_.*?_.*", into = "gene_name", remove = F) %>% group_by(gene_name) %>% count()
  
  n_filtered <- n_seqs_per_target %>% dplyr::inner_join(n_filtered_per_target, by = "gene_name") %>% mutate(filtered = n.x - n.y)
  
  # n_filtered <- length(names(oligo_seqs)) - length(oligo_names)
  
  pos <- do.call("rbind", lapply(strsplit(oligo_names, "_"), 
                                 function(cont) as.numeric(cont[c(length(cont)-1, length(cont))])))
  
  blat_tbl_sub <- data.frame(oligo_name = oligo_names, 
                             oligo.start = pos[,1], oligo.end = pos[,2]) %>% 
    tidyr::extract(col = oligo_name, regex = "Probe_Antisense_(.*)_.*?_.*", into = "gene_name", remove = F)
  
  # Create porbeset with unique coverage along the targets
  target_lists <- lapply(split(blat_tbl_sub, f = blat_tbl_sub$gene_name), function(x) {rownames(x) <- c(); x})
  probe_sets <- list()
  probe_names <- c()
  
  for(j in seq_along(target_lists)){
    
    print(names(target_lists)[j])
    
    blat_tbl <- target_lists[[j]]
    
    probe_names[j] <- names(target_lists)[j]
    
    dists <- t(apply(blat_tbl[,-c(1:2)], 1, function(elem1){
      
      (as.numeric(elem1[1]) - as.numeric(blat_tbl[,4]))
      
    }))
    
    colnames(dists) <- rownames(dists) <- blat_tbl$oligo_name
    
    probe_set_paths <- list()
    
    start_oligo <- blat_tbl[1,]
    x_start <- 1
    
    while(x_start <=  nrow(blat_tbl) & blat_tbl[x_start,]$oligo.start < start_oligo$oligo.end){
      
      oligo_path_indices <- c()
      
      i <- x_start
      oligo_path_indices <- i
      
      while(i < nrow(blat_tbl)){
        
        indices_dist <- which(dists[,i] >= 1)
        
        if(length(indices_dist) > 0){
          
          i <- min(indices_dist)
          oligo_path_indices <- c(oligo_path_indices, i)
          
        }else{
          
          break;
        }
      }
      probe_set_paths[[x_start]] <- as.numeric(rownames(blat_tbl)[oligo_path_indices])
      x_start <- x_start + 1
    }
    
    length_paths <- unlist(lapply(probe_set_paths, length))
    
    probe_sets[[j]] <- lapply(probe_set_paths[which(length_paths == max(length_paths))], function(probe_set_idx){
      
      blat_tbl[probe_set_idx,]
      
    })
    
    # Plot probesets
    #---------------
    
    counter_probe_set <- 0
    target_gene_name <- unique(probe_sets[[j]][[1]]$gene_name)
    
    lapply(probe_sets[[j]], function(single_probe_set){
      
      counter_probe_set <<- counter_probe_set + 1
      
      pos <- as.numeric(unlist(apply(single_probe_set[,c(-1,-2)], 1, function(cont) seq(from = cont[1], to = cont[2]))))
      
      plot_dir <- file.path(super_plot_dir, target_gene_name)
      
      if(!file.exists(plot_dir)){
        dir.create(plot_dir)
      }
      
      target_boundaries <- c(0, target_seqs %>% filter(grepl(x = target_name, pattern = target_gene_name, fixed = T)) %>% pull(target_length))
      # target_boundaries <- c(0, target_seqs %>% filter(target_name == target_gene_name) %>% pull(target_length))
      
      pdf(file.path(plot_dir, paste0(target_gene_name, "_", counter_probe_set, ".pdf")), width = 10, height = 4)
      barplot(beside = T, table(factor(pos, levels = seq(target_boundaries[1],target_boundaries[2], 1))), col = "black",
              main = paste0("Oligo alignment ", target_gene_name, "(", dim(single_probe_set)[1],")")) 
      dev.off()
      
      gaps_vec <- gaps_path(x = single_probe_set[,3:4], x_end = target_boundaries[2])
      
      coverage <- (as.numeric(max_match) * nrow(single_probe_set)) / target_boundaries[2]
      
      measures_tbl[row_count, ] <<- list(target_gene_name, counter_probe_set, match_threshold, n_filtered %>% filter(gene_name == target_gene_name) %>% pull(filtered), coverage, 
                                         length(gaps_vec > 0), sum(gaps_vec[gaps_vec > 0]), min(gaps_vec), max(gaps_vec), 
                                         mean(gaps_vec), sd(gaps_vec))
      row_count <<- row_count + 1
      
      
    })
  }
  
  # Write probe sets into fasta files
  counter_probe_set <- 0
  fasta_probe_sets <- lapply(probe_sets, function(target){
    
    counter_probe_set <<- counter_probe_set + 1
    
    probe_counter <<- 1
    lapply(target, function(probe_set){
      
      fasta_dir <- file.path(super_fasta_dir, probe_names[counter_probe_set])
      
      if(!file.exists(fasta_dir)){
        dir.create(fasta_dir, recursive = T)
      }
      
      fasta_file_name <- file.path(fasta_dir, paste0(probe_names[counter_probe_set], "_Set_", probe_counter, ".fa"))
      fasta_idx <- which(attr(oligo_seqs, "name") %in% probe_set$oligo_name)
      seqinr::write.fasta(sequences = oligo_seqs[fasta_idx], file.out = fasta_file_name, names = attr(oligo_seqs[fasta_idx], "name"))
      
      probe_counter <<- probe_counter + 1 
    })
  })
}
# 
measures_tbl %>% select(mismatches, filtered, gene_name) %>% unique() %>%
  ggplot2::ggplot(ggplot2::aes(x = mismatches, y = filtered, col = gene_name)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::theme_bw() +
  ggplot2::xlab("allowed BLAT mismatches") +
  ggplot2::ylab("#Filtered oligos")
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_mismatches_vs_filtered_oligos.pdf"), device = "pdf")

measures_tbl %>% select(mismatches, coverage, gene_name) %>% unique() %>%
  ggplot2::ggplot(ggplot2::aes(x = mismatches, y = coverage, col = gene_name)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::theme_bw() +
  ggplot2::xlab("allowed BLAT mismatches") +
  ggplot2::ylab("Transcript coverage") + 
  ggplot2::coord_cartesian(ylim = c(0,1))
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_mismatches_vs_coverage_transcript.pdf"), device = "pdf")

measures_tbl %>% select(mismatches, cum_gaps, gene_name) %>% unique() %>%
  ggplot2::ggplot(ggplot2::aes(x = mismatches, y = cum_gaps, col = gene_name)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~gene_name) +
  ggplot2::xlab("allowed BLAT mismatches") +
  ggplot2::ylab("Cumulative gaps between oligos") 
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_mismatches_vs_cumulative_gap.pdf"), device = "pdf")
