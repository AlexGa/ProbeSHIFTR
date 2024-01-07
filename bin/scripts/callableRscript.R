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

blat_dir <- args[1] # directory containing the BLAT output created from java framework
outDir <- args[2] # output directory containing the final oligo sets
plot_dir <- args[3] # directory containing the analysis plots for different thresholds

oligo_fasta_file <- args[4] # fasta file containing oligo sequences from java program
oligo_length <- as.numeric(args[5]) # length of oligos -> infer from sequence length!
target_fasta_file <- args[6] # fasta file containing the RNA target sequences


blat_dir <- "/Users/alg22/LRIB/Cooperation/CastelloLab/ProbeDesign/OligoDesign_RVFV_targets_antisense_90nt/BLAT_results"
outDir <- "/Users/alg22/LRIB/Cooperation/CastelloLab/ProbeDesign/OligoDesign_RVFV_targets_antisense_90nt/final_oligo_designs"
plot_dir <- "/Users/alg22/LRIB/Cooperation/CastelloLab/ProbeDesign/OligoDesign_RVFV_targets_antisense_90nt/plots_final_oligo_designs"
oligo_fasta_file <- "/Users/alg22/LRIB/Cooperation/CastelloLab/ProbeDesign/OligoDesign_RVFV_targets_antisense_90nt/Unfiltered_oligos_antisense_90nt.fa"
oligo_length <- 90
target_fasta_file <- "/Users/alg22/LRIB/Cooperation/CastelloLab/ProbeDesign/references/RVFV_targets.fa"
check_within_target <- as.logical("true")
anno_files <- unlist(strsplit(x = "/Users/alg22/LRIB/ProbeDesigner/data/hg38/Homo_sapiens.GRCh38.106.gtf.gz;/Users/alg22/LRIB/Cooperation/CastelloLab/ProbeDesign/OligoDesign_references/Aedes_aegypti_lvpagwg.AaegL5.57_newSeq_Prefix.gtf.gz", split = ";"))


# boolean value to check within 
check_within_target = as.logical(args[7])

anno_files <- NULL # check anno_file

if(length(args) > 7){
  
  anno_files <- unlist(strsplit(x = args[8], split = ";")) # gtf annotation compressed or uncompressed
  
}

# Not yet implemented
ignored_homologous_regions <- c()

if(length(args) > 8){
  
  ignored_homologous_regions <- args[9] # vector -> put into list format, each element represents one target sequence
  
}

# Read input data 
#----------------

# GTF annotation file
anno <- NULL
if(!is.null(anno_files)){
  
  if(all(sapply(anno_files, file.exists))){
    
    print("Using transcript annotation file and filter based on mature spliced transcripts...")
    
    anno_list <- sapply(anno_files, function(single_anno_file){
      
      if(!file.exists(single_anno_file)){
        
        stop(paste0("The annotation file ", single_anno_file, " does not exist. Please check the path."))
        
      }
      
      print(paste0("Extracting transcript information from ",  single_anno_file, " and filter based on mature spliced transcripts..."))
      single_anno <- rtracklayer::import(single_anno_file)
      as.data.frame(single_anno[single_anno$type %in% "exon" & single_anno$gene_biotype %in% "protein_coding"][,c("type", "transcript_id")])
    }, simplify = F)
    
    anno_df <- do.call("rbind", anno_list)
    rownames(anno_df) <- 1:nrow(anno_df)
    
    anno <- GenomicRanges::makeGRangesFromDataFrame(anno_df, keep.extra.columns = T)
    
    rm(anno_df)
    rm(anno_list)
  }
}

# Read target RNA sequences from fasta 
target_seqs <- seqinr::read.fasta(target_fasta_file)
target_seqs <- data.frame(target_name = attr(target_seqs, "name"), 
                          target_length = unlist(lapply(target_seqs, length)))

target_anno <- GenomicRanges::makeGRangesFromDataFrame(
  data.frame(
    chr = target_seqs$target_name,
    start = 1,
    end = target_seqs$target_length,
    strand = ".",
    transcript_id = target_seqs$target_name,
    type = "exon"
  ),
  keep.extra.columns = TRUE
)

# Read oligo sequences 
oligo_seqs <- seqinr::read.fasta(oligo_fasta_file)

# Files fro BLAT comparisons 
files <- dir(blat_dir)

if(length(files) > 0){
  blat_df_big <- do.call("rbind", sapply(files, function(f){
    readr::read_delim(file.path(blat_dir, f), skip = 6, col_names = F, show_col_types = FALSE, na = "")
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

# Remove self-matches
blat_df_big <- blat_df_big %>% 
               dplyr::mutate(query.chr = gsub(x = query.name, pattern = "Probe_(Antisense|Sense)_(.*)_(.*)?_(.*)", repl = "\\2"),
                             query.start = as.numeric(gsub(x = query.name, pattern = "Probe_(Antisense|Sense)_(.*)_(.*)?_(.*)", 
                                                           repl = "\\3")),
                             query.end = as.numeric(gsub(x = query.name, pattern = "Probe_(Antisense|Sense)_(.*)_(.*)?_(.*)", 
                                                         repl = "\\4")) + 1) %>%
                dplyr::filter(!(query.chr == target.name & query.start == target.start & query.end == target.end)) %>% 
                dplyr::select(-query.chr, -query.start, -query.end)

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
  
  if(check_within_target){
    
    anno_target <- GenomicRanges::makeGRangesFromDataFrame(
      rbind(as.data.frame(target_anno),
            as.data.frame(anno))
    )
  }else{
    
    anno_target <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(anno))
    
  }
  
  ixn <- IRanges::findOverlaps(query = probe_anno,
                               subject = anno_target)

  blat_df <- probe_anno[unique(ixn@from)] %>% as.data.frame() %>% unique()

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
colnames(measures_tbl) <- c("gene_name", "probe_id", "mismatches", "filtered", 
                            "coverage", "#gaps", "cum_gaps", "min_gaps", "max_gaps", 
                            "mean_gaps", "sd_gaps")

row_count <- 1

regexp_for_split <- c()

if(grepl(x = attr(oligo_seqs, "name")[1], pattern = "Probe_Antisense")){
  
  regexp_for_split <- "Probe_Antisense_(.*)_.*?_.*"
  
}else if(grepl(x = attr(oligo_seqs, "name")[1], pattern = "Probe_Sense")){
  
  regexp_for_split <- "Probe_Sense_(.*)_.*?_.*"
  
}else{
  
  stop("Please check the FASTA IDs of your target sequences.")
  
}

rm(anno)
rm(probe_anno)
rm(ixn)

off_targets_all_thresholds <- list()

match_seqs <- seq(from = min(blat_df$matches), to = oligo_length, length.out = 10)
match_seqs <- seq(from = min(blat_df$matches), to = oligo_length, by = 5)

oligo_names <- names(oligo_seqs)

for(match_cut_ind in seq_along(match_seqs)){
  
  match_threshold <- match_seqs[match_cut_ind]
  
  print(paste0("Caculate oligo design based on threshold: ", match_threshold))
  
  max_match <- oligo_length
  
  blat_df_sub <- blat_df %>% dplyr::filter(matches >= match_threshold) %>%
                             dplyr::mutate(gene = genesFromOligos(oligo_name))
  
  off_target_df <- blat_df %>%
                   dplyr::mutate(gene = genesFromOligos(oligo_name)) %>% 
                   group_by(oligo_name, gene) %>% 
                   dplyr::count() %>% arrange(-n)
  
  oligo_names <- names(oligo_seqs)
  
  n_seqs_per_target <- oligo_names %>% as_tibble() %>% 
                                       tidyr::extract(col = value, regex = regexp_for_split[1], 
                                                      into = "gene_name", remove = F) %>% 
                                       group_by(gene_name) %>%  dplyr::count()
  
  if(dim(blat_df_sub)[1] != 0){
    
    blat_removed <- blat_df_sub %>% 
                    group_by(oligo_name, gene) %>% 
                    summarise(N = n()) %>% filter(N > 1)
    
    oligo_names <- oligo_names[!oligo_names %in% blat_removed$oligo_name]
  }
  
  n_filtered_per_target <- oligo_names %>% as_tibble() %>% 
                           tidyr::extract(col = value, regex = regexp_for_split[1], into = "gene_name", remove = F) %>% 
                           group_by(gene_name) %>% dplyr::count()
  
  n_filtered <- n_seqs_per_target %>% dplyr::inner_join(n_filtered_per_target, by = "gene_name") %>% mutate(filtered = n.x - n.y)
  
  # n_filtered <- length(names(oligo_seqs)) - length(oligo_names)
  
  pos <- do.call("rbind", lapply(strsplit(oligo_names, "_"), 
                                 function(cont) as.numeric(cont[c(length(cont)-1, length(cont))])))
  
  blat_tbl_sub <- data.frame(oligo_name = oligo_names, 
                             oligo.start = pos[,1], oligo.end = pos[,2]) %>% 
    tidyr::extract(col = oligo_name, regex = regexp_for_split[1], into = "gene_name", remove = F)
  
  # Create porbeset with unique coverage along the targets
  target_lists <- lapply(split(blat_tbl_sub, f = blat_tbl_sub$gene_name), function(x) {rownames(x) <- c(); x})
  probe_sets <- list()
  off_targets_per_set <- list()
  off_targets_per_set_set <- list()
  
  probe_names <- c()
  
  if(length(target_lists) == 0){
    next;
  }
  
  super_plot_dir <- paste0(outDir, "/match_threshold_",match_threshold)
  super_fasta_dir <- paste0(outDir, "/match_threshold_",match_threshold)
  
  if(!file.exists(super_fasta_dir)){
    dir.create(super_fasta_dir, recursive = T)
  }
  
  if(!file.exists(super_plot_dir)){
    dir.create(super_plot_dir, recursive = T)
  }
  
  print("Start preparing paths")
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
    
    while(x_start <=  nrow(blat_tbl) & blat_tbl[x_start,]$oligo.start < start_oligo$oligo.end & length(probe_set_paths) < 10){
      
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
    
    off_targets_per_set[[j]] <- lapply(probe_sets[[j]], function(tab){
                                      off_targets_per_set <- off_target_df %>% filter(oligo_name %in% (tab %>% pull(oligo_name))) %>% as.data.frame()
                                      stats <- c(mean_off = sum(off_targets_per_set$n)/nrow(tab), 
                                        max_off = max(off_targets_per_set$n),
                                        oligos_w_off = nrow(off_targets_per_set))
                                      stats
                                })
    
    off_targets_per_set_set[[j]] <- lapply(probe_sets[[j]], function(tab){
                                    off_targets_per_set <- off_target_df %>% dplyr::rename(gene_name = gene) %>%
                                                           filter(oligo_name %in% (tab %>% pull(oligo_name))) %>% 
                                                           as.data.frame()
                                    off_targets_per_set
                                })
    # Plot probesets
    #---------------
    
    counter_probe_set <- 0
    target_gene_name <- unique(probe_sets[[j]][[1]]$gene_name)
    set_ind <- 0

    lapply(probe_sets[[j]], function(single_probe_set){

      counter_probe_set <<- counter_probe_set + 1
      set_ind <<- set_ind + 1

      pos <- as.numeric(unlist(apply(single_probe_set[,c(-1,-2)], 1, function(cont) seq(from = cont[1], to = cont[2]))))

      plot_dir <- file.path(super_plot_dir, target_gene_name)

      if(!file.exists(plot_dir)){
        dir.create(plot_dir)
      }

      target_boundaries <- c(0, target_seqs %>% filter(grepl(x = target_name, pattern = target_gene_name, fixed = T)) %>% pull(target_length))
      # target_boundaries <- c(0, target_seqs %>% filter(target_name == target_gene_name) %>% pull(target_length))

      pdf(file.path(plot_dir, paste0(target_gene_name, "_", counter_probe_set, ".pdf")), width = 10, height = 4)
      barplot(beside = T, table(factor(pos, levels = seq(target_boundaries[1],target_boundaries[2], 1))), col = "black",
              main = paste0("Oligo alignment ", target_gene_name, "(", dim(single_probe_set)[1],")\n",
                            paste0("expected off-targets: ", format(off_targets_per_set[[j]][[set_ind]][1], digits = 2),
                                   " | Max off-targets: ", off_targets_per_set[[j]][[set_ind]][2],
                                   " | #oligos with off-targets: ", off_targets_per_set[[j]][[set_ind]][3])

                            ))
      dev.off()

      gaps_vec <- gaps_path(x = single_probe_set[,3:4], x_end = target_boundaries[2])

      coverage <- (as.numeric(max_match) * nrow(single_probe_set)) / target_boundaries[2]

      measures_tbl[row_count, ] <<- list(target_gene_name, counter_probe_set, match_threshold, n_filtered %>% filter(gene_name == target_gene_name) %>% pull(filtered), coverage,
                                         length(gaps_vec > 0), sum(gaps_vec[gaps_vec > 0]), min(gaps_vec), max(gaps_vec),
                                         mean(gaps_vec), sd(gaps_vec))
      row_count <<- row_count + 1


    })
    ### HIER WEITERMACHEN -> coverage der off targets pro oligo ###
    
    off_target_df_sets <- list()
    
    for(probe_off_index in seq_along(probe_sets[[j]])){
      
      sub_plot_dir <- file.path(super_plot_dir, target_gene_name)
      
      if(!file.exists(sub_plot_dir)){
        dir.create(sub_plot_dir)
      }
      
      single_probe_set_off_cov <- probe_sets[[j]][[probe_off_index]] %>% 
                                  dplyr::full_join(off_targets_per_set_set[[j]][[probe_off_index]]) %>%
                                  dplyr::mutate(n = replace(n, is.na(n), 0))
      
      oligo_set_seqs <- sapply(oligo_seqs[single_probe_set_off_cov$oligo_name], function(seq_elem) seqinr::c2s(unlist(seq_elem)))
      
      off_target_df_sets[[probe_off_index]] <- tibble(oligo_name = names(oligo_set_seqs), sequence = oligo_set_seqs) %>% 
            inner_join(single_probe_set_off_cov, by = "oligo_name") %>% 
            mutate(Set_ID = probe_off_index) %>%
            dplyr::rename(target_seq = gene_name,
                          `#off targets` = n) %>%
            dplyr::select(-oligo.start, -oligo.end) %>% 
            dplyr::relocate(target_seq, Set_ID, oligo_name, `#off targets`, sequence)
      
      # pos <- as.numeric(unlist(apply(probe_sets[[j]][,c(-1,-2)], 1, function(cont) seq(from = cont[1], to = cont[2]))))
      
      pos <- as.numeric(unlist(apply(single_probe_set_off_cov[,c(-1,-2)], 1, function(cont) rep(seq(from = cont[1], to = cont[2]), each = cont[3]))))
      
      target_boundaries <- c(0, target_seqs %>% 
                                filter(grepl(x = target_name, pattern = target_gene_name, fixed = T)) %>% 
                                pull(target_length)
                             )
  
     pdf(file.path(sub_plot_dir, paste0(target_gene_name, "_", probe_off_index, "_with_off_target.pdf")), width = 10, height = 4)
     barplot(beside = T, 
             table(factor(pos,
                          levels = seq(target_boundaries[1],target_boundaries[2], 1))),
             border = "red4", add = F,
             main = paste0("Oligo alignment ", target_gene_name, "(", dim(single_probe_set_off_cov)[1],")\n",
                           paste0("expected off-targets: ", format(off_targets_per_set[[j]][[probe_off_index]][1], digits = 2),
                                  " | Max off-targets: ", off_targets_per_set[[j]][[probe_off_index]][2],
                                  " | #oligos with off-targets: ", off_targets_per_set[[j]][[probe_off_index]][3])
                           
             ), ylim = c(-1, max(table(factor(pos,
                                              levels = seq(target_boundaries[1],target_boundaries[2], 1)))))) 
     
     pos <- as.numeric(unlist(apply(single_probe_set_off_cov[,c(-1,-2)], 1, function(cont) seq(from = cont[1], to = cont[2]))))
     barplot(beside = T, 
             ifelse(table(factor(pos, levels = seq(target_boundaries[1],target_boundaries[2], 1)))  == 0, yes = 0, no = -1), 
             add = T) 
     dev.off()
     
    }
    
    names(off_target_df_sets) <- paste0("Set_", seq_along(probe_sets[[j]]))
    writexl::write_xlsx(off_target_df_sets, path = file.path(sub_plot_dir, paste0(target_gene_name, "_sets.xlsx")))
  }
  
  names(off_targets_per_set) <- names(target_lists)
  
  
  off_targets_all_thresholds[[match_cut_ind]] <- off_targets_per_set
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

names(off_targets_all_thresholds) <- match_seqs



off_targets_all_thresholds_list <- lapply(off_targets_all_thresholds, function(list_cutoff){

  cutoff_genes_list <- list()
  for(i in seq_along(list_cutoff)){
    cutoff_genes_list[[i]] <- do.call("rbind", list_cutoff[[i]]) %>% as_tibble() %>% mutate(gene_name = names(list_cutoff)[i])
  }
  cutoff_genes_list
})

off_targets_stats_list <- list()

for (i in seq_along(off_targets_all_thresholds_list)){
  off_targets_stats_list[[i]] <- do.call("rbind", off_targets_all_thresholds_list[[i]]) %>% as_tibble() %>% mutate(mismatches = names(off_targets_all_thresholds_list)[i])
}

off_targets_stats <- do.call("rbind", off_targets_stats_list)

off_targets_stats %>% group_by(gene_name, mismatches) %>% summarise(mean_off = mean(mean_off)) %>%
  ggplot2::ggplot(ggplot2::aes(x = as.numeric(mismatches), y = mean_off, col = gene_name)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::xlab("minimal allowed BLAT match length") +
  ggplot2::ylab("Mean off targets")
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_min_match_length_vs_mean_off_targets.pdf"), device = "pdf")


off_targets_stats  %>% group_by(gene_name, mismatches) %>% summarise(max_off = max(max_off)) %>%
  ggplot2::ggplot(ggplot2::aes(x = as.numeric(mismatches), y = max_off, col = gene_name)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::xlab("minimal allowed BLAT match length") +
  ggplot2::ylab("Maximum number off targets over all sets")
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_min_match_length_vs_max_off_targets.pdf"), device = "pdf")


off_targets_stats %>% group_by(gene_name, mismatches) %>% summarise(oligos_w_off = mean(oligos_w_off)) %>%
  ggplot2::ggplot(ggplot2::aes(x = as.numeric(mismatches), y = oligos_w_off, col = gene_name)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::xlab("minimal allowed BLAT match length") +
  ggplot2::ylab("Mean Number oligos with off targets per set")
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_min_match_length_vs_max_off_targets.pdf"), device = "pdf")

measures_tbl %>% dplyr::select(mismatches, filtered, gene_name) %>% unique() %>%
  ggplot2::ggplot(ggplot2::aes(x = mismatches, y = filtered, col = gene_name)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::theme_bw() +
  ggplot2::xlab("minimal allowed BLAT match length") +
  ggplot2::ylab("#Filtered oligos")
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_min_match_length_vs_filtered_oligos.pdf"), device = "pdf")

measures_tbl %>% dplyr::select(mismatches, coverage, gene_name) %>% unique() %>%
  ggplot2::ggplot(ggplot2::aes(x = mismatches, y = coverage, col = gene_name)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::theme_bw() +
  ggplot2::xlab("minimal allowed BLAT match length") +
  ggplot2::ylab("Transcript coverage") + 
  ggplot2::coord_cartesian(ylim = c(0,1))
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_min_match_length_vs_coverage_transcript.pdf"), device = "pdf")

measures_tbl %>% dplyr::select(mismatches, cum_gaps, gene_name) %>% unique() %>%
  ggplot2::ggplot(ggplot2::aes(x = mismatches, y = cum_gaps, col = gene_name)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~gene_name) +
  ggplot2::xlab("minimal allowed BLAT match length") +
  ggplot2::ylab("Cumulative gaps between oligos") 
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_min_match_length_vs_cumulative_gap.pdf"), device = "pdf")
