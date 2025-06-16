suppressWarnings({
  library(dplyr)
  library(readr)
  library(data.table)
})

gaps_path <- function(x, x_end){
  
  x_e <- c(-1, x[,2])
  x_s <- c(x[,1], x_end)
  
  return(x_s - x_e - 1)
}

genesFromOligos <- function(x) {
  unlist(lapply(strsplit(x,"_"),function(l) paste(l[3:(length(l)-2)],collapse="_")))
}

color_theme_schema <- ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom",
                 legend.title = ggplot2::element_blank(),
                 legend.text = ggplot2::element_text(size = 11),
                 axis.text = ggplot2::element_text(colour = "black", size = 12),
                 axis.title = ggplot2::element_text(colour = "black", size = 14, face = "bold"),
                 strip.background = ggplot2::element_rect(fill = "white"),
                 strip.text = ggplot2::element_text(size = 12, color = "black"),
                 plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

args <- commandArgs(trailingOnly = TRUE)

blat_dir <- args[1] # directory containing the BLAT output created from java framework
outDir <- args[2] # output directory containing the final oligo sets
plot_dir <- args[3] # directory containing the analysis plots for different thresholds
oligo_fasta_file <- args[4] # fasta file containing oligo sequences from java program
oligo_length <- as.numeric(args[5]) # length of oligos -> inferred from sequence length
target_fasta_file <- args[6] # fasta file containing the RNA target sequences

# boolean value to check within 
check_within_target = as.logical(args[7])

# fragment size 
fragment_size <- as.numeric(args[8])

anno_files <- NULL # check anno_file

if(length(args) > 8){
  
  anno_files <- unlist(strsplit(x = args[9], split = ";")) # gtf annotation compressed or uncompressed
  
}

# Not yet implemented
#--------------------
# ignored_homologous_regions <- c()

if(length(args) > 9){
  
  ignored_homologous_regions <- args[10] # vector -> put into list format, each element represents one target sequence
  
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
      
      
      if(any(grepl(pattern = "gene_biotype", x = single_anno %>% as_tibble() %>% colnames()))){
        
        print(paste0("Using only protein_coding transcripts based on gene_biotype entry from ", single_anno_file))
        single_anno <- as.data.frame(single_anno[single_anno$type %in% "exon" & single_anno$gene_biotype %in% "protein_coding"][,c("type", "transcript_id")])
        
      }else{
        
        print(paste0("Using all transcripts provided by ", single_anno_file))
        single_anno <- as.data.frame(single_anno[single_anno$type %in% "exon"][,c("type", "transcript_id")])
        
      }
      single_anno
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
    
    readr::read_delim(file.path(blat_dir, f), 
                      skip = 6, 
                      col_names = F, 
                      show_col_types = FALSE, 
                      na = "")
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
  dplyr::mutate(query.chr = gsub(x = query.name, pattern = "Probe_(Antisense|Sense)_(.+)_([0-9]+)_([0-9]+)$", repl = "\\2"),
                query.start = as.numeric(gsub(x = query.name, pattern = "Probe_(Antisense|Sense)_(.+)_([0-9]+)_([0-9]+)$", 
                                              repl = "\\3")),
                query.end = as.numeric(gsub(x = query.name, pattern = "Probe_(Antisense|Sense)_(.+)_([0-9]+)_([0-9]+)$", 
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
  # using only mature spliced sequences
  
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
  blat_df <- data.frame(id = 1:nrow(blat_df), blat_df)
  rm(ixn)
  
}else{
  
  blat_df <- probe_anno %>% 
    as.data.frame() %>%
    unique()
}


## Not yet implemented

# ignored_blat_results <- c()
# 
# for(target_gene in names(ignored_homologous_regions)){
#   
#   homo_genes <- paste0(c(target_gene, ignored_homologous_regions[[target_gene]]), collapse = "|")
#   
#   ignored_regions <- anno %>% as.data.frame() %>% filter(grepl(pattern = homo_genes, 
#                                                                x = gene_name, ignore.case = T), 
#                                                          type %in% "exon") %>% 
#     GenomicRanges::makeGRangesFromDataFrame()
#   
#   blat_df_sub_query <- blat_df %>% filter(grepl(pattern = target_gene, x = oligo_name))
#   
#   blat_df_sub_query_g_ranges <- GenomicRanges::makeGRangesFromDataFrame(with(blat_df_sub_query, 
#                                                                              data.frame(chr = seqnames, 
#                                                                                         strand = strand,
#                                                                                         start = start, 
#                                                                                         end = end,
#                                                                                         id = id)),keep.extra.columns = T)
#   
#   ixn <- IRanges::findOverlaps(query = blat_df_sub_query_g_ranges, 
#                                subject = ignored_regions, ignore.strand = TRUE)
#   
#   # remove unspecific within ignored (multiple in the region)
#   
#   homo_specific_ignored <- blat_df_sub_query_g_ranges[unique(ixn@from)]
#   # homo_specific_ignored <- homo_specific_ignored %>% group_by(oligo_name) %>% 
#   #                          mutate(N = n()) %>% 
#   #                         dplyr::select(oligo_name, N, gene) %>% filter(N <= 1)
#   
#   ignored_blat_results <- c(ignored_blat_results, homo_specific_ignored$id)
# }
# 
# blat_df <- blat_df %>% filter(!id %in% ignored_blat_results)
## end of -- Not yet implemented

measures_tbl <- data.frame(matrix(ncol = 11, nr = 0))
colnames(measures_tbl) <- c("gene_name", "probe_id", "mismatches", "filtered", 
                            "coverage", "#gaps", "cum_gaps", "min_gaps", "max_gaps", 
                            "mean_gaps", "sd_gaps")

row_count <- 1

regexp_for_split <- c()

if(grepl(x = attr(oligo_seqs, "name")[1], pattern = "Probe_Antisense")){
  
  regexp_for_split <- c("Probe_Antisense_(.+)_[0-9]+_[0-9]+$", "Probe_Antisense_(.+)_([0-9]+)_([0-9]+)$")
  
}else if(grepl(x = attr(oligo_seqs, "name")[1], pattern = "Probe_Sense")){
  
  regexp_for_split <- c("Probe_Sense_(.+)_[0-9]+_[0-9]+$", "Probe_Sense_(.+)_([0-9]+)_([0-9]+)$")
  
}else{
  
  stop("Please check the FASTA IDs of your target sequences.")
  
}

rm(anno)
rm(probe_anno)


all_oligos_df <- data.frame(oligo_name = names(oligo_seqs)) %>% as_tibble() %>% 
  tidyr::extract(col = oligo_name, regex = regexp_for_split[2], into = c("target_name", "oligo.start", "oligo.end"), remove = F) %>% 
  mutate(oligo.start = as.numeric(oligo.start), 
         oligo.end = as.numeric(oligo.end))


off_targets_all_thresholds <- list()
match_seqs <- c()

if(oligo_length < 50){
  match_seqs <- seq(from = min(blat_df$matches), to = oligo_length, length.out = 5)
}else{
  match_seqs <- seq(from = min(blat_df$matches), to = oligo_length, by = 10)
}


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
  
  blat_tbl_sub <- blat_tbl_sub %>% filter(gene_name %in% target_seqs$target_name)
  
  # Create porbeset with unique coverage along the targets
  target_lists <- lapply(split(blat_tbl_sub, f = blat_tbl_sub$gene_name), function(x) {rownames(x) <- c(); x})
  probe_sets <- list()
  sets_for_gap_fillings <- list()
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
      off_targets_per_oligo <- off_target_df %>% filter(oligo_name %in% (tab %>% pull(oligo_name))) %>% as.data.frame()
      
      if(dim(off_targets_per_oligo)[1] == 0){
        off_targets_per_oligo <- data.frame("oligo_name" = tab$oligo_name, 
                                            "gene" = tab$gene,
                                            "n" = 1)
      }
      
      stats <- c(mean_off = sum(off_targets_per_oligo$n)/nrow(tab), 
                 max_off = max(off_targets_per_oligo$n),
                 oligos_w_off = nrow(off_targets_per_oligo))
      stats
    })
    
    off_targets_per_set_set[[j]] <- lapply(probe_sets[[j]], function(tab){
      off_target_df %>% dplyr::rename(gene_name = gene) %>%
        filter(oligo_name %in% (tab %>% pull(oligo_name))) %>% 
        as.data.frame()
    })
    
    sets_for_gap_fillings[[j]] <- list()
    
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
      
      target_boundaries <- c(0, target_seqs %>% filter(target_name == target_gene_name) %>% pull(target_length))
      
      gaps_vec <- gaps_path(x = single_probe_set[,3:4], x_end = target_boundaries[2])
      
      gaps_df <- data.frame(start = c(0, single_probe_set[,4]), end = c(single_probe_set[,3], target_boundaries[2]-1)) %>% mutate(gap = end - start)
      
      # find gaps longer than fragment size and search for oligos to fill them
      
      gaps2fill <- gaps_df %>% filter(gap > fragment_size)
      
      
      oligos2fill <- unlist(apply(gaps2fill, 1, function(gap_elem){
        
        pos2fill <- all_oligos_df %>% filter(target_name %in% target_gene_name, oligo.start >= as.numeric(gap_elem[1]), oligo.end <= as.numeric(gap_elem[2]))
        dists <- t(apply(as.matrix(pos2fill[,-c(1:2)]), 1, function(elem1){
          
          (as.numeric(elem1[2]) - as.numeric(pos2fill %>% pull(oligo.start)) +1)
          
        }))
        
        rownames(dists) <- colnames(dists) <- pos2fill$oligo_name
        oligo_set <- NULL
        
        if(length(dists) != 0){
          
          final_end <- -1
          row_i <- 1
          oligo_set <- rownames(dists)[row_i]
          
          while(final_end < gap_elem[2] - oligo_length - 1){
            
            is_dist_less_zero <- dists[row_i,] <= 0
            
            if(any(is_dist_less_zero)){
              
              row_i <- min(which(is_dist_less_zero))
              
            }else{
              
              # Check if there is an overlapping oligo sequence which could be used
              row_i_new <- which.min(dists[row_i,])
              
              # Check if only possible oligo that could be used was already picked up
              # if so we reached the end of possible oligo seqs and loop ends
              
              if(rownames(dists)[row_i] == rownames(dists)[row_i_new]){
                
                break
                
              }else{
                
                row_i <- row_i_new
                
              }
            }
            
            oligo_set <- c(oligo_set, rownames(dists)[row_i])
            
            final_end <- pos2fill %>% filter(oligo_name == rownames(dists)[row_i]) %>% pull(oligo.end)
          }
        }
        
        oligo_set
      }))
      
      sets_for_gap_fillings[[j]][[set_ind]] <<- oligos2fill
      
      coverage <- (as.numeric(max_match) * nrow(single_probe_set)) / target_boundaries[2]
      
      measures_tbl[row_count, ] <<- list(target_gene_name, counter_probe_set, match_threshold, n_filtered %>% filter(gene_name == target_gene_name) %>% pull(filtered), coverage,
                                         length(gaps_vec > 0), sum(gaps_vec[gaps_vec > 0]), min(gaps_vec), max(gaps_vec),
                                         mean(gaps_vec), sd(gaps_vec))
      row_count <<- row_count + 1
    })
    
    off_target_df_sets <- list()
    #sets_for_gap_fillings
    off_target_df_sets_with_gap_fillings <- list()
    
    probe_result_dir <- file.path(super_plot_dir, target_gene_name)
    
    for(probe_off_index in seq_along(probe_sets[[j]])){
      
      sub_plot_dir <- file.path(probe_result_dir, "plots")
      sub_fasta_dir <- file.path(probe_result_dir, "fasta")
      
      if(!file.exists(sub_plot_dir)){
        dir.create(sub_plot_dir, recursive = TRUE)
      }
      
      if(!file.exists(sub_fasta_dir)){
        dir.create(sub_fasta_dir, recursive = TRUE)
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
      
      
      fasta_file_name <- file.path(sub_fasta_dir, paste0(target_gene_name, "_Set_", probe_off_index, ".fa"))
      seqinr::write.fasta(sequences = oligo_seqs[single_probe_set_off_cov$oligo_name], 
                          file.out = fasta_file_name, 
                          names = names(oligo_seqs[single_probe_set_off_cov$oligo_name]))
      
      
      pos <- as.numeric(unlist(apply(single_probe_set_off_cov[,c(-1,-2)], 1, function(cont) rep(seq(from = cont[1], to = cont[2]), each = cont[3]))))
      
      target_boundaries <- c(0, target_seqs %>% 
                               filter(target_name == target_gene_name) %>% 
                               pull(target_length)
      )
      
      title_text <- paste0(target_gene_name, " [", dim(single_probe_set_off_cov)[1]," oligos]")
      
      subtitle_text <- paste0("expected off-targets: ", format(off_targets_per_set[[j]][[probe_off_index]][1], digits = 2),
                              " | max off-targets: ", off_targets_per_set[[j]][[probe_off_index]][2],
                              " | #oligos with off-targets: ", off_targets_per_set[[j]][[probe_off_index]][3])
      
      off_target_plot <- as.data.frame(table(factor(pos,levels = seq(target_boundaries[1],target_boundaries[2], 1)))) %>% 
        ggplot2::ggplot(ggplot2::aes(x = as.numeric(Var1), y =  Freq)) +
        ggplot2::geom_bar(stat = "identity", fill = "red4", col = "red4") + 
        ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), name = "") +
        ggplot2::scale_y_continuous(expand = c(0, 0, 0.05, 0), name = "# off-targets") +
        ggplot2::labs(title = title_text,
                      subtitle = subtitle_text) +
        color_theme_schema +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       plot.margin = ggplot2::margin(b = -0.5,unit = "cm"))
      
      
      pos_cov <- as.numeric(unlist(apply(single_probe_set_off_cov[,c(-1,-2)], 1, function(cont) rep(seq(from = cont[1], to = cont[2]), each = ifelse(cont[3] == 0, 1, cont[3])))))
      
      coverage_plot <- as.data.frame(table(factor(pos_cov,levels = seq(target_boundaries[1],target_boundaries[2], 1)))) %>% dplyr::mutate(Freq = dplyr::if_else(Freq == 0, true = 0, false = 1)) %>% 
        ggplot2::ggplot(ggplot2::aes(x = as.numeric(Var1), y =  Freq)) +
        ggplot2::geom_bar(stat = "identity", fill = "black", col = "black") + 
        ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), name = "target sequence position", breaks = scales::pretty_breaks(n = 10)) +
        ggplot2::scale_y_continuous(expand = c(0, 0, 0, 0), name = "", breaks = c(0, 1), labels = c("", "")) +
        color_theme_schema +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       plot.margin = ggplot2::margin(t = 0))
      
      comb_plot <- cowplot::plot_grid(off_target_plot, coverage_plot, rel_heights = c(0.85, 0.15), ncol = 1, align = "v")
      ggplot2::ggsave(plot = comb_plot, filename = file.path(sub_plot_dir, paste0(target_gene_name, "_", probe_off_index, "_coverage_and_off_targets.pdf")), width = 10, height = 4, device = "pdf")
      
      
      if(length(sets_for_gap_fillings[[j]]) == length(probe_sets[[j]])){
        
        single_probe_set_off_cov_with_filled_oligos <- rbind(all_oligos_df %>% 
                                                               filter(oligo_name %in% sets_for_gap_fillings[[j]][[probe_off_index]]) %>% 
                                                               dplyr::rename(gene_name = target_name), probe_sets[[j]][[probe_off_index]]) %>% 
          dplyr::arrange(oligo.start) %>%
          dplyr::inner_join(off_target_df, by = c("oligo_name" = "oligo_name")) %>%
          dplyr::mutate(n = replace(n, is.na(n), 0)) %>% dplyr::select(-gene)
        
        
        fasta_file_name <- file.path(sub_fasta_dir, paste0(target_gene_name, "_Set_", probe_off_index, "_filled_gaps_with_off_target_oligos.fa"))
        seqinr::write.fasta(sequences = oligo_seqs[single_probe_set_off_cov_with_filled_oligos$oligo_name], 
                            file.out = fasta_file_name, 
                            names = names(oligo_seqs[single_probe_set_off_cov_with_filled_oligos$oligo_name]))
        
        oligo_set_seqs_with_filled_oligos <- sapply(oligo_seqs[single_probe_set_off_cov_with_filled_oligos$oligo_name], function(seq_elem) seqinr::c2s(unlist(seq_elem)))
        
        
        off_target_df_sets_with_gap_fillings[[probe_off_index]] <- tibble(oligo_name = names(oligo_set_seqs_with_filled_oligos), sequence = oligo_set_seqs_with_filled_oligos) %>% 
          inner_join(single_probe_set_off_cov_with_filled_oligos, by = "oligo_name") %>% 
          mutate(Set_ID = probe_off_index) %>%
          dplyr::rename(target_seq = gene_name,
                        `#off targets` = n) %>%
          dplyr::select(-oligo.start, -oligo.end) %>% mutate("Filled to cover fragment size" = dplyr::if_else(oligo_name %in% sets_for_gap_fillings[[j]][[probe_off_index]], true = "true", false = "false")) %>%
          dplyr::relocate(target_seq, Set_ID, oligo_name, `Filled to cover fragment size`, `#off targets`, sequence)
        
        pos_filled <- as.numeric(unlist(apply(single_probe_set_off_cov_with_filled_oligos %>% 
                                                filter(oligo_name %in% sets_for_gap_fillings[[j]][[probe_off_index]]) %>% 
                                                dplyr::select(oligo.start, oligo.end, n), 1, function(cont) rep(seq(from = cont[1], to = cont[2]), each = cont[3]))))
        
        off_target_plot_with_filled <- as.data.frame(table(factor(pos,levels = seq(target_boundaries[1],target_boundaries[2], 1)))) %>% 
          ggplot2::ggplot(ggplot2::aes(x = as.numeric(Var1), y =  Freq)) +
          ggplot2::geom_bar(stat = "identity", fill = "red4", col = "red4") + 
          ggplot2::geom_bar(data = as.data.frame(table(factor(pos_filled,levels = seq(target_boundaries[1],target_boundaries[2], 1)))),
                            stat = "identity", fill = "darkgreen", col = "darkgreen", alpha = 0.1) + 
          ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), name = "") +
          ggplot2::scale_y_continuous(expand = c(0, 0, 0.05, 0), name = "# off-targets") +
          ggplot2::labs(title = title_text,
                        subtitle = subtitle_text) +
          color_theme_schema +
          ggplot2::theme(panel.grid = ggplot2::element_blank(),
                         axis.text.x = ggplot2::element_blank(),
                         axis.ticks.x = ggplot2::element_blank(),
                         plot.margin = ggplot2::margin(b = -0.5,unit = "cm"))
        
        comb_plot <- cowplot::plot_grid(off_target_plot_with_filled, coverage_plot, rel_heights = c(0.85, 0.15), ncol = 1, align = "v")
        ggplot2::ggsave(plot = comb_plot, filename = file.path(sub_plot_dir, paste0(target_gene_name, "_", probe_off_index, "_coverage_and_off_targets_with_filled_gaps.pdf")), width = 10, height = 4, device = "pdf")
        
      }
      
      
    }
    
    names(off_target_df_sets) <- paste0("Set_", seq_along(probe_sets[[j]]))
    writexl::write_xlsx(off_target_df_sets, path = file.path(probe_result_dir, paste0(target_gene_name, "_sets.xlsx")))
    
    if(length(off_target_df_sets_with_gap_fillings) != 0){
      names(off_target_df_sets_with_gap_fillings) <- paste0("Set_", seq_along(probe_sets[[j]]))
      writexl::write_xlsx(off_target_df_sets_with_gap_fillings, path = file.path(probe_result_dir, paste0(target_gene_name, "_sets_including_filled_gaps.xlsx")))
    }
    
  }
  
  names(off_targets_per_set) <- names(target_lists)
  off_targets_all_thresholds[[match_cut_ind]] <- off_targets_per_set
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
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::xlab("minimal allowed BLAT match length") +
  ggplot2::ylab("Mean off targets") + 
  ggplot2::scale_colour_brewer(palette = "Paired") +
  color_theme_schema  +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_min_match_length_vs_mean_off_targets.pdf"), device = "pdf", dpi = 300)


off_targets_stats  %>% group_by(gene_name, mismatches) %>% summarise(max_off = max(max_off)) %>%
  ggplot2::ggplot(ggplot2::aes(x = as.numeric(mismatches), y = max_off, col = gene_name)) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::xlab("minimal allowed BLAT match length") +
  ggplot2::ylab("Maximum number off targets over all sets") +
  ggplot2::scale_colour_brewer(palette = "Paired") +
  color_theme_schema +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_min_match_length_vs_max_off_targets.pdf"), device = "pdf", dpi = 300)


off_targets_stats %>% group_by(gene_name, mismatches) %>% summarise(oligos_w_off = mean(oligos_w_off)) %>%
  ggplot2::ggplot(ggplot2::aes(x = as.numeric(mismatches), y = oligos_w_off, col = gene_name)) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::xlab("minimal allowed BLAT match length") +
  ggplot2::ylab("Mean Number oligos with off targets per set") +
  ggplot2::scale_colour_brewer(palette = "Paired") +
  color_theme_schema +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_min_match_length_vs_max_off_targets.pdf"), device = "pdf", dpi = 300)

measures_tbl %>% dplyr::select(mismatches, filtered, gene_name) %>% unique() %>%
  ggplot2::ggplot(ggplot2::aes(x = mismatches, y = filtered, col = gene_name)) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::xlab("minimal allowed BLAT match length") +
  ggplot2::ylab("#Filtered oligos") +
  ggplot2::scale_colour_brewer(palette = "Paired") +
  color_theme_schema +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_min_match_length_vs_filtered_oligos.pdf"), device = "pdf", dpi = 300)

measures_tbl %>% dplyr::select(mismatches, coverage, gene_name) %>% unique() %>%
  ggplot2::ggplot(ggplot2::aes(x = mismatches, y = coverage, col = gene_name)) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::xlab("minimal allowed BLAT match length") +
  ggplot2::ylab("Transcript coverage") + 
  ggplot2::coord_cartesian(ylim = c(0,1)) +
  ggplot2::scale_colour_brewer(palette = "Paired") +
  color_theme_schema +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_min_match_length_vs_coverage_transcript.pdf"), device = "pdf", dpi = 300)

measures_tbl %>% dplyr::select(mismatches, cum_gaps, gene_name) %>% unique() %>%
  ggplot2::ggplot(ggplot2::aes(x = mismatches, y = cum_gaps, col = gene_name)) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line(linewidth = 0.75) +
  # ggplot2::facet_wrap(~gene_name, ncol = 2) +
  ggplot2::xlab("minimal allowed BLAT match length") +
  ggplot2::ylab("Cumulative gaps between oligos") +
  ggplot2::scale_colour_brewer(palette = "Paired") +
  color_theme_schema +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)) 
ggplot2::ggsave(filename = file.path(plot_dir, "comparison_min_match_length_vs_cumulative_gap.pdf"), device = "pdf", dpi = 300)
