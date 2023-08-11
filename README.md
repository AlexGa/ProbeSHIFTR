# ProbeSHIFTR
Oligosequence probe design for Selective RNase-H-mediated interactome framing for target RNA regions

## Requirements

- JavaSE-14
- R (tested on R-4.2.1)
- required R packages:
	- readr
	- data.table
	- rtracklayer
	- GenomicRanges
	- IRanges
	- tidyr
	- seqinr
	- ggplot2
- install BLAT
- java, R, BLAT and all required  packages can be found in `environmental.yaml`

## Installation

Use conda to create an environment for ProbeSHIFTR:

`conda env create -f environment.yaml`

Acitvate ProbeSHIFTR environment:

`conda activate probeshiftr`


## Usage

```
java -jar ProbeSHIFTR.jar [-bf <arg>] [-blat <arg>] [-d <arg>] [-D <arg>] 
	   					  [-g <arg>] [-h] [-l <arg>] [-match <arg>] [-n <arg>] 
	   					  [-nm <arg>] [-o <arg>] [-pbl <arg>] [-r <arg>] 
	   					  [-rbl <arg>] [-rmf <arg>] [-rscript <arg>] 
	   					  [-score <arg>] [-t <arg>]
 -bf,--bases2filter <arg>               bases to filter for polybases in
                                        oligos (Default: ACGT)
 -blat,--blat-path <arg>                path to BLAT executable (Default:
                                        assumed to be in the environmental
                                        variable PATH)
 -d,--database-fasta <arg>              database sequence file for BLAT
                                        searches (e.g. genome.fasta)
 -D,--database-dir <arg>                directory containing fasta files
                                        for BLAT searches
 -g,--gtf-file <arg>                    gtf/gff file containing regions to
                                        ignore if designed oligo shows
                                        overlap
 -h,--help                              print this message
 -l,--oligo-length <arg>                length of antisense oligos
 -match,--min-match <arg>               minimal matches for BLAT searches
                                        (Default: 1)
 -n,--include-n <arg>                   include N in repeat filtering
                                        (Default: true)
 -nm,--n-repeats <arg>                  relative freuqnecies of N repeats
                                        in oligos (Default: 0.001)
 -o,--output-dir <arg>                  directory to store oligo designs
                                        and temporary files (e.g. BLAT
                                        results)
 -pbl,--polybase-length <arg>           relative length of polybases
                                        within oligo (Default: 0.8)
 -r,--max-repeats <arg>                 maximal percetage of repeats
                                        (Default: 0.07)
 -rbl,--repeat-length-polybases <arg>   length of polybase repreats within
                                        sequences (Default: 15)
 -rmf,--repeat-masking-format <arg>     repeats masked with upper or lower
                                        cases or no masking in data
                                        (lower, upper; Default: lower)
 -rscript,--rscript-path <arg>          path to Rscript executable
                                        (Default: assumed to be in the
                                        environmental variable PATH)
 -score,--min-score <arg>               mininmal score for BLAT searches
                                        (Default: 10)
 -t,--target-fasta <arg>                fasta containing target sequence
                                        for antisense oligos
```
