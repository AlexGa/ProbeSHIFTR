---
editor_options: 
  markdown: 
    wrap: 72
---

# ProbeSHIFTR

Oligonucleotide probe designer for **S**elective RNase-**H**-mediated
**I**nteractome **F**raming for **T**arget **R**NA Regions (SHIFTR)

### Design of oligo sequences for the publication

**SHIFTR enables the unbiased identification of proteins bound to
specific RNA regions in live cells**.

Jens Aydin, Alexander Gabel, Sebastian Zielinski, Sabina Ganskih, Nora
Schmidt, Christina R. Hartigan, Monica Schenone, Steven A. Carr, and
Mathias Munschauer.

## Requirements

-   JavaSE-14
-   R (tested on R-4.2.1)
-   required R packages:
    -   readr
    -   data.table
    -   rtracklayer
    -   GenomicRanges
    -   IRanges
    -   tidyr
    -   seqinr
    -   ggplot2
-   install BLAT
-   java, R, BLAT and all required packages can be found in
    `environmental.yaml`

## Installation

Use conda to create an environment for ProbeSHIFTR:

`conda env create -f environment.yaml`

Acitvate ProbeSHIFTR environment:

`conda activate probeshiftr`

## Examples

### Create 25nt long oligos for 7sk and u1

---

#### Before running the examples

Decompress the sequence file in seq_db.

```         
gunzip test/seq_db/fasta/chr1.fa.gz
```

------------------------------------------------------------------------

#### 1.) Perform BLAT comparisons against a single sequence

```         
java -jar jar/ProbeSHIFTR.jar -d test/seq_db/fasta/chr1.fa 
                              -t test/test_input/test_targets.fa 
                              -l 25
                              -o test/test_output
```

------------------------------------------------------------------------

#### 2.) Perform BLAT comparisons against a directory of uncompressed fasta files

```         
java -jar jar/ProbeSHIFTR.jar -D test/seq_db/fasta/ 
                              -t test/test_input/test_targets.fa 
                              -l 25 
                              -o test/test_output
```

------------------------------------------------------------------------

**Note:** In the upper two commands there are no `u1` oligos left
because intergenic regions are considered as off-targets. When comparing
against large genomic sequences without any annotation of exonic
regions, off-targets are defined in intergenic regions which might be
too general for certain applications (target sequences of RNA-binding
proteins).

------------------------------------------------------------------------

#### 3.) Perform BLAT comparisons against exonic regions (recommended)

If the fasta files in the `seq_db` are genome sequences, such as
chr1 of HG38, you can add annotation data in `gtf` format. Thus, the
filtering of off-targets takes only hits into account that overlap with
exonic regions. If the gtf contains `gene_biotype "protein_coding";`,
off-targets are only defined by oligos showing similarities to
protein-coding exons.

```         
java -jar jar/ProbeSHIFTR.jar -D test/seq_db/fasta/ 
                              -t test/test_input/test_targets.fa 
                              -l 25 
                              -g test/seq_dq/gtf/hg38.ncbiRefSeq_chr1.gtf.gz
                              -o test/test_output
```

------------------------------------------------------------------------

## Usage

```         
usage: java -jar ProbeSHIFTR.jar [-bf <arg>] [-blat <arg>] [-cw] [-d
       <arg>] [-D <arg>] [-f <arg>] [-g <arg>] [-h] [-irk] [-l <arg>]
       [-log] [-match <arg>] [-n <arg>] [-ncpf] [-nm <arg>] [-o <arg>] [-p
       <arg>] [-pbl <arg>] [-r <arg>] [-rbl <arg>] [-rmf <arg>] [-rscript
       <arg>] [-s] [-score <arg>] [-t <arg>] [-v]
 -bf,--bases2filter <arg>               bases to filter for polybases in
                                        oligos (Default: ACGT)
 -blat,--blat-path <arg>                path to BLAT executable (Default:
                                        assumed to be in the environmental
                                        variable PATH)
 -cw,--check-within-target              check for homologies within target
                                        sequence (Default: false)
 -d,--database-fasta <arg>              database sequence file for BLAT
                                        searches (e.g. genome.fasta)
 -D,--database-dir <arg>                directory containing fasta files
                                        for BLAT searches
 -f,--fragment-size <arg>               if gaps between oligos are longer
                                        than the fragment size then oligos
                                        with too many off-targets are used
                                        to fill these gaps (Default: 200)
 -g,--gtf-files <arg>                   gtf/gff annotation files to
                                        consider protein-coding exonic
                                        off-target regions (separated by
                                        semicolon without spaces e.g.
                                        anno1.gtf;anno2.gtf)
 -h,--help                              print this message
 -irk,--include-repetitive-kmers        include repetitive kmers (Default:
                                        false)
 -l,--oligo-length <arg>                length of antisense oligos
 -log,--log-file                        write log file to evaluate
                                        sequence complexity filters
                                        (Default: false)
 -match,--min-match <arg>               minimal matches for BLAT searches
                                        (Default: 1)
 -n,--include-n <arg>                   include N in repeat filtering
                                        (Default: true)
 -ncpf,--no-complexity-filter           set of complexity filter (Default:
                                        set on)
 -nm,--n-repeats <arg>                  relative freuqnecies of N repeats
                                        in oligos (Default: 0.001)
 -o,--output-dir <arg>                  directory to store oligo designs
                                        and temporary files (e.g. BLAT
                                        results)
 -p,--processes <arg>                   number of parallel
                                        threads/processes to run BLAT
                                        comparisons
 -pbl,--polybase-length <arg>           relative length of polybases
                                        within oligo (Default: 0.8)
 -r,--max-repeats <arg>                 maximal percetage of repeats
                                        (Default: 0.07)
 -rbl,--repeat-length-polybases <arg>   length of polybase repreats within
                                        sequences (Default: 15)
 -rmf,--repeat-masking-format <arg>     are repeats masked wth lower cases
                                        (Default: true)
 -rscript,--rscript-path <arg>          path to Rscript executable
                                        (Default: assumed to be in the
                                        environmental variable PATH)
 -s,--get-sense                         create sense oligos (Default:
                                        false)
 -score,--min-score <arg>               mininmal score for BLAT searches
                                        (Default: 10)
 -t,--target-fasta <arg>                fasta containing target
                                        sequence(s)
 -v,--verbose                           Print debugging output (Default:
                                        false)
```
