## With this script, we can obtain a fasta file from an abundance table. We will
## translate the leaf names (OTU IDs) to their 16S sequence [using rep_seqs].
## Then, for each sample and each leaf name, we include N entries for each 
## sequence in the output fasta file, N being the absolute abundance of that
## leaf in that sample

## Files
OTU_table_file    <- "tomate_23sep_otu_table.csv"
fasta_temp_folder <- "fasta_temp_folder"
rep_seqs <- "rep_set.fasta"
output <- "tomate_subset.fa"

## Parameters
remove_last_column <- TRUE
sep <- "\t"

#---------------------------------------------------------------------------

library("Biostrings")
library("stringr")

system(paste("mkdir", fasta_temp_folder))

# read table
a <- read.csv(OTU_table_file, sep=sep, row.names=1)
if (remove_last_column) {
  a <- a[-ncol(a)] # delete taxonomy column
}

# translate otu number to sequence...
s = readDNAStringSet(rep_seqs)
names(s) <- lapply(names(s), function(seq) {str_split(seq, " ")[[1]][1]}) %>% unlist
rownames(a) <- lapply(rownames(a), function(otu) {s[[otu]] %>% as.character})


# dada2 to fasta (re-replication, necessary for OTU assignment, before rarification)
for (i in seq(ncol(a))) {
  ids <- paste0(colnames(a)[i],
                '_',
                seq(colSums(a)[i]))

  seqs.re_replicated <- rep(rownames(a),times=a[[i]])

  writeFasta(object=ShortRead(sread=DNAStringSet(seqs.re_replicated),
                              id=BStringSet(ids)),
             file=paste0(fasta_temp_folder,"/",
                         str_remove(colnames(a)[i],
                                    ".fastq.gz"),"_dada2cleaned",".fasta"),
             width=width)
}

system(paste0('cat ', fasta_temp_folder, '/* > ', output))
system(paste("rm -rf", fasta_temp_folder))

