
library(dplyr)
library(ggplot2)
library(seqinr)
library(ape)
library(phangorn)
library(msa)
library(ggtree)
library(DECIPHER)

setwd("C:/Users/Alícya/Desktop/Fabricio/Projetos/Phylogenetic_diversity/New/markers/2_COI/local_database/5_blast/blast_thresholds/leray")

blast. <- read.table("blast_result_pid_96.0.txt", sep="", header = F)

colnames(blast.) <- c("qseqid", "seqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Identify unique values in the column
unique_values <- unique(blast.$qseqid)

# Create a sequence of strings for each unique value
replacement_values <- paste0("Otu", seq(unique(blast.$qseqid)))

blast.$qseqid <- ifelse(blast.$qseqid %in% unique_values, replacement_values[match(blast.$qseqid, unique_values)], blast.$qseqid)

#Separate only the max value of match per OTU
blast.max <- blast. %>% group_by(qseqid) %>% top_n(1,pident)

#Copy the "blast.max" object
blast_ <- blast.max

#Organizing the names of species
blast.max$seqid <- gsub("\\d+_", "", blast.max$seqid)
blast.max$seqid <- sapply(blast.max$seqid, function(cell) {
  if (!grepl("cf\\.|aff\\.", cell)) {
    sub("^(.*?_[^_]+)_.*", "\\1", cell)
  } else {
    cell
  }
})

# Read the database file
fishes_database <- read.fasta(file = "C:/Users/Alícya/Desktop/Fabricio/Projetos/Phylogenetic_diversity/New/markers/2_COI/local_database/5_blast//fishes_database.fasta", as.string = TRUE)

selected_rows <- blast_[duplicated(blast_$qseqid) | duplicated(blast_$qseqid, fromLast = TRUE), ]
selected_rows <- blast.[blast.$qseqid %in% selected_rows$qseqid, ]

# Assuming 'dataframe' is your dataframe and 'category' is the column you want to use for grouping
categories <- unique(selected_rows$qseqid)

# Create a list to store the tree objects
result_list <- list()
names_list <- list()


#Loop to select | align | run neighbour joining tree | plot check figures
for (cat in categories) {
  # Subset the dataframe based on the current category
  subset_df <- selected_rows[selected_rows$qseqid == cat, ]
  
  # Select sequences and pi based on the current subset
  seqs <- unlist(subset_df[, 2])
  pi <- unlist(subset_df[, 3])
  
  # Find matching headers
  filtered_sequences <- fishes_database[names(fishes_database) %in% seqs]
  
  # Create a matrix with sequences
  sequences <- sapply(filtered_sequences, function(x) x[[1]])
  names(sequences) <-  paste0(seqs, sep = " ; ", sep = "pi = ", pi)
  
  # Perform multiple sequence alignment
  alignment <- msa(sequences, type = "DNA")
  
  # Compute distance matrix
  dist_matrix <- dist.dna(as.DNAbin(as.matrix(alignment)), gamma = TRUE, model = "JC69")
  
  # Determine the number of values in the current category
  num_values <- nrow(subset_df)
  
  # Choose the tree construction method based on the number of values
  if (num_values < 3) {
    nj_tree <- upgma(dist_matrix)
  } else {
    nj_tree <- nj(dist_matrix)
  }
  
  # Plot the tree
  #plot_tree  <- plot(nj_tree)
  
  # Save the plot in pdf format with a filename based on the category
  plot_filename <- paste0("check_", gsub("/", "_", cat), ".png")
  png(plot_filename, width = 800, height = 600)
  plot(nj_tree)
  dev.off()
  
  # Store the tree object in the result list
  result_list[[cat]] <- nj_tree
  
  # Store the names in the names_list
  names_list[[cat]] <- seqs
  
}

#Select the name of selected sequence based on names_list object

Selected_Sequences <- c("47027_Apistogramma_cacatuoides_voucher_KW11T030_c", #Species doesn't occurs here - We use Apistogramma sp.
                        "21_Corydoras_melanistius_GIBI", # Select base in a pi, because only one of the two options have identification at species level
                        "46930_Hemigrammus_ocellifer_voucher_KW11-T322_cyt", # Select base in a pi, because only one of the two options have identification at species level
                        "49_Hyphessobrycon_rosaceus_GIBI", # One of the two options occurs only in another different basin 
                        "48064_Jupiaba_anterior_voucher_LBPV-57143_cytochr", # Based in Species resolution we have two options (Jupiaba anterior and Astianax anterior - old specie name)
                        "52_Knodus_savannensis_GIBI", # Olny one with specie level identification
                        "56_Metynnis_maculatus_GIBI", # Select based in major occurence (4 - 10) 
                        "75_Serrasalmus_eigenmanni_GIBI", # Select base in a pi, because only one of the two options have identification at species level 
                        "20115_Tatia_punctata_voucher_SU08-1280_cytochrome") # One of the two options occurs only in another different basin

new_correct_rows <- blast.[grepl(paste(Selected_Sequences, collapse="|"), blast.$seqid, ignore.case = TRUE), ] 

#Change the output names
new_correct_rows$seqid <- gsub("\\d+_", "", new_correct_rows$seqid)
new_correct_rows$seqid <- sapply(new_correct_rows$seqid, function(cell) {
  if (!grepl("cf\\.|aff\\.", cell)) {
    sub("^(.*?_[^_]+)_.*", "\\1", cell)
  } else {
    cell
  }
})


#Separate only the max value of match per OTU
blast_max <- blast.max[!duplicated(blast.max$qseqid) & !duplicated(blast.max$qseqid, fromLast = TRUE), ]

blast_final <- data.frame(rbind(blast_max,new_correct_rows))

print("Molecular Identification Finished")