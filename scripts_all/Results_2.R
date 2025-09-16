# Load required libraries
library(Biostrings)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Read the multiple alignment FASTA file.
fasta_file <- "/home/kousis/Downloads/output.fasta"
alignment <- readDNAStringSet(fasta_file)

# Extract headers from the FASTA file.
headers <- names(alignment)

# Extract genomic coordinates from headers using regex.
# Assumes headers are in the format: "SampleName.bam|chr02(15499676-15499966)"
coord_pattern <- "\\((\\d+)-(\\d+)\\)"
coords <- str_match(headers, coord_pattern)
start_coord_all <- as.numeric(coords[,2])
end_coord_all   <- as.numeric(coords[,3])

# Create a directory to store the results if it doesn't exist.
output_dir <- "MAF_Results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Identify unique regions based on start and end coordinates.
region_ids <- unique(paste(start_coord_all, end_coord_all, sep = "-"))

# Loop over each unique region.
for (rid in region_ids) {
  # Parse the region coordinates.
  region_parts <- unlist(strsplit(rid, "-"))
  current_start <- as.numeric(region_parts[1])
  current_end   <- as.numeric(region_parts[2])
  
  # Filter indices corresponding to the current region.
  idx_region <- which(start_coord_all == current_start & end_coord_all == current_end)
  
  # Subset the alignment for the current region.
  alignment_region <- alignment[idx_region]
  
  # Extract region information from the first header for the plot title.
  # Example header: "Zoum_E.bam|chr02(15499676-15499966)"
  first_header <- names(alignment_region)[1]
  header_parts <- str_split(first_header, "\\|")[[1]]
  region_info_raw <- header_parts[2]  # e.g., "chr02(15499676-15499966)"
  # Format the region info as "chr02: 15499676-15499966"
  region_title <- sub("(chr\\S+)\\((\\d+-\\d+)\\)", "\\1: \\2", region_info_raw)
  
  # Process sample names: remove everything after ".bam" so that only the sample name remains.
  sample_names <- sapply(names(alignment_region), function(x) sub("\\.bam.*", "", x))
  names(alignment_region) <- sample_names
  
  # Split each sequence into individual characters and create a matrix.
  seq_chars <- lapply(as.character(alignment_region), function(x) unlist(strsplit(x, "")))
  seq_matrix <- do.call(rbind, seq_chars)
  rownames(seq_matrix) <- names(alignment_region)
  
  # Convert the matrix to a data frame.
  df <- as.data.frame(seq_matrix)
  df$Sample <- rownames(df)
  
  # Reshape the data frame from wide to long format.
  df_long <- df %>%
    gather(key = "PosIndex", value = "Base", -Sample) %>%
    mutate(PosIndex = as.numeric(gsub("V", "", PosIndex)),
           GenomePos = current_start + PosIndex - 1)
  
  # Define a custom color mapping for common bases and IUPAC codes.
  base_colors <- c("A" = "#E41A1C", "T" = "#377EB8", "C" = "#4DAF4A", "G" = "#984EA3",
                   "Y" = "#FF7F00", "R" = "#FFFF33", "K" = "#A65628", "M" = "#F781BF",
                   "S" = "#999999", "W" = "#66C2A5", "N" = "gray", "-" = "black")
  
  # Create the tile plot visualization.
  p <- ggplot(df_long, aes(x = GenomePos, y = Sample, fill = Base)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = base_colors, na.value = "gray") +
    labs(title = paste("Multiple Alignment:", region_title),
         x = "Genomic Position",
         y = "Sample") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  # Define the output PDF file name.
  # For example, use the region title with spaces replaced by underscores.
  output_file <- file.path(output_dir, paste0(gsub("[ :]", "_", region_title), ".pdf"))
  
  # Save the plot as a PDF.
  pdf(output_file, width = 25, height = 15)
  print(p)
  dev.off()
  
  message("Saved plot for region ", region_title, " as ", output_file)
}
