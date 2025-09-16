# Load required libraries
library(Biostrings)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Read the multiple alignment FASTA file.
fasta_file <- "/home/kousis/Downloads/output_2.fasta"
alignment <- readDNAStringSet(fasta_file)
# Extract headers
headers <- names(alignment)
headers
# Extract genomic coordinates from headers using regex.
# Assumes headers are in the format: "SampleName.bam|chr02(15499676-15499966)"
coord_pattern <- "\\((\\d+)-(\\d+)\\)"
coords <- str_match(headers, coord_pattern)
start_coord_all <- as.numeric(coords[,2])
end_coord_all   <- as.numeric(coords[,3])

# Create a directory to store the results if it doesn't exist.
output_dir <- "MAF_Results_1"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

## Part 1: Save alignment plots for each region (as done previously)
# Identify unique regions based on start and end coordinates.
region_ids <- unique(paste(start_coord_all, end_coord_all, sep = "-"))

for (rid in region_ids) {
  # Parse region coordinates.
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
  # Format as "chr02: 15499676-15499966"
  region_title <- sub("(chr\\S+)\\((\\d+-\\d+)\\)", "\\1: \\2", region_info_raw)
  
  # Process sample names to keep only the part before ".bam"
  sample_names <- sapply(names(alignment_region), function(x) sub("\\.bam.*", "", x))
  names(alignment_region) <- sample_names
  
  # Split each sequence into individual characters and create a matrix.
  seq_chars <- lapply(as.character(alignment_region), function(x) unlist(strsplit(x, "")))
  seq_matrix <- do.call(rbind, seq_chars)
  rownames(seq_matrix) <- names(alignment_region)
  
  # Convert matrix to data frame.
  df <- as.data.frame(seq_matrix)
  df$Sample <- rownames(df)
  
  # Reshape the data frame from wide to long format.
  df_long <- df %>%
    gather(key = "PosIndex", value = "Base", -Sample) %>%
    mutate(PosIndex = as.numeric(gsub("V", "", PosIndex)),
           GenomePos = current_start + PosIndex - 1)
  
  # Define custom color mapping for bases.
  base_colors <- c("A" = "#E41A1C", "T" = "#377EB8", "C" = "#4DAF4A", "G" = "#984EA3",
                   "Y" = "#FF7F00", "R" = "#FFFF33", "K" = "#A65628", "M" = "#F781BF",
                   "S" = "#999999", "W" = "#66C2A5", "N" = "gray", "-" = "black")
  
  # Create the tile plot.
  p <- ggplot(df_long, aes(x = GenomePos, y = Sample, fill = Base)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = base_colors, na.value = "gray") +
    labs(title = paste("Multiple Alignment:", region_title),
         x = "Genomic Position",
         y = "Sample") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  # Define output file name.
  output_file <- file.path(output_dir, paste0(gsub("[ :]", "_", region_title), ".pdf"))
  
  # Save the plot as a PDF.
  pdf(output_file, width = 25, height = 15)
  print(p)
  dev.off()
  
  message("Saved alignment plot for region ", region_title, " as ", output_file)
}


## Part 2: Create the diversity heatmap across regions

# This data frame will hold diversity info for every region.
diversity_df <- data.frame()

# Loop over each unique region.
for (rid in region_ids) {
  # Parse region coordinates.
  region_parts <- unlist(strsplit(rid, "-"))
  current_start <- as.numeric(region_parts[1])
  current_end   <- as.numeric(region_parts[2])
  
  # Filter alignment for the current region.
  idx_region <- which(start_coord_all == current_start & end_coord_all == current_end)
  alignment_region <- alignment[idx_region]
  
  # Extract region title from the first header.
  first_header <- names(alignment_region)[1]
  header_parts <- str_split(first_header, "\\|")[[1]]
  region_info_raw <- header_parts[2]
  region_title <- sub("(chr\\S+)\\((\\d+-\\d+)\\)", "\\1: \\2", region_info_raw)
  
  # Process sample names.
  sample_names <- sapply(names(alignment_region), function(x) sub("\\.bam.*", "", x))
  names(alignment_region) <- sample_names
  
  # Split each sequence into individual characters and create a matrix.
  seq_chars <- lapply(as.character(alignment_region), function(x) unlist(strsplit(x, "")))
  seq_matrix <- do.call(rbind, seq_chars)
  
  # Determine region length (number of columns).
  region_length <- ncol(seq_matrix)
  
  # For each position, calculate the diversity (number of unique letters)
  # and record the consensus letter if diversity == 1.
  diversity <- apply(seq_matrix, 2, function(col) length(unique(col)))
  consensus <- apply(seq_matrix, 2, function(col) {
    u <- unique(col)
    if(length(u) == 1) u else NA
  })
  
  # Build a data frame for this region.
  df_region <- data.frame(Region = region_title,
                          Position = 1:region_length,
                          Diversity = diversity,
                          Consensus = consensus,
                          stringsAsFactors = FALSE)
  
  diversity_df <- rbind(diversity_df, df_region)
}

# Order regions as they appeared.
diversity_df$Region <- factor(diversity_df$Region, levels = unique(diversity_df$Region))

# Separate positions with uniform bases (diversity == 1) and diverse positions.
df_uniform <- diversity_df %>% filter(Diversity == 1)
df_diverse <- diversity_df %>% filter(Diversity > 1)

# Create the heatmap plot.
# For positions with diversity > 1, we map a gradient (e.g., low diversity in light blue, high in dark blue).
# For uniform positions (diversity == 1), we fill white and overlay the consensus letter.
p_heatmap <- ggplot() +
  geom_tile(data = df_diverse, aes(x = Position, y = Region, fill = Diversity), color = "black") +
  geom_tile(data = df_uniform, aes(x = Position, y = Region), fill = "white", color = "black") +
  geom_text(data = df_uniform, aes(x = Position, y = Region, label = Consensus), size = 3) +
  scale_fill_gradient(low = "#cde0c9", high = "darkgreen", name = "Diversity\n(count)") +
  labs(title = "Alignment Diversity Heatmap",
       x = "Position in Region",
       y = "Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Save the heatmap as a PDF.
heatmap_file <- file.path(output_dir, "diversity_heatmap.pdf")
pdf(heatmap_file, width = 30, height = 20)
print(p_heatmap)
dev.off()

message("Saved diversity heatmap as ", heatmap_file)
