# Load required libraries
suppressPackageStartupMessages({
  library(qiime2R)
  library(tidyverse)
  library(biomformat)
})

# Main function
process_data <- function(feature_table_path, taxonomy_path, metadata_path, 
                         output_dir, abundance_threshold) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Read in the QIIME 2 artifacts and metadata
  feature_table <- read_qza(feature_table_path)
  taxonomy_data <- read_qza(taxonomy_path)
  metadata <- read_q2metadata(metadata_path)
  
  # Convert feature table to data frame
  feature_table_df <- as.data.frame(feature_table$data) %>%
    rownames_to_column("FeatureID")
  
  # Extract taxonomy information
  taxonomy_df <- data.frame(
    FeatureID = taxonomy_data$data[,1],
    Taxonomy = taxonomy_data$data[,2],
    stringsAsFactors = FALSE
  )
  
  # Reshape feature table to long format
  feature_long <- feature_table_df %>%
    gather(key = "SampleID", value = "Abundance", -FeatureID) %>%
    # Join with taxonomy
    left_join(taxonomy_df, by = "FeatureID")
  
  # Aggregate abundances by taxonomy and sample
  species_level_abundances <- feature_long %>%
    group_by(SampleID, Taxonomy) %>%
    summarise(Abundance = sum(Abundance)) %>%
    ungroup()
  
  # Apply abundance threshold filter per sample at species level
  filtered_species <- species_level_abundances %>%
    group_by(SampleID) %>%
    mutate(
      TotalReads = sum(Abundance),
      RelativeAbundance = Abundance / TotalReads,
      # Set abundance to 0 if relative abundance is < threshold IN THIS SAMPLE
      Abundance = if_else(RelativeAbundance < abundance_threshold, 0, Abundance)
    ) %>%
    ungroup()
  
  # Get feature-level abundances with the sample-specific filtering applied
  filtered_features <- feature_long %>%
    # Join with the filtered species data to get the filtered abundances
    left_join(
      filtered_species %>% 
        select(SampleID, Taxonomy, Abundance) %>%
        rename(FilteredAbundance = Abundance),
      by = c("SampleID", "Taxonomy")
    ) %>%
    # If the species was filtered out in this sample (FilteredAbundance = 0),
    # set the feature abundance to 0
    mutate(
      Abundance = if_else(FilteredAbundance == 0, 0, Abundance)
    ) %>%
    select(FeatureID, SampleID, Abundance, Taxonomy)
  
  # Remove any features that now have zero abundance across all samples
  filtered_features <- filtered_features %>%
    group_by(FeatureID) %>%
    filter(sum(Abundance) > 0) %>%
    ungroup()
  
  # Convert back to wide format for QIIME2
  filtered_feature_table <- filtered_features %>%
    select(FeatureID, SampleID, Abundance) %>%
    pivot_wider(
      names_from = SampleID,
      values_from = Abundance,
      values_fill = 0
    )
  
  # Create feature table in BIOM format
  filtered_feature_matrix <- as.matrix(filtered_feature_table[,-1])
  rownames(filtered_feature_matrix) <- filtered_feature_table$FeatureID
  
  # Create BIOM object
  biom_obj <- make_biom(data = filtered_feature_matrix)
  
  # Write BIOM file
  write_biom(biom_obj, file.path(output_dir, "filtered_feature_table.biom"))
  
  # Write taxonomy file
  filtered_taxonomy <- taxonomy_df %>%
    filter(FeatureID %in% filtered_feature_table$FeatureID)
  
  filtered_taxonomy_output <- data.frame(
    FeatureID = filtered_taxonomy$FeatureID,
    Taxonomy = filtered_taxonomy$Taxonomy
  )
  
  write.table(filtered_taxonomy_output, 
              file.path(output_dir, "filtered_taxonomy.tsv"), 
              sep="\t", 
              quote=FALSE, 
              row.names=FALSE,
              col.names = c("Feature ID", "Taxon"))
  
  # Print summary statistics
  cat("\nSummary Statistics:\n")
  cat("Original number of features:", nrow(feature_table_df), "\n")
  cat("Number of features after filtering:", nrow(filtered_feature_table), "\n")
  cat("Number of features removed:", nrow(feature_table_df) - nrow(filtered_feature_table), "\n")
  
  # Calculate and print the number of unique taxa before and after filtering
  n_taxa_before <- length(unique(taxonomy_df$Taxonomy))
  n_taxa_after <- length(unique(filtered_features$Taxonomy))
  cat("\nNumber of unique taxa before filtering:", n_taxa_before, "\n")
  cat("Number of unique taxa after filtering:", n_taxa_after, "\n")
  cat("Number of taxa removed:", n_taxa_before - n_taxa_after, "\n")
  
  # Calculate and print the number of reads before and after filtering
  total_reads_before <- sum(feature_long$Abundance)
  total_reads_after <- sum(filtered_features$Abundance)
  cat("\nTotal reads before filtering:", total_reads_before, "\n")
  cat("Total reads after filtering:", total_reads_after, "\n")
  cat("Percentage of reads retained:", round(total_reads_after/total_reads_before * 100, 2), "%\n")
}

