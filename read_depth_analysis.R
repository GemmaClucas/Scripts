analyze_read_depths <- function(feature_table_path, 
                                taxonomy_path, 
                                metadata_path, 
                                output_folder,
                                save_output = TRUE) {
  
  # Load required packages
  require(tidyverse)
  require(qiime2R)
  require(knitr)
  
  # Read in the QIIME 2 artifacts and metadata
  feature_table <- read_qza(feature_table_path)
  taxonomy_data <- read_qza(taxonomy_path)
  metadata <- read_q2metadata(metadata_path)
  
  # Convert feature table to data frame
  feature_table_df <- as.data.frame(feature_table$data) %>%
    rownames_to_column("FeatureID")
  
  # Clean up taxonomy data
  taxonomy_clean <- taxonomy_data$data %>%
    parse_taxonomy() %>%
    rownames_to_column("FeatureID")
  
  # Reshape feature table to long format
  feature_long <- feature_table_df %>%
    gather(key = "SampleID", value = "Abundance", -FeatureID)
  
  # Merge all data together
  complete_data <- feature_long %>%
    left_join(taxonomy_clean, by = "FeatureID") %>%
    left_join(metadata, by = "SampleID")
  
  # Calculate average read depths by Type and Plate
  read_depth_comparison <- complete_data %>%
    group_by(SampleID, Type, Plate) %>%
    summarise(TotalReads = sum(Abundance)) %>%
    group_by(Type, Plate) %>%
    summarise(
      MeanReads = round(mean(TotalReads), 2),
      MedianReads = round(median(TotalReads), 2),
      SDReads = round(sd(TotalReads), 2),
      n = n()
    )
  
  # Create empty list to store results for each plate
  plate_results <- list()
  
  # message("Plate counts at this point:")
  # print(table(metadata$Plate, useNA = "ifany"))
  # 
  # message("Rows with missing Plate:")
  # print(metadata[is.na(metadata$Plate), , drop = FALSE])
  # 
  # message("Plate counts in complete_data:")
  # print(table(complete_data$Plate, useNA = "ifany"))
  # 
  # message("Rows in complete_data with missing Plate:")
  # print(
  #   complete_data %>%
  #     filter(is.na(Plate)) %>%
  #     select(SampleID, Type, Plate, everything())
  # )
  
  # Process each plate separately
  for(current_plate in unique(complete_data$Plate)) {
    # Get all expected types from metadata for this plate
    expected_types <- metadata %>%
      filter(Plate == current_plate) %>%
      pull(Type) %>%
      unique()
    
    # Get actual total reads for blanks and mock (this only sums real data)
    ext_blank_total_reads <- complete_data %>%
      filter(Plate == current_plate, Type == "EXTBLANK") %>%
      group_by(SampleID) %>%
      summarise(TotalReads = sum(Abundance)) %>%
      pull(TotalReads) %>%
      sum()
    
    fld_blank_total_reads <- complete_data %>%
      filter(Plate == current_plate, Type == "FLDBLANK") %>%
      group_by(SampleID) %>%
      summarise(TotalReads = sum(Abundance)) %>%
      pull(TotalReads) %>%
      sum()
    
    pcr_blank_total_reads <- complete_data %>%
      filter(Plate == current_plate, Type == "PCRBLANK") %>%
      group_by(SampleID) %>%
      summarise(TotalReads = sum(Abundance)) %>%
      pull(TotalReads) %>%
      sum()
    
    # Get true numbers from metadata
    true_ext_blank_number <- metadata %>% 
      filter(Plate == current_plate, Type == "EXTBLANK") %>% 
      nrow()
    
    true_fld_blank_number <- metadata %>% 
      filter(Plate == current_plate, Type == "FLDBLANK") %>% 
      nrow()
    
    true_pcr_blank_number <- metadata %>% 
      filter(Plate == current_plate, Type == "PCRBLANK") %>% 
      nrow()
    
    # Filter data for current plate
    plate_data <- read_depth_comparison %>% filter(Plate == current_plate)
    
    # Calculate adjusted averages using true totals and true numbers
    adjusted_ext_blank_avg <- ext_blank_total_reads / true_ext_blank_number
    adjusted_fld_blank_avg <- fld_blank_total_reads / true_fld_blank_number
    adjusted_pcr_blank_avg <- pcr_blank_total_reads / true_pcr_blank_number
    sample_avg <- plate_data$MeanReads[plate_data$Type == "SAMPLE"]
    
    # message("Current plate: ", current_plate)
    # message("Expected types: ", paste(expected_types, collapse = ", "))
    # message("Number of expected types: ", length(expected_types))
    
    # Create base comparison dataframe with all expected types
    adjusted_comparison <- data.frame(
      Type = expected_types,
      Plate = current_plate,
      MeanReads = 0,
      MedianReads = 0,
      SDReads = 0,
      n = 0,
      PercentOfSample = 0
    )
    
    # Update with actual data where available
    for(type in unique(plate_data$Type)) {
      idx <- adjusted_comparison$Type == type
      if(any(idx)) {
        adjusted_comparison[idx, "MeanReads"] <- plate_data$MeanReads[plate_data$Type == type]
        adjusted_comparison[idx, "MedianReads"] <- plate_data$MedianReads[plate_data$Type == type]
        adjusted_comparison[idx, "SDReads"] <- plate_data$SDReads[plate_data$Type == type]
        adjusted_comparison[idx, "n"] <- plate_data$n[plate_data$Type == type]
      }
    }
    
    # Update blank rows with adjusted means and true n
    if("EXTBLANK" %in% expected_types) {
      adjusted_comparison$MeanReads[adjusted_comparison$Type == "EXTBLANK"] <- round(adjusted_ext_blank_avg, 2)
      adjusted_comparison$n[adjusted_comparison$Type == "EXTBLANK"] <- true_ext_blank_number
    }
    if("FLDBLANK" %in% expected_types) {
      adjusted_comparison$MeanReads[adjusted_comparison$Type == "FLDBLANK"] <- round(adjusted_fld_blank_avg, 2)
      adjusted_comparison$n[adjusted_comparison$Type == "FLDBLANK"] <- true_fld_blank_number
    }
    if("PCRBLANK" %in% expected_types) {
      adjusted_comparison$MeanReads[adjusted_comparison$Type == "PCRBLANK"] <- round(adjusted_pcr_blank_avg, 2)
      adjusted_comparison$n[adjusted_comparison$Type == "PCRBLANK"] <- true_pcr_blank_number
    }
    
    # Calculate percentages for all types
    adjusted_comparison$PercentOfSample <- round((adjusted_comparison$MeanReads / sample_avg) * 100, 2)
    adjusted_comparison$PercentOfSample[adjusted_comparison$Type == "SAMPLE"] <- 100
    
    # Store results for this plate
    plate_results[[as.character(current_plate)]] <- adjusted_comparison
  }
  
  
  # Calculate combined statistics across all plates
  all_plates_summary <- complete_data %>%
    group_by(SampleID, Type) %>%
    summarise(TotalReads = sum(Abundance)) %>%
    group_by(Type) %>%
    summarise(
      MeanReads = round(mean(TotalReads), 2),
      MedianReads = round(median(TotalReads), 2),
      SDReads = round(sd(TotalReads), 2),
      n = n()
    )
  
  # Get true numbers from metadata across all plates
  true_ext_blank_number <- metadata %>% 
    filter(Type == "EXTBLANK") %>% 
    nrow()
  
  true_fld_blank_number <- metadata %>% 
    filter(Type == "FLDBLANK") %>% 
    nrow()
  
  true_pcr_blank_number <- metadata %>% 
    filter(Type == "PCRBLANK") %>% 
    nrow()
  
  # Calculate total reads for each type across all plates
  ext_blank_total_reads <- complete_data %>%
    filter(Type == "EXTBLANK") %>%
    group_by(SampleID) %>%
    summarise(TotalReads = sum(Abundance)) %>%
    pull(TotalReads) %>%
    sum()
  
  fld_blank_total_reads <- complete_data %>%
    filter(Type == "FLDBLANK") %>%
    group_by(SampleID) %>%
    summarise(TotalReads = sum(Abundance)) %>%
    pull(TotalReads) %>%
    sum()
  
  pcr_blank_total_reads <- complete_data %>%
    filter(Type == "PCRBLANK") %>%
    group_by(SampleID) %>%
    summarise(TotalReads = sum(Abundance)) %>%
    pull(TotalReads) %>%
    sum()
  
  # Calculate adjusted averages
  adjusted_ext_blank_avg <- ext_blank_total_reads / true_ext_blank_number
  adjusted_fld_blank_avg <- fld_blank_total_reads / true_fld_blank_number
  adjusted_pcr_blank_avg <- pcr_blank_total_reads / true_pcr_blank_number
  sample_avg <- all_plates_summary$MeanReads[all_plates_summary$Type == "SAMPLE"]
  
  # Create base comparison dataframe with all types
  all_plates_comparison <- data.frame(
    Type = unique(metadata$Type),
    MeanReads = 0,
    MedianReads = 0,
    SDReads = 0,
    n = 0,
    PercentOfSample = 0
  )
  
  # Update with actual data where available
  for(type in unique(all_plates_summary$Type)) {
    idx <- all_plates_comparison$Type == type
    if(any(idx)) {
      all_plates_comparison[idx, "MeanReads"] <- all_plates_summary$MeanReads[all_plates_summary$Type == type]
      all_plates_comparison[idx, "MedianReads"] <- all_plates_summary$MedianReads[all_plates_summary$Type == type]
      all_plates_comparison[idx, "SDReads"] <- all_plates_summary$SDReads[all_plates_summary$Type == type]
      all_plates_comparison[idx, "n"] <- all_plates_summary$n[all_plates_summary$Type == type]
    }
  }
  
  # Update blank rows with adjusted means and true n
  all_plates_comparison$MeanReads[all_plates_comparison$Type == "EXTBLANK"] <- round(adjusted_ext_blank_avg, 2)
  all_plates_comparison$n[all_plates_comparison$Type == "EXTBLANK"] <- true_ext_blank_number
  
  all_plates_comparison$MeanReads[all_plates_comparison$Type == "FLDBLANK"] <- round(adjusted_fld_blank_avg, 2)
  all_plates_comparison$n[all_plates_comparison$Type == "FLDBLANK"] <- true_fld_blank_number
  
  all_plates_comparison$MeanReads[all_plates_comparison$Type == "PCRBLANK"] <- round(adjusted_pcr_blank_avg, 2)
  all_plates_comparison$n[all_plates_comparison$Type == "PCRBLANK"] <- true_pcr_blank_number
  
  # Calculate percentages for all types
  all_plates_comparison$PercentOfSample <- round((all_plates_comparison$MeanReads / sample_avg) * 100, 2)
  all_plates_comparison$PercentOfSample[all_plates_comparison$Type == "SAMPLE"] <- 100
  
  # # Display plate-specific results first
  # for(plate in names(plate_results)) {
  #   cat(sprintf("\n### Read Depth Summary - Plate %s\n", plate))
  #   print(knitr::kable(plate_results[[plate]], 
  #                      caption = sprintf("Summary of read depths by sample type - Plate %s", plate),
  #                      col.names = c("Type", "Plate", "Mean Reads", "Median Reads", "SD Reads", "n", "% of Sample Reads"),
  #                      align = c('l', 'l', 'r', 'r', 'r', 'r', 'r')))
  #   cat("\n")
  # }
  
  # Display combined results
  cat("\n### Read Depth Summary - All Plates Combined\n")
  print(knitr::kable(all_plates_comparison,
                     caption = "Summary of read depths by sample type - All Plates Combined",
                     col.names = c("Type", "Mean Reads", "Median Reads", "SD Reads", "n", "% of Sample Reads"),
                     align = c('l', 'r', 'r', 'r', 'r', 'r')))
  
  
# Modify the save output section to use the output_folder parameter
  if(save_output) {
    for(plate in names(plate_results)) {
      write.csv(plate_results[[plate]],
                file.path(output_folder, sprintf("adjusted_read_depth_comparison_plate_%s.csv", plate)),
                row.names = FALSE)
    }
    write.csv(all_plates_comparison,
              file.path(output_folder, "adjusted_read_depth_comparison_all_plates.csv"),
              row.names = FALSE)
  }
  
  # Return results as a list
  return(list(
    plate_results = plate_results,
    all_plates_comparison = all_plates_comparison
  ))
}