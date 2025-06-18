# visualize_nextclade_json.R

## Use:
# Rscript viz_nextclade_cli.R \
# --nextclade-input-dir {nextclade_dir} \
# --json-file {nextclade.json} \
# --plotly-output-dir {plotly_dir} \
# --ggplotly-output-dir {ggplotly_dir}

library(argparse)
library(jsonlite)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(plotly)
library(scales)
library(seqinr)

# Custom null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# Create a parser
parser <- ArgumentParser(description="Visualize Nextclade JSON outputs")

# Add arguments
parser$add_argument("--nextclade-input-dir", type="character", required=TRUE,
                    help="Directory containing Nextclade output files (including nextclade.json and cds_translation files)")
parser$add_argument("--json-file", type="character", required=TRUE,
                    help="Path to the main nextclade.json file")
parser$add_argument("--plotly-output-dir", type="character", required=TRUE,
                    help="Output directory for Plotly HTML files")
parser$add_argument("--ggplotly-output-dir", type="character", required=TRUE,
                    help="Output directory for GGPlotly HTML files")

# Parse arguments
args <- parser$parse_args()

# Critical paths:
base_result_dir <- args$nextclade_input_dir
json_file_path <- args$json_file
plotly_output_dir <- args$plotly_output_dir
ggplotly_output_dir <- args$ggplotly_output_dir

# Ensure output directories exist
if(!dir.exists(plotly_output_dir)) dir.create(plotly_output_dir, recursive = TRUE)
if(!dir.exists(ggplotly_output_dir)) dir.create(ggplotly_output_dir, recursive = TRUE)

cat("Nextclade input directory:", base_result_dir, "\n")
cat("Main JSON file:", json_file_path, "\n")
cat("Plotly output directory:", plotly_output_dir, "\n")
cat("GGPlotly output directory:", ggplotly_output_dir, "\n")


cds_translation_files <- list.files(base_result_dir, pattern = "^nextclade\\.cds_translation.*\\.fasta$", full.names = TRUE)
if (length(cds_translation_files) == 0) {
  stop(paste("No 'nextclade.cds_translation.*.fasta' files found in directory:", base_result_dir))
}
genetic_feature_info <- list()
for (fasta_file_path in cds_translation_files) {
  filename <- basename(fasta_file_path)
  feature_name <- sub("^nextclade\\.cds_translation\\.(.*?)\\.fasta$", "\\1", filename)
  tryCatch({
    fasta_content <- read.fasta(fasta_file_path, seqtype = "AA", as.string = TRUE)
    if (length(fasta_content) > 0 && !is.null(fasta_content[[1]][1])) {
      genetic_feature_info[[feature_name]] <- nchar(fasta_content[[1]][1])
    } else {
      genetic_feature_info[[feature_name]] <- NA 
    }
  }, error = function(e) { genetic_feature_info[[feature_name]] <- NA })
}
genetic_feature_info <- genetic_feature_info[!is.na(genetic_feature_info)]
genetic_features_to_plot <- names(genetic_feature_info)
if (length(genetic_features_to_plot) == 0) stop("No genetic features with valid lengths identified.")
message(paste("Identified genetic features to plot:", paste(genetic_features_to_plot, collapse=", ")))

# --- Import main nextclade.json file ---
json_file_path <- file.path(base_result_dir, "nextclade.json")
if (!file.exists(json_file_path)) {
  stop(paste("Error: The main JSON file", json_file_path, "was not found."))
}
# Use simplifyDataFrame=FALSE and simplifyVector=FALSE for robust list structures
raw_json_data <- fromJSON(json_file_path, 
                          flatten = FALSE, 
                          simplifyVector = FALSE, 
                          simplifyDataFrame = FALSE)

if (is.null(raw_json_data$results) || !is.list(raw_json_data$results)) {
  stop("The 'results' field is missing or not a list in the JSON data.")
}

# Extract all sequence names that were successfully processed (are in results)
all_processed_seqnames <- sapply(raw_json_data$results, function(res) {
  if (is.list(res) && "seqName" %in% names(res)) res$seqName else NA_character_
})
all_processed_seqnames <- unique(all_processed_seqnames[!is.na(all_processed_seqnames)])
if (length(all_processed_seqnames) == 0) {
  stop("No sequence names found in JSON results.")
}


# Helper function to create mutation labels (1-indexed for display)
create_mutation_label <- function(type, pos_0idx, ref = NULL, alt = NULL, ins_seq = NULL) {
  pos_1idx <- pos_0idx + 1
  if (type == "Substitution") return(paste0(ref, pos_1idx, alt))
  if (type == "Deletion") return(paste0(ref, pos_1idx, "-"))
  if (type == "Insertion") return(paste0(pos_1idx, ":", ins_seq)) # pos for insertion is AA *before*
  if (type == "Unknown") return("X")
  return(NA_character_)
}


# --- Loop through each genetic feature to create and save plots ---
for (current_feature_name in genetic_features_to_plot) {
  
  message(paste0("\n--- Processing genetic feature: ", current_feature_name, " ---"))
  
  current_feature_ref_length <- genetic_feature_info[[current_feature_name]]
  feature_mutations_list <- list() # To store mutations for THIS feature from all sequences
  
  # Iterate through each sequence in the JSON results
  for (seq_idx in seq_along(raw_json_data$results)) {
    seq_result <- raw_json_data$results[[seq_idx]]
    
    if (!is.list(seq_result) || !"seqName" %in% names(seq_result)) next # Skip malformed
    current_seq_name_from_json <- seq_result[["seqName"]]
    
    # 1. AA Substitutions for current feature
    if (!is.null(seq_result$aaSubstitutions) && length(seq_result$aaSubstitutions) > 0) {
      subs_list <- lapply(seq_result$aaSubstitutions, function(sub) {
        if (is.list(sub) && (sub$cdsName %||% NA_character_) == current_feature_name) {
          tibble(
            seqName = current_seq_name_from_json,
            position = as.numeric(sub$pos %||% NA_real_), # 0-indexed
            mutation_type = "Substitution",
            label = create_mutation_label("Substitution", as.numeric(sub$pos %||% NA_real_), 
                                          as.character(sub$refAa %||% NA_character_), 
                                          as.character(sub$qryAa %||% NA_character_))
          )
        } else { NULL }
      })
      feature_mutations_list <- append(feature_mutations_list, Filter(Negate(is.null), subs_list))
    }
    
    # 2. AA Deletions for current feature
    if (!is.null(seq_result$aaDeletions) && length(seq_result$aaDeletions) > 0) {
      dels_list <- lapply(seq_result$aaDeletions, function(del) {
        if (is.list(del) && (del$cdsName %||% NA_character_) == current_feature_name) {
          tibble(
            seqName = current_seq_name_from_json,
            position = as.numeric(del$pos %||% NA_real_), # 0-indexed
            mutation_type = "Deletion",
            label = create_mutation_label("Deletion", as.numeric(del$pos %||% NA_real_), 
                                          as.character(del$refAa %||% NA_character_))
          )
        } else { NULL }
      })
      feature_mutations_list <- append(feature_mutations_list, Filter(Negate(is.null), dels_list))
    }
    
    # 3. AA Insertions for current feature
    if (!is.null(seq_result$aaInsertions) && length(seq_result$aaInsertions) > 0) {
      ins_list <- lapply(seq_result$aaInsertions, function(ins) {
        # Note: JSON for aaInsertions uses 'cds' not 'cdsName' in your example
        if (is.list(ins) && (ins$cds %||% NA_character_) == current_feature_name) { 
          tibble(
            seqName = current_seq_name_from_json,
            position = as.numeric(ins$pos %||% NA_real_), # 0-indexed AA pos *before* insertion
            mutation_type = "Insertion",
            label = create_mutation_label("Insertion", as.numeric(ins$pos %||% NA_real_), 
                                          ins_seq = as.character(ins$ins %||% NA_character_))
          )
        } else { NULL }
      })
      feature_mutations_list <- append(feature_mutations_list, Filter(Negate(is.null), ins_list))
    }
    
    # 4. Unknown AA Ranges for current feature
    if (!is.null(seq_result$unknownAaRanges) && length(seq_result$unknownAaRanges) > 0) {
      for (unknown_entry in seq_result$unknownAaRanges) { # Each entry is for a CDS
        if (is.list(unknown_entry) && (unknown_entry$cdsName %||% NA_character_) == current_feature_name) {
          if (!is.null(unknown_entry$ranges) && length(unknown_entry$ranges) > 0) {
            
            expanded_unknowns_for_cds <- bind_rows(lapply(unknown_entry$ranges, function(r_spec) {
              # r_spec is like: list(range = list(begin=0, end=227), character="X")
              if(is.list(r_spec) && is.list(r_spec$range) && 
                 !is.null(r_spec$range$begin) && !is.null(r_spec$range$end)) {
                
                range_begin <- as.numeric(r_spec$range$begin %||% NA_real_)
                range_end <- as.numeric(r_spec$range$end %||% NA_real_)
                
                if(!is.na(range_begin) && !is.na(range_end) && range_end > range_begin){
                  return(tibble(position = seq(range_begin, range_end - 1))) # 0-indexed positions
                }
              }
              return(NULL)
            }))
            
            if (nrow(expanded_unknowns_for_cds) > 0) {
              unknowns_to_add <- expanded_unknowns_for_cds %>%
                mutate(
                  seqName = current_seq_name_from_json,
                  mutation_type = "Unknown",
                  label = "X"
                )
              feature_mutations_list <- append(feature_mutations_list, list(unknowns_to_add))
            }
          }
          break # Found the current feature, no need to check other unknownAaRanges for this sequence
        }
      }
    }
  } # End loop over sequences for current feature
  
  if (length(feature_mutations_list) == 0) {
    message(paste("No mutations or unknown AA ranges found for feature:", current_feature_name, ". Skipping plot generation for this feature."))
    next 
  }
  
  combined_mutations_all_types_feature <- bind_rows(feature_mutations_list) %>%
    filter(!is.na(position)) # Ensure position is valid
  
  if (nrow(combined_mutations_all_types_feature) == 0) {
    message(paste("No valid mutations after filtering for feature:", current_feature_name))
    next
  }
  
  # Prepare data for plotting (prioritization, tooltips, etc.)
  combined_mutations_final_feature <- combined_mutations_all_types_feature %>%
    mutate(priority = case_when(mutation_type %in% c("Substitution", "Deletion", "Insertion") ~ 1, 
                                mutation_type == "Unknown" ~ 2, TRUE ~ 3)) %>%
    arrange(seqName, position, priority) %>%
    distinct(seqName, position, .keep_all = TRUE) %>% # Keeps highest priority for overlapping positions
    select(-priority) %>%
    mutate(
      mutation_fill_char = case_when(
        mutation_type == "Unknown" ~ "X",
        mutation_type == "Deletion" ~ "-",
        mutation_type == "Substitution" ~ str_sub(label, -1),
        TRUE ~ NA_character_ # For insertions, which are plotted as points
      )
    )
  
  proximity_threshold <- 3
  aggregated_tooltips_list_feature <- vector("list", length = nrow(combined_mutations_final_feature))
  if (nrow(combined_mutations_final_feature) > 0) {
    combined_mutations_final_feature_sorted <- combined_mutations_final_feature %>% arrange(seqName, position)
    for (k in 1:nrow(combined_mutations_final_feature_sorted)) {
      current_seq_val <- combined_mutations_final_feature_sorted$seqName[k]
      current_pos_val <- combined_mutations_final_feature_sorted$position[k] # 0-indexed
      
      nearby_muts_df <- combined_mutations_final_feature_sorted %>%
        filter(seqName == current_seq_val, 
               abs(position - current_pos_val) <= proximity_threshold) %>%
        arrange(position) %>%
        # Tooltip: "0-indexed_pos mutation_type DisplayLabel(1-indexed)"
        mutate(tooltip_item_str = paste(position, mutation_type, label)) 
      
      aggregated_tooltips_list_feature[[k]] <- paste(nearby_muts_df$tooltip_item_str, collapse = "<br>")
    }
    combined_mutations_final_feature <- combined_mutations_final_feature_sorted %>%
      mutate(aggregated_tooltip = unlist(aggregated_tooltips_list_feature))
  } else {
    combined_mutations_final_feature$aggregated_tooltip <- character(0)
  }
  
  data_for_tiles_feature <- combined_mutations_final_feature %>%
    filter(mutation_type %in% c("Substitution", "Deletion", "Unknown"))
  
  data_for_insertion_points_feature <- combined_mutations_final_feature %>%
    filter(mutation_type == "Insertion")
  
  alphabet_colors <- c(
    "A"="#F0A0FF", "B"="#0075DC", "C"="#993F00", "D"="#4C005C", "E"="#20fbab",
    "F"="#005C31", "G"="#2BCE48", "H"="#FFCC99", "I"="#c75835", "J"="#94FFB5",
    "K"="#8F7C00", "L"="#9DCC00", "M"="#C20088", "N"="#003380", "O"="#FFA405",
    "P"="#FFA8BB", "Q"="#426600", "R"="#FF0010", "S"="#5EF1F2", "T"="#00998F",
    "U"="#E0FF66", "V"="#740AFF", "W"="#990000", "X"="#b2babb", "Y"="#FFE100",
    "Z"="#FF5005", "-"="#000000" 
  )
  insertion_color <- "darkorchid" 
  
  max_x_axis_limit_feature <- if (nrow(combined_mutations_final_feature) > 0 && any(!is.na(combined_mutations_final_feature$position))) {
    max(c(current_feature_ref_length -1, combined_mutations_final_feature$position), na.rm = TRUE) # Max with 0-indexed ref length
  } else {
    current_feature_ref_length -1 # 0-indexed
  }
  max_x_axis_limit_feature <- max(max_x_axis_limit_feature, current_feature_ref_length -1, 0, na.rm = TRUE) # Ensure at least 0 for 0-indexed
  
  # --- ggplotly version for current feature ---
  message(paste("Generating ggplotly version for:", current_feature_name))
  # (Plotting code for ggplotly remains largely the same, just use *_feature variables)
  # Ensure x-axis is 0-indexed, and labels reflect this or convert for display if needed
  p_ggplot_feature <- ggplot() +
    geom_tile(
      data = data_for_tiles_feature,
      aes(x = position, y = seqName, fill = mutation_fill_char, text = aggregated_tooltip),
      width = 5, height = 0.9 # Note: width is in data units (0-indexed positions)
    ) +
    geom_point(
      data = data_for_insertion_points_feature,
      aes(x = position, y = seqName, text = aggregated_tooltip),
      shape = 17, color = insertion_color, size = 2.5
    ) +
    scale_fill_manual(values = alphabet_colors, na.value = "transparent") +
    scale_x_continuous(
      name = paste0("AA Position in ", current_feature_name, " (0-indexed)"),
      limits = c(-0.5, max_x_axis_limit_feature + 10.5), # Adjusted for 0-indexed and tile width
      breaks = scales::pretty_breaks(n = max(5, round((max_x_axis_limit_feature+1)/150))),
      expand = c(0.01, 0.01)
    ) +
    scale_y_discrete(name = NULL, limits = rev(all_processed_seqnames)) +
    labs(title = paste0(current_feature_name, " Amino Acid Alterations (ggplotly)")) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = 'none', axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ggplotly_plot_feature <- ggplotly(p_ggplot_feature, tooltip = "text")
  
  # --- Pure plotly version for current feature ---
  message(paste("Generating pure plotly version for:", current_feature_name))
  # (Plotting code for plotly remains largely the same, use *_feature variables)
  x_breaks_plotly_feature <- scales::pretty_breaks(n = max(5, round((max_x_axis_limit_feature+1)/150)))(c(0, max_x_axis_limit_feature + 10))
  fig_feature <- plot_ly()
  if (nrow(data_for_tiles_feature) > 0) {
    tile_data_for_plotly_feature <- data_for_tiles_feature %>%
      mutate(marker_color = alphabet_colors[mutation_fill_char]) %>%
      filter(!is.na(marker_color)) 
    if(nrow(tile_data_for_plotly_feature) > 0) { 
      fig_feature <- fig_feature %>% add_trace(
        data = tile_data_for_plotly_feature,
        x = ~position, y = ~factor(seqName, levels = all_processed_seqnames), type = 'scatter', mode = 'markers',
        marker = list(symbol = 'line-ns-open', color = ~I(marker_color), size = 15, opacity = 1, line = list(width=0)),
        text = ~aggregated_tooltip, hoverinfo = 'text', showlegend = FALSE, name = paste0(current_feature_name, " Alterations")
      )
    }
  }
  if (nrow(data_for_insertion_points_feature) > 0) {
    fig_feature <- fig_feature %>% add_trace(
      data = data_for_insertion_points_feature,
      x = ~position, y = ~factor(seqName, levels = all_processed_seqnames), type = 'scatter', mode = 'markers',
      marker = list(symbol = 'triangle-up', color = insertion_color, size = 12),
      text = ~aggregated_tooltip, hoverinfo = 'text', showlegend = FALSE, name = paste0(current_feature_name, " Insertions")
    )
  }
  num_sequences_display_feature <- length(all_processed_seqnames)
  plot_dynamic_height_feature <- max(400, 150 + num_sequences_display_feature * 22)
  
  fig_feature <- fig_feature %>% layout(
    height = plot_dynamic_height_feature,
    title = list(text = paste0(current_feature_name, " Amino Acid Alterations (plotly)"), x = 0.5),
    xaxis = list(title = paste0("AA Position in ", current_feature_name, " (0-indexed)"),
                 range = c(-0.5, max_x_axis_limit_feature + 0.5), # Adjusted for 0-indexed
                 tickvals = x_breaks_plotly_feature, ticktext = x_breaks_plotly_feature,
                 tickangle = 45, zeroline = FALSE, showgrid = TRUE, gridcolor = "rgb(229,229,229)"),
    yaxis = list(title = "", categoryorder = "array", categoryarray = rev(all_processed_seqnames),
                 zeroline = FALSE, showgrid = TRUE, gridcolor = "rgb(229,229,229)", dtick = 1),
    legend = list(showlegend = FALSE), plot_bgcolor = "white", paper_bgcolor = "white", hovermode = 'closest', autosize = FALSE
  )
  
  # Export plots
  ggplotly_filename <- file.path(ggplotly_output_dir, paste0(current_feature_name, "_mutations_ggplotly.html"))
  plotly_filename <- file.path(plotly_output_dir, paste0(current_feature_name, "_mutations_plotly.html"))
  
  htmlwidgets::saveWidget(ggplotly_plot_feature, ggplotly_filename, selfcontained = TRUE)
  htmlwidgets::saveWidget(fig_feature, plotly_filename, selfcontained = TRUE)
  message(paste("HTML plots saved for", current_feature_name))
  
  files_dir_to_remove <- sub("\\.html$", "_files", plotly_filename)
  if (dir.exists(files_dir_to_remove)) {
    unlink(files_dir_to_remove, recursive = TRUE, force = TRUE)
    message(paste("Removed directory:", files_dir_to_remove))
  }
  
  files_dir_to_remove <- sub("\\.html$", "_files", ggplotly_filename)
  if (dir.exists(files_dir_to_remove)) {
    unlink(files_dir_to_remove, recursive = TRUE, force = TRUE)
    message(paste("Removed directory:", files_dir_to_remove))
  }
  
} # End of loop for genetic_features

cat("Writing success flag to 'success.txt'", "\n")

writeLines(paste("The rscipt viz_nextclade_cli.R has finished without errors!"), paste0(base_result_dir, "/success.txt"))

cat("R script finished successfully.\n")