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
library(purrr)
library(slider)

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

# --- Identify genetic features and their lengths (from FASTA files) ---
# (This section is efficient and remains unchanged)
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
# (This section remains unchanged)
if (!file.exists(json_file_path)) {
  stop(paste("Error: The main JSON file", json_file_path, "was not found."))
}
message("Reading main JSON file...")
raw_json_data <- fromJSON(json_file_path, 
                          flatten = FALSE, 
                          simplifyVector = FALSE, 
                          simplifyDataFrame = FALSE)
if (is.null(raw_json_data$results) || !is.list(raw_json_data$results)) {
  stop("The 'results' field is missing or not a list in the JSON data.")
}
message("JSON file loaded.")


# Extract all sequence names that were successfully processed
all_processed_seqnames <- sapply(raw_json_data$results, function(res) res$seqName %||% NA_character_)
all_processed_seqnames <- unique(all_processed_seqnames[!is.na(all_processed_seqnames)])
if (length(all_processed_seqnames) == 0) stop("No sequence names found in JSON results.")


# Helper function to create mutation labels (1-indexed for display)
create_mutation_label <- function(type, pos_0idx, ref = NULL, alt = NULL, ins_seq = NULL) {
  pos_1idx <- pos_0idx + 1
  
  dplyr::case_when(
    type == "Substitution" ~ paste0(ref, pos_1idx, alt),
    type == "Deletion"     ~ paste0(ref, pos_1idx, "-"),
    type == "Insertion"    ~ paste0(pos_1idx, ":", ins_seq),
    type == "Unknown"      ~ "X",
    TRUE                   ~ NA_character_ # Default case for anything else
  )
}


message("Parsing all mutations from JSON into a single dataframe...")

all_mutations_df <- map_dfr(raw_json_data$results, function(seq_result) {
  
  if (!is.list(seq_result) || is.null(seq_result[["seqName"]])) return(NULL)
  current_seq_name <- seq_result[["seqName"]]
  
  # 1. AA Substitutions
  subs_df <- map_dfr(seq_result$aaSubstitutions, ~tibble(
    cdsName = .x$cdsName %||% NA_character_,
    position = as.numeric(.x$pos %||% NA_real_),
    mutation_type = "Substitution",
    ref = as.character(.x$refAa %||% NA_character_),
    alt = as.character(.x$qryAa %||% NA_character_)
  ))
  
  # 2. AA Deletions
  dels_df <- map_dfr(seq_result$aaDeletions, ~tibble(
    cdsName = .x$cdsName %||% NA_character_,
    position = as.numeric(.x$pos %||% NA_real_),
    mutation_type = "Deletion",
    ref = as.character(.x$refAa %||% NA_character_)
  ))
  
  # 3. AA Insertions
  ins_df <- map_dfr(seq_result$aaInsertions, ~tibble(
    cdsName = .x$cds %||% NA_character_, # Note: JSON uses 'cds' here
    position = as.numeric(.x$pos %||% NA_real_),
    mutation_type = "Insertion",
    ins_seq = as.character(.x$ins %||% NA_character_)
  ))
  
  # 4. Unknown AA Ranges
  unknowns_df <- map_dfr(seq_result$unknownAaRanges, function(unknown_entry) {
    feature_name <- unknown_entry$cdsName %||% NA_character_
    map_dfr(unknown_entry$ranges, function(r_spec) {
      range_begin <- as.numeric(r_spec$range$begin %||% NA_real_)
      range_end <- as.numeric(r_spec$range$end %||% NA_real_)
      if (is.na(range_begin) || is.na(range_end) || range_end <= range_begin) return(NULL)
      tibble(
        cdsName = feature_name,
        position = seq(range_begin, range_end - 1), # 0-indexed positions
        mutation_type = "Unknown"
      )
    })
  })
  
  # Combine all mutation types for the current sequence
  bind_rows(subs_df, dels_df, ins_df, unknowns_df) %>%
    mutate(seqName = current_seq_name, .before = 1)
  
}) %>%
  filter(!is.na(position)) # Basic cleaning


# Ensure optional columns exist, even if no such mutations are found in the JSON.
# This prevents errors if a dataset is missing e.g., all insertions.
if (!"ref" %in% names(all_mutations_df)) {
  all_mutations_df <- all_mutations_df %>% mutate(ref = NA_character_)
}
if (!"alt" %in% names(all_mutations_df)) {
  all_mutations_df <- all_mutations_df %>% mutate(alt = NA_character_)
}
if (!"ins_seq" %in% names(all_mutations_df)) {
  all_mutations_df <- all_mutations_df %>% mutate(ins_seq = NA_character_)
}


message("Pre-processing and creating tooltips for all mutations...")

if (nrow(all_mutations_df) > 0) {
  
  proximity_threshold <- 3
  
  all_mutations_processed_df <- all_mutations_df %>%
    mutate(label = create_mutation_label(mutation_type, position, ref, alt, ins_seq)) %>%
    mutate(priority = case_when(mutation_type %in% c("Substitution", "Deletion", "Insertion") ~ 1, 
                                mutation_type == "Unknown" ~ 2, TRUE ~ 3)) %>%
    arrange(seqName, cdsName, position, priority) %>%
    distinct(seqName, cdsName, position, .keep_all = TRUE) %>%
    mutate(
      mutation_fill_char = case_when(
        mutation_type == "Deletion" ~ "-",
        mutation_type == "Substitution" ~ alt,
        TRUE ~ NA_character_
      ),
      tooltip_item_str = paste(position, mutation_type, label)
    ) %>%
    group_by(seqName, cdsName) %>%
    arrange(position, .by_group = TRUE) %>%
    mutate(
      aggregated_tooltip = slider::slide_index_chr(
        .x = tooltip_item_str, .i = position,
        .f = ~paste(.x, collapse = "<br>"),
        .before = proximity_threshold, .after = proximity_threshold
      )
    ) %>%
    ungroup() %>%
    select(-priority, -ref, -alt, -ins_seq, -tooltip_item_str)
  
  
  # Create a separate, summarized dataframe for the Unknown ranges
  message("Consolidating consecutive 'Unknown' ranges for efficient plotting...")
  unknown_ranges_df <- all_mutations_df %>%
    filter(mutation_type == "Unknown") %>%
    group_by(seqName, cdsName) %>%
    arrange(position, .by_group = TRUE) %>%
    # Identify blocks of consecutive positions by checking if the gap is > 1
    mutate(block_id = cumsum(c(1, diff(position)) > 1)) %>%
    group_by(seqName, cdsName, block_id) %>%
    # For each block, get the start and end
    summarise(
      start_pos = min(position),
      end_pos = max(position),
      .groups = 'drop'
    ) %>%
    # Create a simple tooltip for the entire range
    mutate(
      aggregated_tooltip = paste0("Unknown range<br>Positions: ", start_pos, " - ", end_pos)
    ) %>%
    select(-block_id)
  
} else {
  all_mutations_processed_df <- tibble()
  unknown_ranges_df <- tibble() # Ensure it exists but is empty
}
message("Data parsing and pre-processing complete.")

# --- Loop through each genetic feature to CREATE PLOTS from pre-processed data ---
for (current_feature_name in genetic_features_to_plot) {
  
  message(paste0("\n--- Generating plots for genetic feature: ", current_feature_name, " ---"))
  
  feature_df <- all_mutations_processed_df %>% filter(cdsName == current_feature_name)
  
  ## NEW ##
  # Filter the consolidated unknown ranges for the current feature
  data_for_unknown_rects_feature <- unknown_ranges_df %>% filter(cdsName == current_feature_name)
  
  if (nrow(feature_df) == 0 && nrow(data_for_unknown_rects_feature) == 0) {
    message(paste("No alterations found for feature:", current_feature_name, ". Skipping."))
    next 
  }
  
  # Tiles are now only for substitutions and deletions
  data_for_tiles_feature <- feature_df %>%
    filter(mutation_type %in% c("Substitution", "Deletion"))
  
  data_for_insertion_points_feature <- feature_df %>%
    filter(mutation_type == "Insertion")
  
  # --- Plotting Configuration ---
  alphabet_colors <- c(
    "A"="#F0A0FF", "B"="#0075DC", "C"="#993F00", "D"="#4C005C", "E"="#20fbab",
    "F"="#005C31", "G"="#2BCE48", "H"="#FFCC99", "I"="#c75835", "J"="#94FFB5",
    "K"="#8F7C00", "L"="#9DCC00", "M"="#C20088", "N"="#003380", "O"="#FFA405",
    "P"="#FFA8BB", "Q"="#426600", "R"="#FF0010", "S"="#5EF1F2", "T"="#00998F",
    "U"="#E0FF66", "V"="#740AFF", "W"="#990000", "X"="#b2babb", "Y"="#FFE100",
    "Z"="#FF5005", "-"="#000000" 
  )
  insertion_color <- "darkorchid"
  unknown_color <- alphabet_colors[["X"]] # Use the color defined for 'X'
  
  current_feature_ref_length <- genetic_feature_info[[current_feature_name]]
  max_pos_from_muts <- max(c(feature_df$position, data_for_unknown_rects_feature$end_pos), na.rm = TRUE, -Inf)
  max_x_axis_limit_feature <- max(current_feature_ref_length - 1, max_pos_from_muts, na.rm = TRUE)
  max_x_axis_limit_feature <- max(max_x_axis_limit_feature, 0, na.rm = TRUE)
  
  # --- ggplotly version ---
  
  # Initialize the base plot with scales and themes that apply to everything.
  p_ggplot_feature <- ggplot() +
    # Use the complete list of sequences for the Y-axis limits
    scale_y_discrete(name = NULL, limits = rev(all_processed_seqnames)) +
    scale_x_continuous(
      name = paste0("AA Position in ", current_feature_name, " (0-indexed)"),
      limits = c(-0.5, max_x_axis_limit_feature + 10.5),
      breaks = scales::pretty_breaks(n = max(5, round((max_x_axis_limit_feature+1)/150))),
      expand = c(0.01, 0.01)
    ) +
    # Add scales for aesthetics that will be used by the layers
    scale_color_manual(values = alphabet_colors, na.value = "transparent", drop = FALSE) +
    labs(title = paste0(current_feature_name, " Amino Acid Alterations (ggplotly)")) +
    theme_bw(base_size = 10) +
    theme(legend.position = 'none', axis.title.y = element_blank(), 
          plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Add the 'Unknown' ranges layer ONLY if there is data for it
  if (nrow(data_for_unknown_rects_feature) > 0) {
    p_ggplot_feature <- p_ggplot_feature +
      geom_segment(
        data = data_for_unknown_rects_feature,
        aes(x = start_pos, xend = end_pos, y = seqName, yend = seqName, text = aggregated_tooltip),
        color = unknown_color, 
        linewidth = 6 
      )
  }
  
  # Add the substitutions/deletions layer ONLY if there is data for it
  if (nrow(data_for_tiles_feature) > 0) {
    p_ggplot_feature <- p_ggplot_feature +
      geom_point(
        data = data_for_tiles_feature,
        aes(x = position, y = seqName, color = mutation_fill_char, text = aggregated_tooltip),
        shape = 15, 
        size = 2
      )
  }
  
  # Add the insertions layer ONLY if there is data for it
  if (nrow(data_for_insertion_points_feature) > 0) {
    p_ggplot_feature <- p_ggplot_feature +
      geom_point(
        data = data_for_insertion_points_feature,
        aes(x = position, y = seqName, text = aggregated_tooltip),
        shape = 17, 
        color = insertion_color, 
        size = 2.5
      )
  }
  
  # Now, convert the safely constructed ggplot object to a plotly object
  ggplotly_plot_feature <- ggplotly(p_ggplot_feature, tooltip = "text")
  
  # --- Pure plotly version ---
  message("Generating pure plotly version...")
  x_breaks_plotly_feature <- scales::pretty_breaks(n = max(5, round((max_x_axis_limit_feature+1)/150)))(c(0, max_x_axis_limit_feature + 10))
  fig_feature <- plot_ly()
  
  # Use add_segments for the consolidated "Unknown" ranges. This is very efficient.
  if (nrow(data_for_unknown_rects_feature) > 0) {
    fig_feature <- fig_feature %>% add_segments(
      data = data_for_unknown_rects_feature,
      x = ~start_pos, xend = ~end_pos,
      y = ~factor(seqName, levels = all_processed_seqnames), yend = ~factor(seqName, levels = all_processed_seqnames),
      line = list(color = unknown_color, width = 15), # Draw a thick line
      hoverinfo = 'text', text = ~aggregated_tooltip,
      showlegend = FALSE, name = "Unknown Range"
    )
  }
  
  if (nrow(data_for_tiles_feature) > 0) {
    tile_data_for_plotly_feature <- data_for_tiles_feature %>%
      mutate(marker_color = alphabet_colors[mutation_fill_char]) %>%
      filter(!is.na(marker_color))
    fig_feature <- fig_feature %>% add_trace(
      data = tile_data_for_plotly_feature,
      x = ~position, y = ~factor(seqName, levels = all_processed_seqnames), type = 'scatter', mode = 'markers',
      marker = list(symbol = 'line-ns-open', color = ~I(marker_color), size = 15, opacity = 1, line = list(width=0)),
      text = ~aggregated_tooltip, hoverinfo = 'text', showlegend = FALSE, name = "Alterations"
    )
  }
  if (nrow(data_for_insertion_points_feature) > 0) {
    fig_feature <- fig_feature %>% add_trace(
      data = data_for_insertion_points_feature,
      x = ~position, y = ~factor(seqName, levels = all_processed_seqnames), type = 'scatter', mode = 'markers',
      marker = list(symbol = 'triangle-up', color = insertion_color, size = 12),
      text = ~aggregated_tooltip, hoverinfo = 'text', showlegend = FALSE, name = "Insertions"
    )
  }
  
  plot_dynamic_height_feature <- max(400, 150 + length(all_processed_seqnames) * 22)
  fig_feature <- fig_feature %>% layout(
    height = plot_dynamic_height_feature,
    title = list(text = paste0(current_feature_name, " Amino Acid Alterations (plotly)"), x = 0.5),
    xaxis = list(title = paste0("AA Position in ", current_feature_name, " (0-indexed)"),
                 range = c(-0.5, max_x_axis_limit_feature + 0.5),
                 tickvals = x_breaks_plotly_feature, ticktext = x_breaks_plotly_feature,
                 tickangle = 45, zeroline = FALSE, showgrid = TRUE, gridcolor = "rgb(229,229,229)"),
    yaxis = list(title = "", categoryorder = "array", categoryarray = rev(all_processed_seqnames),
                 zeroline = FALSE, showgrid = TRUE, gridcolor = "rgb(229,229,229)", dtick = 1),
    legend = list(showlegend = FALSE), plot_bgcolor = "white", paper_bgcolor = "white", hovermode = 'closest'
  )
  
  # --- Export plots ---
  ggplotly_filename <- file.path(ggplotly_output_dir, paste0(current_feature_name, "_mutations_ggplotly.html"))
  plotly_filename <- file.path(plotly_output_dir, paste0(current_feature_name, "_mutations_plotly.html"))
  
  htmlwidgets::saveWidget(ggplotly_plot_feature, ggplotly_filename, selfcontained = TRUE)
  htmlwidgets::saveWidget(fig_feature, plotly_filename, selfcontained = TRUE)
  message(paste("HTML plots saved for", current_feature_name))
  
  unlink(sub("\\.html$", "_files", plotly_filename), recursive = TRUE, force = TRUE)
  unlink(sub("\\.html$", "_files", ggplotly_filename), recursive = TRUE, force = TRUE)
}

cat("Writing success flag to 'success.txt'", "\n")
writeLines(paste("The rscipt viz_nextclade_cli.R has finished without errors!"), paste0(base_result_dir, "/success.txt"))

message("R script has finished.\n")