#' Plot EcoTyper Abundance
#'
#' @description
#' Reads EcoTyper output directories and generates highly customizable bar plots
#' representing cell type abundances. It supports generating a combined relative
#' abundance plot or individual absolute abundance plots per cell type.
#'
#' @param folder_path A string. Path to the root directory containing the cell type
#'   subfolders (e.g., output from EcoTyper). Default is \code{"files"}.
#' @param plot_type A string. Either \code{"combined"} to generate a single relative
#'   abundance plot for all cell types, or \code{"per_cell_type"} to generate an
#'   absolute abundance plot for each cell type separately.
#' @param filter_assigned A logical. If \code{TRUE}, the data is subsetted to include
#'   only samples present in the \code{state_assignment.txt} file. Default is \code{TRUE}.
#' @param sample_mapping A named character vector to map substrings in sample names to
#'   presentable group labels for faceting. E.g., \code{c("AMR" = "American")}. Default is \code{NULL}.
#' @param cell_types A character vector specifying which cell types to include in the plot.
#'   If \code{NULL} (default), all available cell types in the directory are used.
#' @param output_dir A string. Directory path where the generated PNG plots will be saved.
#'   Default is the current working directory (\code{"."}).
#'
#' @return Invisible \code{NULL}. The function is called primarily for its side effect
#'   of saving plot images to the disk.
#'
#' @import ggplot2
#' @importFrom magrittr `%>%`
#' @importFrom utils read.delim
#' @importFrom rlang .data
#'
#' @export
plot_abundance <- function(folder_path = "files",
                           plot_type = c("combined", "per_cell_type"),
                           filter_assigned = TRUE,
                           sample_mapping = NULL,
                           cell_types = NULL,
                           output_dir = ".") {

  # 1. Input Validation
  plot_type <- match.arg(plot_type)

  if (!dir.exists(folder_path)) {
    stop("The specified folder_path '", folder_path, "' does not exist.", call. = FALSE)
  }

  if (!is.null(sample_mapping) && is.null(names(sample_mapping))) {
    stop("The 'sample_mapping' parameter must be a named vector. Example: c('Pattern' = 'Label')", call. = FALSE)
  }

  # 2. Setup Directories and Subfolders
  all_dirs <- list.dirs(path = folder_path, full.names = FALSE, recursive = FALSE)

  # Remove root directory representation and strictly ignore the 'Ecotypes' folder
  subfolders <- all_dirs[all_dirs != "" & all_dirs != "Ecotypes"]

  if (length(subfolders) == 0) {
    stop("No cell type subfolders found inside '", folder_path, "'.", call. = FALSE)
  }

  # 3. Filter specific cell types if requested by the user
  if (!is.null(cell_types)) {
    missing_cells <- setdiff(cell_types, subfolders)
    if (length(missing_cells) > 0) {
      warning("The following requested cell types were not found and will be ignored: ",
              paste(missing_cells, collapse = ", "), call. = FALSE)
    }

    subfolders <- intersect(subfolders, cell_types)

    if (length(subfolders) == 0) {
      stop("None of the requested 'cell_types' were found in the directory.", call. = FALSE)
    }
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  # 4. Helper Function to Process Each Cell Type
  process_cell_type <- function(cell_t) {
    abundance_file <- file.path(folder_path, cell_t, "state_abundances.txt")
    assign_file <- file.path(folder_path, cell_t, "state_assignment.txt")

    if (!file.exists(abundance_file)) {
      warning("Abundance file not found for: ", cell_t, ". Skipping...", call. = FALSE)
      return(NULL)
    }

    data <- utils::read.delim(abundance_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # Filter by assigned samples if requested
    if (filter_assigned && file.exists(assign_file)) {
      data_assign <- utils::read.delim(assign_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      valid_samples <- intersect(as.character(data_assign$ID), colnames(data))
      data <- data[, valid_samples, drop = FALSE]
    }

    if (ncol(data) == 0) return(NULL)

    # Create normalized long-format dataframe
    df <- data.frame(
      samples = colnames(data),
      abundance = colSums(data),
      cell_type = cell_t,
      stringsAsFactors = FALSE
    )

    return(df)
  }

  # 5. Aggregate Data
  df_list <- lapply(subfolders, process_cell_type)
  combined_df <- do.call(rbind, df_list)
  rownames(combined_df) <- NULL

  if (is.null(combined_df) || nrow(combined_df) == 0) {
    stop("No valid data could be compiled from the provided directory.", call. = FALSE)
  }

  # 6. Dynamic Group Annotation based on string matching
  if (!is.null(sample_mapping)) {
    combined_df$group <- "Other"

    for (pattern in names(sample_mapping)) {
      matches <- grepl(pattern, combined_df$samples)
      combined_df$group[matches] <- sample_mapping[[pattern]]
    }

    # Convert to factor to preserve the order provided by the user in the mapping
    group_levels <- c(unname(sample_mapping), "Other")
    combined_df$group <- factor(combined_df$group, levels = unique(group_levels))
  }

  # Disable scientific notation for axes
  old_scipen <- getOption("scipen")
  options(scipen = 999)
  on.exit(options(scipen = old_scipen)) # Ensures it resets when function finishes

  # 7. Plot Generation
  if (plot_type == "combined") {

    p <- ggplot2::ggplot(combined_df, aes(fill = .data$cell_type, y = .data$abundance, x = .data$samples)) +
      geom_bar(position = "fill", stat = "identity") +
      labs(
        title = "Cell type abundance per sample",
        subtitle = "EcoTyper Output",
        x = "Sample",
        y = "Abundance (%)",
        fill = "Cell Type"
      ) +
      scale_fill_brewer(palette = "Paired") +
      theme(axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5, hjust = 1, face = "bold"))

    if (!is.null(sample_mapping)) {
      p <- p + facet_grid(~ .data$group, scales = "free", space = "free_x")
    }

    out_file <- file.path(output_dir, "Cell_type_abundance_combined.png")
    ggplot2::ggsave(out_file, plot = p, width = 2560, height = 1440, units = "px")
    message("Successfully saved combined plot: ", out_file)

  } else if (plot_type == "per_cell_type") {

    unique_cells <- unique(combined_df$cell_type)

    for (c_t in unique_cells) {
      df_sub <- combined_df[combined_df$cell_type == c_t, ]
      num_col <- length(unique(df_sub$samples))

      p <- ggplot2::ggplot(df_sub, aes(fill = .data$cell_type, y = .data$abundance, x = .data$samples)) +
        geom_bar(position = "stack", stat = "identity") +
        labs(
          title = paste0("Cell type abundance per sample - ", c_t),
          subtitle = "EcoTyper Output",
          x = paste0("Sample (n = ", num_col, ")"),
          y = "Abundance",
          fill = "Cell Type"
        ) +
        scale_fill_brewer(palette = "Paired") +
        theme(axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5, hjust = 1, face = "bold"))

      if (!is.null(sample_mapping)) {
        p <- p + facet_grid(~ .data$group, scales = "free", space = "free_x")
      }

      out_file <- file.path(output_dir, paste0("Cell_type_abundance_", c_t, ".png"))
      ggplot2::ggsave(out_file, plot = p, width = 2560, height = 1440, units = "px")
    }
    message("Successfully saved ", length(unique_cells), " individual plots in: ", output_dir)
  }

  return(invisible(NULL))
}
