#' @title plotTracks
#'
#' @description Plot tracks for a given genomic region of interest (ROI).
#'
#' @details plotTracks is used to take the ChIP, repeat, and gene annotation data using \code{calcTracks} and plot tracks 
#' for a given region of interest (ROI).
#'
#' @param trackData A list of data frames generated with calcTracks.
#' @param bigwigNames  Character vector containing the names to describe the \code{bigwigFiles} you are using (for example: "H3K9me3_replicate1"). 
#' @param color_map Named vector indicating the plotting color for each group of tracks.(for example: c("H3K9me3" = "red")).
#'
#' @return A plot of genomic tracks of the data in \code{bigwigFiles} and the repeat and gene annotation.
#'
#' @examples
#' #see vignette
#'
#' @importFrom dplyr group_by summarise filter pull
#' @importFrom ggplot2 ggplot geom_area facet_wrap scale_fill_manual scale_x_continuous scale_y_continuous annotate
#' @importFrom ggplot2 labs theme element_blank element_text element_rect margin arrow geom_text aes geom_segment
#' @importFrom cowplot theme_cowplot
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom stats na.omit
#' @importFrom S4Vectors mcols
#'
#' @export
plotTracks <- function(trackData, bigwigNames, color_map){ 
  
  combined_coverage <- trackData$combined_coverage
  reps_df <- trackData$reps_df
  genes_df <- trackData$genes_df
  exons_df <- trackData$exons_df
  
  # Determine x-axis limits (shared across all plots)
  x_limits <- range(combined_coverage$index)

  # Calculate maximum score per group
  group_max_scores <- combined_coverage %>%
  group_by(group) %>%
  summarise(max_score = max(na.omit(score)))

  #empty list for plots
  group_plots <- list()
  
  # Reorder the factor levels of the file_name column
  combined_coverage$file_name <- factor(combined_coverage$file_name, levels = bigwigNames)
  
  # Loop through each group to create individual group plots
  for (grp in unique(combined_coverage$group)) {
    # Filter data for the current group
    group_data <- combined_coverage %>% dplyr::filter(group == grp)
    
    # Get the maximum score for the group
    max_score <- group_max_scores %>% dplyr::filter(group == grp) %>% pull(max_score)
    
    # Initialize an empty vector for grpcolor
    grpcolor <- c()
    
    for (g in grp) {
      grpcolor <- c(grpcolor, rep(color_map[g], length(unique(group_data$file_name))))  
      }
    
  
    # Create a coverage plot for the group
    group_plot <- ggplot(group_data) +
      geom_area(aes(x = index, y = score, fill = group)) +
      facet_wrap(~ file_name, ncol = 1) +  # Facet for individual tracks
      scale_fill_manual(values = grpcolor) +
      scale_x_continuous(limits = x_limits) +
      scale_y_continuous(limits = c(0, max_score)) +  # Scale to group max
      labs(x = "Genomic Position", y = paste0("Cpm in ", grp)) +
      
      
      # Adjust themes
      theme_cowplot() +
      theme(
        panel.grid = element_blank(),  # Remove grid lines
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none",
        panel.spacing = unit(0.5, "lines"),
        axis.text.y = element_text(size = 10),  
        axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),
        axis.text.x = element_blank(),  # Hide x-axis text
        axis.ticks.x = element_blank(),  # Hide x-axis ticks
        axis.title.x =  element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        plot.margin = margin(t = 10, b = 0, l = 10, r = 10)
      ) 
    
    # Add the group's plot to the list
    group_plots[[grp]] <- group_plot
  }
  
  # Create the repeat annotation plot (shared across all groups)
  annotation_plot <- ggplot() +
    geom_segment(data = reps_df, 
                 aes(x = start, xend = end, y = ypos, yend = ypos),
                 # arrow = arrow(length = unit(0.4, "inches"), type = "closed"),
                 arrow = arrow(length = unit(0.1, "inches")),
                 lineend="butt",linejoin="mitre",
                 color = "black", size = 8) +
    geom_text(data = reps_df,
              aes(x = (start + end) / 2, y = ypos, label = repeat_name),
              color = "white", size = 4.5, hjust = 0.5) +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = c(-8,8)) +
    labs(x = unique(combined_coverage$seqnames), y = NULL) +
    
    
    theme_cowplot() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 2, b = 10, l = 10, r = 10)
    )
  
  
  # Create the gene annotation plot (shared across all groups)
  gene_annotation_plot <- ggplot() +
    geom_segment(data = genes_df, 
                 aes(x = start, xend = end, y = ypos, yend = ypos),
                 # arrow = arrow(length = unit(0.4, "inches"), type = "closed"),
                 arrow = arrow(length = unit(0.05, "inches")),
                 lineend="butt",linejoin="mitre",
                 color = "black", size = 1) + 
    geom_segment(data = exons_df, 
                 aes(x = start, xend = end, y = ypos, yend = ypos),
                 # arrow = arrow(length = unit(0.4, "inches"), type = "closed"),
                 arrow = arrow(length = unit(0.1, "inches")),
                 lineend="butt",linejoin="mitre",
                 color = "black", size = 6) +
    geom_text(data = genes_df,
              aes(x = (start + end) / 2, y = ypos + 3, label = gene_name),
              color = "black", size = 4.5, hjust = 0.5) +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = c(-8,9)) +
    labs(x = unique(combined_coverage$seqnames), y = NULL) +
    
    # Add scale bar (200bp) 
    annotate("segment",
             x = max(x_limits) - 200, xend = max(x_limits),
             y = 7, yend =  7,
             size = 0.8, color = "black") +
    annotate("text",
             x = max(x_limits) - 100, y = 6,
             label = "200bp", size = 3, hjust = 0.5) +
    
    theme_cowplot() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 2, b = 10, l = 10, r = 10)
    )
  
  # Combine all group plots and the annotation plot using patchwork
  final_plot <- wrap_plots(group_plots, ncol = 1) / annotation_plot + gene_annotation_plot +
    plot_layout(heights = c(rep(5, length(group_plots)), 2, 2))  # Adjust heights for better spacing
  
  # Display the final plot
  return(final_plot)


}
