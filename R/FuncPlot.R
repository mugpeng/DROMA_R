create_plot_with_common_axes <- function(p, x_title = "Common X-Axis Title",
                                         y_title = "Common Y-Axis Title") {
  # Create a function that will generate the plot when called
  function() {
    # Convert patchwork to a grob
    p_grob <- patchworkGrob(p)

    # Create a new plotting area
    grid.newpage()

    # Draw the patchwork
    grid.draw(p_grob)

    # Add common x-axis title
    grid.text(x_title,
              x = 0.5, y = 0.02,
              gp = gpar(fontsize = 18, fontface = "bold"))

    # Add common y-axis title (rotated)
    grid.text(y_title,
              x = 0.01, y = 0.5,
              rot = 90,
              gp = gpar(fontsize = 18, fontface = "bold"))

    # Return the grob for potential further use
    invisible(p_grob)
  }
}
