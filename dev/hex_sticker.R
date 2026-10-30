
# ============================================================
# Create a hexSticker for the R package "repair"
# Using: hexSticker + ggplot2 + ggimage
# ============================================================

# Install packages if needed
?

# Load packages
library(hexSticker)
library(ggplot2)
library(ggimage)
library(dplyr)
library(tidyr)
library(tibble)


# ============================================================
# 1. CUSTOMIZATION SETTINGS
# ============================================================

# ---- Package title ----
package_name <- "{repair}"

# ---- Image file ----
# Replace this with the path to your image.
# Examples: "my_image.png", "figures/logo.png"
image_file <- "inst/images/handyman.png"


# ---- Colors ----
background <- rgb(50, 50, 50, maxColorValue = 255)     # Hexagon background
border_color <-
  title_color <-
  axis_color <-
  rgb(113, 159, 198, maxColorValue = 255)                 # Hexagon border
line_color <- rgb(208, 110, 104, maxColorValue = 255)     # Spaghetti plot lines


# ---- Fonts ----
title_font <- "mono"         # Try "serif", "sans", or "mono"
title_fontface <- "bold"


# ============================================================
# 2. CREATE EXAMPLE SPAGHETTI DATA
# ============================================================

# This is example data to demonstrate the plot.
# Replace this section with your own data if desired.

set.seed(123)

n_subjects <- 18

plot_data <- expand_grid(
  subject = factor(paste0("Subject_", 1:n_subjects)),
  time = 1:4
) %>%
  mutate(
    # Individual trajectories with different slopes
    value = 2 +
      as.numeric(subject) * 0.10 +
      time * runif(n_subjects, 0.15, 0.95)[as.numeric(subject)] +
      rnorm(n(), 0, 0.12)
  )


# ============================================================
# 3. CREATE THE SPAGHETTI / FANNING SLOPES PLOT
# ============================================================

spaghetti_plot <- ggplot(
  plot_data,
  aes(
    x = time,
    y = value,
    group = subject
  )
) +

  # Individual subject trajectories
  geom_smooth(
    color = line_color,
    linewidth = 0.3,
    alpha = 0.35,
    method = "lm",
    se = FALSE
  ) +

  # Individual observations
  # geom_point(
  #   color = line_color,
  #   size = 1.5,
  #   alpha = 0.8
  # ) +

  # Customize axes if desired
  scale_x_continuous(
    breaks = 1:4,
    labels = c("T1", "T2", "T3", "T4")
  ) +

  # Customize plot labels
  labs(
    x = NULL,
    y = NULL
  ) +

  # Minimal plot appearance
  theme_classic(base_size = 10) +

  theme(
    # Remove background
    panel.background = element_rect(
      fill = NA,
      color = NA
    ),

    plot.background = element_rect(
      fill = NA,
      color = NA
    ),

    # Remove gridlines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    # Remove ticks
    axis.ticks = element_blank(),

    # Axis customization
    axis.text = element_text(
      color = axis_color,
      size = 8
    ),
    axis.line = element_line(
      color = axis_color
    ),

    axis.title = element_blank(),

    # Remove unnecessary margins
    plot.margin = margin(2, 2, 2, 2)
  )


# ============================================================
# 4. CREATE THE HEX STICKER
# ============================================================

# IMPORTANT:
# The s_x, s_y, s_width, and s_height parameters
# control the location and dimensions of the plot.
#
# p_x and p_y control the package title position.

repair_sticker <- sticker(

  # Plot to insert
  subplot = spaghetti_plot,

  # Position of spaghetti plot
  s_x = 0.98,
  s_y = 0.77,

  # Size of spaghetti plot
  s_width = 1.25,
  s_height = 0.95,

  # Package name
  package = package_name,

  # Title position
  p_x = 1,
  p_y = 1.45,

  # Title appearance
  p_color = title_color,
  p_family = title_font,
  p_fontface = title_fontface,
  p_size = 18,

  # Hexagon appearance
  h_fill = background,
  h_color = border_color,
  h_size = 1.5,

  # No spotlight
  spotlight = FALSE,

  # Output filename
  filename = "inst/images/repair_sticker.png",

  # Image resolution
  dpi = 300
)


# ============================================================
# 5. ADD YOUR IMAGE OVER THE PLOT
# ============================================================

# Replace YOUR_IMAGE_FILE.png with your actual image.
#
# This section creates a second sticker layer so the image
# can be positioned independently over the spaghetti plot.
#
# IMPORTANT:
# Adjust image_x, image_y, and image_size to position it.

image_x <- 1.0
image_y <- 0.9
image_size <- 0.35

# Only run this section after supplying a valid image path.

if (file.exists(image_file)) {

  repair_sticker_img <- repair_sticker +

    geom_image(
      aes(
        x = image_x,
        y = image_y,
        image = image_file
      ),
      size = image_size,
      asp = 1
    )

  # Reapply sticker theme and save final version
  repair_sticker_img <- repair_sticker_img +
    theme_sticker(size = 1.5)

  ggsave(
    filename = "inst/images/repair_sticker_final.png",
    plot = repair_sticker_img,
    width = 43.9,
    height = 50.8,
    units = "mm",
    dpi = 300,
    bg = "transparent"
  )

} else {

  message(
    "No image found. Add your image path to image_file ",
    "and rerun the image overlay section."
  )

}


# ============================================================
# 6. PREVIEW
# ============================================================

plot(repair_sticker_img)
