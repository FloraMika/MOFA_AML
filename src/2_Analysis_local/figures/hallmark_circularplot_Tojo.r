# =========================
# Libraries
# =========================
library(tidyverse)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(gtools)
library(stringr)

# =========================
# Load data
# =========================
mfactor_data <- read.xlsx(
  "~/AMLmultiomics/data-drivenFeatures.xlsx",
  sheet = "Blad2"
)

# =========================
# Prepare data
# =========================
mfactor_data$id <- as.numeric(rownames(mfactor_data))
mfactor_data$value <- mfactor_data$value * 100
mfactor_data <- as.data.frame(mfactor_data)

# =========================
# Add empty bars between Hallmarks
# =========================
empty_bar <- 2

to_add <- data.frame(
  matrix(
    NA,
    empty_bar * length(unique(mfactor_data$Hallmark)),
    ncol(mfactor_data)
  )
)

colnames(to_add) <- colnames(mfactor_data)

to_add$Hallmark <- rep(
  unique(mfactor_data$Hallmark),
  each = empty_bar
)

mfactor_data <- rbind(mfactor_data, to_add)

# Sort Hallmarks
mfactor_data <- mfactor_data[gtools::mixedorder(mfactor_data$Hallmark), ]

# Recreate ID
mfactor_data$id <- seq(1, nrow(mfactor_data))

# =========================
# Label positions
# =========================
label_data <- mfactor_data

number_of_bar <- nrow(label_data)

angle <- 90 - 360 * (label_data$id - 0.5) / number_of_bar

label_data$hjust <- ifelse(angle < -90, 1, 0)

label_data$angle <- ifelse(angle < -90, angle + 180, angle)

# =========================
# Clean labels
# =========================
label_data$label_clean <- mfactor_data$label

label_data <- label_data %>%
  mutate(
    label_clean2 = str_replace(
      label_clean,
      "UBAC2;MIR548AN",
      "UBAC2"
    )
  ) %>%
  mutate(
    label_clean2 = str_replace(
      label_clean2,
      "ENSG00000217801.g",
      "RP11-465B22.3"
    )
  )

# =========================
# Text position values
# =========================
label_data <- label_data %>%
  mutate(
    textvalue = if_else(
      is.na(value),
      "NA",
      "100"
    )
  ) %>%
  mutate_at(
    c("textvalue"),
    as.numeric
  )

# =========================
# Base data
# =========================
base_data <- mfactor_data %>%
  group_by(Hallmark) %>%
  summarize(
    start = min(id),
    end = max(id) - empty_bar
  ) %>%
  rowwise() %>%
  mutate(
    title = mean(c(start, end))
  )

base_data <- base_data[
  gtools::mixedorder(base_data$Hallmark),
]

# =========================
# Grid data
# =========================
grid_data <- base_data

grid_data$end <- grid_data$end[
  c(nrow(grid_data), 1:nrow(grid_data)-1)
] + 1

grid_data$start <- grid_data$start - 1

grid_data <- grid_data[-1, ]

# =========================
# Label positions for scale text
# =========================
label_position <- c(
  174.53,
  174.4,
  174.35,
  174.3,
  174.15
)

# =========================
# Final Circular Plot
# =========================
options(repr.plot.width = 39, repr.plot.height = 38)

p <- ggplot(
  mfactor_data,
  aes(
    x = as.factor(id),
    y = value,
    fill = Hallmark
  )
) +

  # Main bars
  geom_bar(
    stat = "identity",
    alpha = 0.5
  ) +

  # Hallmark colors
  scale_fill_manual(
    values = c(
      "Factor1"  = "#6A3D9A",
      "Factor2"  = "#1F78B4",
      "Factor3"  = "#A6CEE3",
      "Factor4"  = "#B2DF8A",
      "Factor5"  = "#33A02C",
      "Factor6"  = "#FFFF99",
      "Factor7"  = "#FDBF6F",
      "Factor8"  = "#FF7F00",
      "Factor9"  = "#FB9A99",
      "Factor10" = "#E31A1C",
      "Factor11" = "#CAB2D6",

      "Factor1N"  = "#FFFFFF",
      "Factor2N"  = "#FFFFFF",
      "Factor3N"  = "#FFFFFF",
      "Factor4N"  = "#FFFFFF",
      "Factor5N"  = "#FFFFFF",
      "Factor6N"  = "#FFFFFF",
      "Factor7N"  = "#FFFFFF",
      "Factor8N"  = "#FFFFFF",
      "Factor9N"  = "#FFFFFF",
      "Factor10N" = "#FFFFFF",
      "Factor11N" = "#FFFFFF"
    )
  ) +

  # Circular grid lines
  geom_segment(
    data = grid_data,
    aes(
      x = end,
      y = 80,
      xend = start,
      yend = 80
    ),
    colour = "black",
    alpha = 1,
    size = 1.8,
    inherit.aes = FALSE
  ) +

  geom_segment(
    data = grid_data,
    aes(
      x = end,
      y = 60,
      xend = start,
      yend = 60
    ),
    colour = "black",
    alpha = 1,
    size = 1.8,
    inherit.aes = FALSE
  ) +

  geom_segment(
    data = grid_data,
    aes(
      x = end,
      y = 40,
      xend = start,
      yend = 40
    ),
    colour = "black",
    alpha = 1,
    size = 1.8,
    inherit.aes = FALSE
  ) +

  geom_segment(
    data = grid_data,
    aes(
      x = end,
      y = 20,
      xend = start,
      yend = 20
    ),
    colour = "black",
    alpha = 1,
    size = 1.8,
    inherit.aes = FALSE
  ) +

  # Radial axis labels
  annotate(
    "text",
    x = label_position,
    y = c(20, 40, 60, 80, 100),
    label = c("0.2", "0.4", "0.6", "0.8", "1.0"),
    color = "black",
    size = 12,
    angle = 0,
    fontface = "bold",
    hjust = 1
  ) +

  ylim(-100, 120) +

  # Theme
  theme_minimal() +

  theme(
    legend.position = c(1.05, 0.3),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 16,
      face = "bold"
    ),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")
  ) +

  coord_polar() +

  # Gene labels
  geom_text(
    data = label_data,
    aes(
      x = id,
      y = textvalue + 15,
      label = label_clean2,
      hjust = hjust
    ),
    color = "black",
    fontface = "bold",
    alpha = 0.6,
    size = 16,
    angle = label_data$angle,
    inherit.aes = FALSE
  ) +

  # Baseline separators
  geom_segment(
    data = base_data,
    aes(
      x = start,
      y = -5,
      xend = end,
      yend = -5
    ),
    color = "black",
    alpha = 0.8,
    size = 1,
    inherit.aes = FALSE
  ) +

  # Sign points
  geom_point(
    data = label_data,
    aes(
      x = id,
      y = textvalue + 5,
      color = sign
    ),
    shape = 19,
    size = 14
  ) +

  # View points
  geom_point(
    data = label_data,
    aes(
      x = id,
      y = textvalue + 11,
      color = view
    ),
    shape = 19,
    size = 14
  ) +

  # Point colors
  scale_color_manual(
    values = c(
      "+" = "black",
      "-" = "grey",

      "DNAm"       = "#A6D854",
      "Mutations"  = "#8DA0CB",
      "Proteomics" = "#66C2A5",
      "RNAseq"     = "#E78AC3",
      "DSRT"       = "#FC8D62"
    )
  )

# Display plot
p

# =========================
# Save plot
# =========================
ggsave(
  plot = p,
  filename = "hallmark_circular_plot.pdf",
  device = "pdf",
  dpi = 600,
  width = 47,
  height = 46,
  bg = "white"
)




