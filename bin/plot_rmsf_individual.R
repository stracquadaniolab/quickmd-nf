#!/usr/bin/env Rscript
library(docopt)
library(ggplot2)
library(dplyr)
#library(bio3d)
doc <- "Usage:
plot_rmsf_individual.R --incsv=<csv> --wtcsv=<wtcsv> --pH=<pH>
Options:
--incsv=<csv>  load input csv file
--wtcsv=<wtcsv> load input WT csv file
--pH=<pH>  pH represented by protein protonation state
"
opts <- docopt(doc)
csv <- opts$`--incsv`
wtcsv <- opts$`--wtcsv`
pH <- opts$`--pH`
#Read in variant and WT data, define output file and title card
data <- read.csv(csv)
data_wt <- read.csv(wtcsv)
stem <- sub(".csv", "", csv)
stem_wt <- sub(".csv", "", wtcsv)
pngout <- paste0(stem, "_", pH, "pH_RMSF.pdf")
titlecard <- paste(stem, " ", pH, "pH RMSF per Residue")

#assign colours
random_color <- function() { rgb(runif(1), runif(1), runif(1))}

#combine data together to be plotted on the same graphic
combined_data <- rbind(
  data.frame(dataset = stem, data),
  data.frame(dataset = stem_wt, data_wt)
)

# Base plot RMSF per residue
line_plot <- ggplot(combined_data, aes(x = Res_ID, y = RMSF, color = dataset)) +
  geom_line() +
  geom_point() +
  coord_cartesian(xlim=c(0, 398), ylim=c(0, 5)) +
  labs(title = titlecard, x = "Res_ID", y = "RMSF", color = "dataset") +
  theme_minimal()
ggsave(pngout, plot = line_plot, width = 20, height = 4, dpi = 300)
