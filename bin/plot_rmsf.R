#!/usr/bin/env Rscript
library(docopt)
library(ggplot2)
library(dplyr)
#library(bio3d)
doc <- "Usage:
plot_rmsf.R --incsv=<csv>
Options:
--incsv=<csv>  load input csv file
"
opts <- docopt(doc)
csv <- as.numeric(opts$`--csv`)
data <- list()
for (x in csv){
  name <- gsub(".csv", "", x)
  data_item <- read.csv(x)
  data_list <- list(name = data_item)
  data <- c(data, data_list)
}
#data <- lapply(csv, read.csv)
plot_data <- bind_rows(data, .id = "Dataset")
#data <- read.csv(csv)
pngout <- "Comparison_allRMSFs.png"

line_plot <- ggplot(plot_data, aes(x = Res_ID, y = RMSF, group = Dataset, color = Dataset)) + 
  geom_line() +  # Add line layer
  geom_point() + # Add points on each data point for emphasis
  labs(title = "RMSF per Residue for different ", x = "Residue Number", y = "RMSF", color = "Group") + # Adding labels
  theme_minimal()  # A minimalistic theme
print(line_plot)
ggsave(pngout, plot = line_plot, width = 6, height = 4, dpi = 300)
