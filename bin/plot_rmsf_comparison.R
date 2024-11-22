#!/usr/bin/env Rscript
library(docopt)
library(ggplot2)
library(dplyr)
#library(bio3d)
doc <- "Usage:
plot_rmsf_individual.R --incsv=<csv> --wtcsv=<wtcsv>
Options:
--incsv=<csv>  load input csv file
--wtcsv=<wtcsv> load input comparison  csv file
"
opts <- docopt(doc)
csv <- opts$`--incsv`
wtcsv <- opts$`--wtcsv`
data <- read.csv(csv)
data_wt <- read.csv(wtcsv)
shortstem <- sub("4.0pH.csv", "", csv)
stem <- sub(".csv", "", csv)
stem_wt <- sub(".csv", "", wtcsv)
pngout <- paste0(shortstem, "_RMSF_pHcomp.pdf")
titlecard <- paste(shortstem, "RMSF per Residue, pH comparison")

random_color <- function() { rgb(runif(1), runif(1), runif(1))}

combined_data <- rbind(
  data.frame(dataset = stem, data),
  data.frame(dataset = stem_wt, data_wt)
)

# Base plot
line_plot <- ggplot(combined_data, aes(x = Res_ID, y = RMSF, color = dataset)) +
  geom_line() +
  geom_point() +
  coord_cartesian(xlim=c(0, 398), ylim=c(0, 5)) +
  labs(title = titlecard, x = "Res_ID", y = "RMSF", color = "dataset") +
  theme_minimal()
ggsave(pngout, plot = line_plot, width = 20, height = 4, dpi = 300)
#line_plot <- ggplot() +
#  geom_line(data = data, aes(x = Res_ID, y = RMSF), color = "Dataset 1") +
#  geom_point(data = data, aes(x = Res_ID, y = RMSF), color = "Dataset 1") +
#  geom_line(data = data_wt, aes(x = Res_ID, y = RMSF), color = "Dataset 2") +
#  geom_point(data = data_wt, aes(x = Res_ID, y = RMSF), color = "Dataset 2") +
#  scale_color_manual(name = "Dataset", values = c("Dataset 1" = "blue", "Dataset 2" = "red")) +
#  labs(title = titlecard, x = "Res_ID", y = "RMSF")
#ggsave(pngout, plot = line_plot, width = 20, height = 4, dpi = 300)

#line_plot <- ggplot(data, aes(x = Res_ID, y = RMSF, )) + 
#  geom_line() +  # Add line layer
#  geom_point() + # Add points on each data point for emphasis
#  labs(title = titlecard, x = "Residue Number", y = "RMSF", color = random_color()) + # Adding labels
#  theme_minimal()  # A minimalistic theme
#print(line_plot)
#ggsave(pngout, plot = line_plot, width = 20, height = 4, dpi = 300)
