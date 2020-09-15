# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# Functions
read_parseBcfStat <- function(path) {
  stat <- read.table(path, header = T, sep = '\t')
}

theme_fj <- 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "grey", fill=NA, size=1))

# Set Arguments
args = commandArgs(trailingOnly = T)
stat_file = args[1]

# Set paths
sn_path   = paste(stat_file, 'SN', sep = '.')
dp_path   = paste(stat_file, 'DP', sep = '.')
qual_path = paste(stat_file, 'QUAL', sep = '.')

# Read tables
sn <- read_parseBcfStat(sn_path) 

dp <- 
  read_parseBcfStat(dp_path) %>%
  select(DP = bin, `Number of sites` = number.of.sites, 
         `Fraction of sites (%)` = fraction.of.sites....)

qual <- 
  read_parseBcfStat(qual_path) %>%
  select(QUAL = Quality, SNPs = number.of.SNPs, Indels = number.of.indels) %>%
  pivot_longer(cols = SNPs:Indels, names_to = "Variant", values_to = "Count")

# Plot
p1 <- 
  ggplot(qual, aes(QUAL, Count, col = Variant)) +
  geom_point(size =1, alpha = 0.4) +
  theme_fj +
  scale_y_log10() +
  scale_x_log10() + 
  annotation_logticks(sides = 'bl', colour = 'grey') +
  theme(legend.position = c(0.2, 0.8),
        legend.title = element_blank())
p2 <- 
  ggplot(dp, aes(DP, `Number of sites`)) +
  geom_point(size =1, alpha = 0.4) +
  theme_fj +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(sides = 'bl', color = 'grey')

# Print to pdf
out = paste(args[1], ".pdf", sep = '')
pdf(out, width = 8, height = 4.5)
plot_grid(p1, p2)
dev.off()

out_mes = paste("Plots written to ", out, sep = '')
print(out_mes)