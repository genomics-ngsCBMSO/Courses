#!Rscript

library(ggplot2)
library(zoo)

# Check if the command line argument is provided
if (length(commandArgs(trailingOnly = TRUE)) < 3 ) {
  stop("Usage: Rscript --vanilla script.R ONT.cov Illumina.cov out.png")
}

# Get the input .cov file paths
file1 <- commandArgs(trailingOnly = TRUE)[1]
file2 <- commandArgs(trailingOnly = TRUE)[2]
output_file <- commandArgs(trailingOnly = TRUE)[3]

# Read the data from the first .cov file
coverage_data <- read.table(file1, header = FALSE)
colnames(coverage_data) <- c("chromosome", "position", "coverage")

# Read the data from the second .cov file
coverage_data2 <- read.table(file2, header = FALSE)
colnames(coverage_data2) <- c("chromosome", "position", "coverage")

# Smooth the coverage data
smoothed_coverage <- rollmean(coverage_data$coverage, k = 200, align = "center", fill = NA)
smoothed_data <- data.frame(position = coverage_data$position, smoothed_coverage = log(smoothed_coverage))

# Smooth the coverage data from the second .cov file
smoothed_coverage2 <- rollmean(coverage_data2$coverage, k = 200, align = "center", fill = NA)
smoothed_data2 <- data.frame(position = coverage_data2$position, smoothed_coverage2 = log(smoothed_coverage2))

# Create the plot with ggplot2
p <- ggplot(smoothed_data, aes(x = position, y = smoothed_coverage)) +
  geom_line(aes(color = "ONT Reads"), size = 0.2) +
  geom_line(data = smoothed_data2, aes(x = position, y = smoothed_coverage2, color = "Illumina Reads"), size = 0.2) +
  labs(x = "Chromosome position", y = "Smooth coverage") +
  scale_color_manual(values = c("ONT Reads" = "blue", "Illumina Reads" = "red")) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 4),
    axis.title = element_text(size = 4),
    axis.text = element_text(size = 4)
  )
png(output_file, width = 1600, height = 600, res = 300)
print(p)
dev.off()  
