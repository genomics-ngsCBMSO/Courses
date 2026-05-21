ctprop_test <- function(df, freqs) {

  wide <- reshape(df[, c("Sample", "CellType", "Proportion")],
                timevar = "Sample",
                idvar = "CellType",
                direction = "wide")
  
  colnames(wide) <- gsub("Proportion.", "", colnames(wide))  

  epsilon <- 1e-9

  wide$log2FC <- log2((wide$Donor1 + epsilon) / (wide$Donor2 + epsilon))
  results <- data.frame()

  for (ct in colnames(freqs)) {
  
      x1 <- freqs["Donor1", ct]
      x2 <- freqs["Donor2", ct]
  
      n1 <- sum(freqs["Donor1", ])
      n2 <- sum(freqs["Donor2", ])
  
      p1 <- x1 / n1
      p2 <- x2 / n2
  
      p_pool <- (x1 + x2) / (n1 + n2)
  
      Z <- (p2 - p1) / sqrt(p_pool * (1 - p_pool) * (1/n1 + 1/n2))
  
      pval <- 2 * pnorm(-abs(Z))
  
     results <- rbind(results,
                   data.frame(CellType = ct,
                              Z = Z,
                              p = pval))
  }
  
  results$p_adj <- p.adjust(results$p, method = "fdr")
  
  results$stars <- ""

  results$stars[results$p_adj < 0.05] <- "*"
  results$stars[results$p_adj < 0.01] <- "**"
  results$stars[results$p_adj < 0.001] <- "***"
  
  final <- merge(wide, results, by = "CellType")
  
  return(final)
}
