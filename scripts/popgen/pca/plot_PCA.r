# Load the wordcloud package for better text positioning
library(wordcloud)

# Load the data
pca <- read.table("mpusilla.snps.eigenval")
eigen <- read.table("mpusilla.snps.eigenvec")

# Calculate variance explained
total_variance <- sum(pca$V1)
pc1_var <- (pca$V1[1] / total_variance) * 100
pc2_var <- (pca$V1[2] / total_variance) * 100
pc3_var <- (pca$V1[3] / total_variance) * 100

# Define group2 samples
group2_samples <- c("RCC1749", "RCC3052")

# Create color vector with transparency (alpha)
colors_transparent <- ifelse(eigen$V2 %in% group2_samples, 
                            rgb(1, 0, 0, alpha=0.7),    # Semi-transparent red
                            rgb(0, 0, 1, alpha=0.7))    # Semi-transparent blue

# Create the PDF
pdf(file="pca.pdf", height=6, width=12)
par(mfrow=c(1,2))

# Calculate expanded axis limits to accommodate labels
x_range <- range(eigen$V3)
y_range_pc2 <- range(eigen$V4)
y_range_pc3 <- range(eigen$V5)

# Expand limits by 20% on each side for even more label space
x_expand <- diff(x_range) * 0.20
y_expand_pc2 <- diff(y_range_pc2) * 0.20
y_expand_pc3 <- diff(y_range_pc3) * 0.20

xlim_expanded <- c(x_range[1] - x_expand, x_range[2] + x_expand)
ylim_pc2_expanded <- c(y_range_pc2[1] - y_expand_pc2, y_range_pc2[2] + y_expand_pc2)
ylim_pc3_expanded <- c(y_range_pc3[1] - y_expand_pc3, y_range_pc3[2] + y_expand_pc3)

# Plot PC1 vs PC2 with expanded limits
plot(eigen$V3, eigen$V4, 
     xlab=paste0("PC1 (", round(pc1_var, 1), "%)"), 
     ylab=paste0("PC2 (", round(pc2_var, 1), "%)"), 
     col=colors_transparent, 
     pch=16, 
     cex=1.2,
     xlim=xlim_expanded,
     ylim=ylim_pc2_expanded)

# Add a slight jitter to overlapping points to separate them visually
jitter_amount <- 0.005  # Adjust this value as needed
eigen_jittered_x <- jitter(eigen$V3, amount=jitter_amount)
eigen_jittered_y_pc2 <- jitter(eigen$V4, amount=jitter_amount)
eigen_jittered_y_pc3 <- jitter(eigen$V5, amount=jitter_amount)

# Redraw with jittered positions and transparency
points(eigen_jittered_x, eigen_jittered_y_pc2, col=colors_transparent, pch=16, cex=1.2)

# Use textplot from wordcloud for non-overlapping labels (use original positions)
textplot(eigen$V3, eigen$V4, 
         words=eigen$V2, 
         cex=0.6,
         new=FALSE,
         show.lines=TRUE,
         xlim=xlim_expanded,
         ylim=ylim_pc2_expanded)

# Plot PC1 vs PC3 with expanded limits
plot(eigen$V3, eigen$V5, 
     xlab=paste0("PC1 (", round(pc1_var, 1), "%)"), 
     ylab=paste0("PC3 (", round(pc3_var, 1), "%)"), 
     col=colors_transparent, 
     pch=16, 
     cex=1.2,
     xlim=xlim_expanded,
     ylim=ylim_pc3_expanded)

# Redraw with jittered positions and transparency
points(eigen_jittered_x, eigen_jittered_y_pc3, col=colors_transparent, pch=16, cex=1.2)

# Use textplot for non-overlapping labels
textplot(eigen$V3, eigen$V5, 
         words=eigen$V2, 
         cex=0.6,
         new=FALSE,
         show.lines=TRUE,
         xlim=xlim_expanded,
         ylim=ylim_pc3_expanded)

dev.off()

# Print some useful info
cat("Total samples:", nrow(eigen), "\n")
cat("Group2 samples found:", sum(eigen$V2 %in% group2_samples), "\n")
cat("PC1 variance explained:", round(pc1_var, 2), "%\n")
cat("PC2 variance explained:", round(pc2_var, 2), "%\n")
cat("PC3 variance explained:", round(pc3_var, 2), "%\n")