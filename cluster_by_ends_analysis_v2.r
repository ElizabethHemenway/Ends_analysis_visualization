#!/usr/bin/env Rscript

### written by: Elizabeth Hemenway
### 11/5/2024
### A script to take the average of each ends analysis window across biological replicates, subtract mutant from wild-type, and cluster based on differences. 
### usage: cluster_by_ends_analysis.R /lab/solexa_gehring/elizabeth/TE_spreading/ends/endsout_afteranalysis/ 3 r3_minus_Col_r3TEs_CG listofdatafiles

library(ggplot2)
library(reshape2)
library(gplots)
library(plotly)
library(plyr)

library(dplyr)
library(pheatmap)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
print(args)
path <- args[1] #where you want output files to go
reps <- as.numeric(args[2]) # 2 or 3 biological reps supported
output <- args[3] #what you want file prefix to be

setwd(path)


#functions 
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# this function takes mutant and wild types separately, merges dataframes, fills any NA cells with the median of the column, and takes the average of each cell across biological replicates.
# then it takes the difference between the mutant average and the wild-type average, and clusters TEs using k means built into pheatmap. 
# read through this section closely, not all settings may be appropriate for your purposes. 
avg_diff_and_cluster <- function(df1, df2, df3, df4, df5, df6, reps=3, k=4, saveas) {
  # Combine the data frames
  if (reps == 2) {
    group1 <- rbind(df1,df2)
    group2 <- rbind(df4,df5)
  } else if (reps == 3) {
    group1 <- rbind(df1,df2,df3)
    group2 <- rbind(df4,df5,df6)
  } else {
    print ("not valid reps")}

  # Loop through each column starting from the second column
  for (i in 2:ncol(group1)) {
  # Replace NA values with the mean of the column
  group1[is.na(group1[, i]), i] <- median(group1[, i], na.rm = TRUE)
  }
  # Loop through each column starting from the second column
  for (i in 2:ncol(group2)) {
  # Replace NA values with the mean of the column
  group2[is.na(group2[, i]), i] <- median(group2[, i], na.rm = TRUE)
  }
  # Calculate mean of TE windows acrosss replicates
  group1_mean <- group1 %>%
    group_by(ID) %>%
    mutate(across(everything(), mean, na.rm = TRUE)) %>%
    distinct(ID, .keep_all = TRUE) %>%
    ungroup()
  group2_mean <- group2 %>%
    group_by(ID) %>%
    mutate(across(everything(), mean, na.rm = TRUE)) %>%
    distinct(ID, .keep_all = TRUE) %>%
    ungroup()
    
  group2_mean_neg <- group2_mean %>%
    mutate(across(where(is.numeric), ~ . * -1))

  # Combine the two data frames
  groups <- bind_rows(group1_mean, group2_mean_neg)
  drop <- c("feature_start", "feature_end", "feature_mid", "feature_mid.1","feature_mid.2", "feature_mid.3") 

  # Remove unwanted columns
  groups <- select(groups, -one_of(drop))

  # Calculate the sum of numeric columns grouped by the first column
  diff <- groups %>%
    group_by(across(1)) %>%  # Group by the first column
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  drop <- c("ID")
  diffdrop <- select(diff, -one_of(drop))

  #generate count matrix with full dataset
  z <- as.matrix((diffdrop))
  #center and scale
  scaledata <- t(scale(t(z)))

  out <- pheatmap(scaledata, color = colorRampPalette(rev(brewer.pal(n = 7, name =
    "RdYlBu")))(100), breaks = NA, border_color = "grey60",
    cellwidth = NA, cellheight = NA, kmeans_k = 4, scale = "none", cluster_rows = TRUE,
    cluster_cols = FALSE, clustering_distance_rows = "euclidean", clustering_distance_cols="euclidean", angle_col=45)
  save_pheatmap_pdf(out, saveas)

  clusterDF <- as.data.frame(factor(out$kmeans$cluster))
  colnames(clusterDF) <- "Cluster"
  clusterout <- cbind(clusterDF$Cluster, diff$ID)
  colnames(clusterout) <- c("Cluster", "ID")  # Rename columns
  return (clusterout)
}

med_diff_and_cluster <- function(df1, df2, df3, df4, df5, df6, reps=3, k=4, saveas) {
  # Combine the data frames
  if (reps == 2) {
    group1 <- rbind(df1,df2)
    group2 <- rbind(df4,df5)
  } else if (reps == 3) {
    group1 <- rbind(df1,df2,df3)
    group2 <- rbind(df4,df5,df6)
  } else {
    print ("not valid reps")}

  # Loop through each column starting from the second column
  for (i in 2:ncol(group1)) {
  # Replace NA values with the median of the column
  group1[is.na(group1[, i]), i] <- median(group1[, i], na.rm = TRUE)
  }
  # Loop through each column starting from the second column
  for (i in 2:ncol(group2)) {
  # Replace NA values with the mean of the column
  group2[is.na(group2[, i]), i] <- median(group2[, i], na.rm = TRUE)
  }
  # Calculate mean of TE windows acrosss replicates
  group1_med <- group1 %>%
    group_by(ID) %>%
    mutate(across(everything(), median, na.rm = TRUE)) %>%
    distinct(ID, .keep_all = TRUE) %>%
    ungroup()
  group2_med <- group2 %>%
    group_by(ID) %>%
    mutate(across(everything(), median, na.rm = TRUE)) %>%
    distinct(ID, .keep_all = TRUE) %>%
    ungroup()
    
  group2_med_neg <- group2_med %>%
    mutate(across(where(is.numeric), ~ . * -1))

  # Combine the two data frames
  groups <- bind_rows(group1_med, group2_med_neg)
  drop <- c("feature_start", "feature_end", "feature_mid", "feature_mid.1","feature_mid.2", "feature_mid.3") 

  # Remove unwanted columns
  groups <- select(groups, -one_of(drop))

  # Calculate the sum of numeric columns grouped by the first column
  diff <- groups %>%
    group_by(across(1)) %>%  # Group by the first column
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  drop <- c("ID")
  diffdrop <- select(diff, -one_of(drop))

  #generate count matrix with full dataset
  z <- as.matrix((diffdrop))
  #center and scale
  scaledata <- t(scale(t(z)))

  out <- pheatmap(scaledata, color = colorRampPalette(rev(brewer.pal(n = 7, name =
    "RdYlBu")))(100), breaks = NA, border_color = "grey60",
    cellwidth = NA, cellheight = NA, kmeans_k = 4, scale = "none", cluster_rows = TRUE,
    cluster_cols = FALSE, clustering_distance_rows = "euclidean", clustering_distance_cols="euclidean", angle_col=45)
  save_pheatmap_pdf(out, saveas)

  clusterDF <- as.data.frame(factor(out$kmeans$cluster))
  colnames(clusterDF) <- "Cluster"
  clusterout <- cbind(clusterDF$Cluster, diff$ID)
  colnames(clusterout) <- c("Cluster", "ID")  # Rename columns
  return (clusterout)
}

med_wt_and_cluster <- function(df4, df5, df6, reps=3, k=8, saveas) {
  # Combine the data frames
  if (reps == 2) {
    group2 <- rbind(df4,df5)
  } else if (reps == 3) {
    group2 <- rbind(df4,df5,df6)
  } else {
    print ("not valid reps")}

  # Loop through each column starting from the second column
  for (i in 2:ncol(group2)) {
  # Replace NA values with the mean of the column
  group2[is.na(group2[, i]), i] <- median(group2[, i], na.rm = TRUE)
  }
  # Calculate mean of TE windows acrosss replicates
  group2_med <- group2 %>%
    group_by(ID) %>%
    mutate(across(everything(), median, na.rm = TRUE)) %>%
    distinct(ID, .keep_all = TRUE) %>%
    ungroup()
    
  drop <- c("ID","feature_start", "feature_end", "feature_mid", "feature_mid.1","feature_mid.2", "feature_mid.3") 

  # Remove unwanted columns
  group2_drop <- select(group2_med, -one_of(drop))

  #generate count matrix with full dataset
  z <- as.matrix((group2_drop))
  #center and scale
  scaledata <- t(scale(t(z)))

  out <- pheatmap(scaledata, color = colorRampPalette(rev(brewer.pal(n = 7, name =
    "RdYlBu")))(100), breaks = NA, border_color = "grey60",
    cellwidth = NA, cellheight = NA, kmeans_k = 8, scale = "none", cluster_rows = TRUE,
    cluster_cols = FALSE, clustering_distance_rows = "euclidean", clustering_distance_cols="euclidean", angle_col=45)
  save_pheatmap_pdf(out, saveas)

  clusterDF <- as.data.frame(factor(out$kmeans$cluster))
  colnames(clusterDF) <- "Cluster"
  clusterout <- cbind(clusterDF$Cluster, group2_med$ID)
  colnames(clusterout) <- c("Cluster", "ID")  # Rename columns
  return (clusterout)
}

#this function is for after you've already clustered the data, and now want to make another heatmap with all of the TEs plotted indiviidually, and with the actual methylation data
avg_diff_and_predef_clust <- function(df1, df2, df3, df4, df5, df6, reps=3) {
  # Combine the data frames
  if (reps == 2) {
    group1 <- rbind(df1,df2)
    group2 <- rbind(df4,df5)
  } else if (reps == 3) {
    group1 <- rbind(df1,df2,df3)
    group2 <- rbind(df4,df5,df6)
  } else {
    print ("not valid reps")}

  # Loop through each column starting from the second column
  for (i in 2:ncol(group1)) {
  # Replace NA values with the mean of the column
  group1[is.na(group1[, i]), i] <- median(group1[, i], na.rm = TRUE)
  }
  # Loop through each column starting from the second column
  for (i in 2:ncol(group2)) {
  # Replace NA values with the mean of the column
  group2[is.na(group2[, i]), i] <- median(group2[, i], na.rm = TRUE)
  }
  # Calculate mean of TE windows acrosss replicates
  group1_mean <- group1 %>%
    group_by(ID) %>%
    mutate(across(everything(), mean, na.rm = TRUE)) %>%
    distinct(ID, .keep_all = TRUE) %>%
    ungroup()
  group2_mean <- group2 %>%
    group_by(ID) %>%
    mutate(across(everything(), mean, na.rm = TRUE)) %>%
    distinct(ID, .keep_all = TRUE) %>%
    ungroup()
    
  group2_mean_neg <- group2_mean %>%
    mutate(across(where(is.numeric), ~ . * -1))

  # Combine the two data frames
  groups <- bind_rows(group1_mean, group2_mean_neg)
  drop <- c("feature_start", "feature_end", "feature_mid", "feature_mid.1","feature_mid.2", "feature_mid.3") 
    
  # Remove unwanted columns
  groups <- select(groups, -one_of(drop))
  # Calculate the sum of numeric columns grouped by the first column
  diff <- groups %>%
    group_by(across(1)) %>%  # Group by the first column
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  drop <- c("ID")
  diffdrop <- select(diff, -one_of(drop))

  #generate count matrix with full dataset
  z <- as.matrix((diff))
  return(z)
  #center and scale
  scaledata <- t(scale(t(z)))
  return (scaledata)
}

med_diff_and_predef_clust <- function(df1, df2, df3, df4, df5, df6, reps=3) {
  # Combine the data frames
  if (reps == 2) {
    group1 <- rbind(df1,df2)
    group2 <- rbind(df4,df5)
  } else if (reps == 3) {
    group1 <- rbind(df1,df2,df3)
    group2 <- rbind(df4,df5,df6)
  } else {
    print ("not valid reps")}

  # Loop through each column starting from the second column
  for (i in 2:ncol(group1)) {
  # Replace NA values with the mean of the column
  group1[is.na(group1[, i]), i] <- median(group1[, i], na.rm = TRUE)
  }
  # Loop through each column starting from the second column
  for (i in 2:ncol(group2)) {
  # Replace NA values with the mean of the column
  group2[is.na(group2[, i]), i] <- median(group2[, i], na.rm = TRUE)
  }
  # Calculate mean of TE windows acrosss replicates
  group1_mean <- group1 %>%
    group_by(ID) %>%
    mutate(across(everything(), median, na.rm = TRUE)) %>%
    distinct(ID, .keep_all = TRUE) %>%
    ungroup()
  group2_mean <- group2 %>%
    group_by(ID) %>%
    mutate(across(everything(), median, na.rm = TRUE)) %>%
    distinct(ID, .keep_all = TRUE) %>%
    ungroup()
    
  group2_mean_neg <- group2_mean %>%
    mutate(across(where(is.numeric), ~ . * -1))

  # Combine the two data frames
  groups <- bind_rows(group1_mean, group2_mean_neg)
  drop <- c("feature_start", "feature_end", "feature_mid", "feature_mid.1","feature_mid.2", "feature_mid.3") 
    
  # Remove unwanted columns
  groups <- select(groups, -one_of(drop))
  # Calculate the sum of numeric columns grouped by the first column
  diff <- groups %>%
    group_by(across(1)) %>%  # Group by the first column
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  drop <- c("ID")
  diffdrop <- select(diff, -one_of(drop))

  #generate count matrix with full dataset
  z <- as.matrix((diff))
  return(z)
  #center and scale
  scaledata <- t(scale(t(z)))
  return (scaledata)
}

med_wt_predef_clust <- function(df4, df5, df6, reps=3) {
  # Combine the data frames
  if (reps == 2) {
    group2 <- rbind(df4,df5)
  } else if (reps == 3) {
    group2 <- rbind(df4,df5,df6)
  } else {
    print ("not valid reps")}

  # Loop through each column starting from the second column
  for (i in 2:ncol(group2)) {
  # Replace NA values with the mean of the column
  group2[is.na(group2[, i]), i] <- median(group2[, i], na.rm = TRUE)
  }
  # Calculate mean of TE windows acrosss replicates
  group2_mean <- group2 %>%
    group_by(ID) %>%
    mutate(across(everything(), median, na.rm = TRUE)) %>%
    distinct(ID, .keep_all = TRUE) %>%
    ungroup()
    
  # Combine the two data frames
  drop <- c("ID","feature_start", "feature_end", "feature_mid", "feature_mid.1","feature_mid.2", "feature_mid.3") 

  # Remove unwanted columns
  group2_drop <- select(group2_med, -one_of(drop))

  #generate count matrix with full dataset
  z <- as.matrix((diff))
  return(z)
  #center and scale
  scaledata <- t(scale(t(z)))
  return (scaledata)
}


if (reps == 2) {
  print("2 replicates each sample")
  m1 <- read.table(args[4], sep="\t", header=1)
  m2 <- read.table(args[5], sep="\t", header=1)
  wt1 <- read.table(args[6], sep="\t", header=1)
  wt2 <- read.table(args[7], sep="\t", header=1)

  m3 <- NA
  wt3 <- NA
  
} else if (reps == 3) {
  print("3 replicates each sample")
  m1 <- read.table(args[4], sep="\t", header=1)
  m2 <- read.table(args[5], sep="\t", header=1)
  m3 <- read.table(args[6], sep="\t", header=1)

  wt1 <- read.table(args[7], sep="\t", header=1)
  wt2 <- read.table(args[8], sep="\t", header=1)
  wt3 <- read.table(args[9], sep="\t", header=1)

} else { 
  print("not an accepted number of replicates")
}

#print(m1)
m_minus_wt_clusters <- med_diff_and_cluster(m1, m2, m3, wt1, wt2, wt3, reps=reps, k=4, saveas=paste0(output, "_clusters_hmap.pdf")) # change function here as needed

m_minus_wt_clust_df <- as.data.frame(m_minus_wt_clusters)
rownames(m_minus_wt_clust_df) <- m_minus_wt_clust_df$ID

clustanno <- select(m_minus_wt_clust_df, -one_of("ID"))
write.csv(clustanno, paste0(output, "_clusters.csv"))

m_minus_wt_diff <- med_diff_and_predef_clust(m1, m2, m3, wt1, wt2, wt3, reps=reps) # change function here as needed


merged_data <- merge(m_minus_wt_diff, m_minus_wt_clust_df, by = "ID")

merged_data <- merged_data[order(merged_data$Cluster), ]
rownames(merged_data) <- merged_data$ID
row_names <- rownames(merged_data)
merged_data <-select(merged_data, -one_of("ID", "Cluster"))


merged_data<-as.matrix((merged_data))

#the actual mC difference data
num_matrix <- matrix(as.numeric(merged_data), nrow = nrow(merged_data))

rownames(num_matrix) <- row_names

scaledata <- t(scale(t(num_matrix)))

#mat_breaks <- seq(-10, 100, length.out = 99) # CG
#mat_breaks <- seq(-10, 100, length.out = 99) # CHG
#mat_breaks <- seq(-10, 10, length.out = 99) # CHH

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
#paletteLength <- 100

#myBreaks <- c(seq(min(num_matrix), 0, length.out=ceiling(paletteLength/2) + 1), 
             # seq(max(num_matrix)/paletteLength, max(num_matrix), length.out=floor(paletteLength/2))) # works well for mCG

paletteLength <- 40

myBreaks <- c(seq(-10, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(10/paletteLength, 10, length.out=floor(paletteLength/2))) # works well for mCHH


map <- pheatmap(num_matrix, color = colorRampPalette(rev(brewer.pal(n = 7, name ="PiYG")))(40),
         breaks = myBreaks, border_color = "grey60",
         cellwidth = NA, cellheight = NA,
         scale = "none", cluster_rows = FALSE, 
         cluster_cols = FALSE, annotation_row=clustanno, show_rownames = F)
save_pheatmap_pdf(map, paste0(output,"_mC_diff_value_cluster_hmap.pdf"))


