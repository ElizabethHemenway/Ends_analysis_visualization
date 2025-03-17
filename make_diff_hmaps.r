#!/usr/bin/env Rscript

### written by: Elizabeth Hemenway with help from chatgpt
### 1/20/25
### A script to take the average of each ends analysis window across biological replicates, subtract mutant from wild-type, and cluster based on differences. 
### usage: ./make_diff_hmaps.r r3_minus_Col_r3TEs_CG listofdatafiles

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
output <- args[2] #what you want file prefix to be

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

med_diff <- function(df1, df2, df3, df4, df5, df6) {
  # Combine the data frames
  group1 <- rbind(df1, df2, df3)
  group2 <- rbind(df4, df5, df6)

  # Compute medians for group1
  group1_med <- group1 %>%
    group_by(ID) %>%
    mutate(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
    distinct(ID, .keep_all = TRUE) %>%
    ungroup()

  # Compute medians for group2
  group2_med <- group2 %>% 
    group_by(ID) %>%
    mutate(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
    distinct(ID, .keep_all = TRUE) %>%
    ungroup()
    
  # Subtract the medians of group2 from group1
  med_diff_df <- group1_med %>%
    left_join(group2_med, by = "ID", suffix = c("_group1", "_group2")) %>%
    mutate(across(ends_with("_group2"), ~ .x * -1)) %>%  # Make group2 values negative
    mutate(across(ends_with("_group1"), ~ if_else(
      is.na(.x) | is.na(get(sub("_group1$", "_group2", cur_column()))), 
	  NA_real_,  # Output NA if either value is NA
	  .x + get(sub("_group1$", "_group2", cur_column()))  # Otherwise, perform the subtraction
	  ))) %>%
    select(ID, ends_with("_group1"))  # Keep only the necessary columns
  # Return the final dataframe with the medians difference
  return(med_diff_df)
}
  
#read in data
m1 <- read.table(args[3], sep="\t", header=1)
m2 <- read.table(args[4], sep="\t", header=1)
m3 <- read.table(args[5], sep="\t", header=1)

wt1 <- read.table(args[6], sep="\t", header=1)
wt2 <- read.table(args[7], sep="\t", header=1)
wt3 <- read.table(args[8], sep="\t", header=1)

#make difference dataframe to plot
m_minus_wt_diff <- med_diff(m1, m2, m3, wt1, wt2, wt3) 

# Remove unwanted columns
#drops <- c("feature_start_group1", "feature_end_group1", "feature_mid_group1", "feature_mid.1_group1","feature_mid.2_group1", "feature_mid.3_group1")     

#m_minus_wt_diff <- select(m_minus_wt_diff, -one_of(drops))

#write.csv(m_minus_wt_diff, paste0(output, "plotted_df.csv"))  ## this is the dataframe I've included in my email. 

# Now, ensure that only numeric columns are included in the matrix - but I want to keep NAs as place holder? 
lastcol=ncol(m_minus_wt_diff)
m_minus_wt_diff[,2:lastcol] <- sapply(m_minus_wt_diff[,2:lastcol], as.numeric)
# Scale the data (center and scale each row - need to do transformations to get to rows and back.) - doesn't work with NAs, but I think I should be scaling prior to clustering?
#scaledata <- t(scale(t(m_minus_wt_diff_mat)))

#generate palette breaks if needed for better visualization
maxdiff <- strtoi(args[9])
print(maxdiff)
paletteLength <- 40

myBreaks <- c(seq(-(maxdiff), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(10/paletteLength, maxdiff, length.out=floor(paletteLength/2)))


# Then be sure only to use these columns as input for pheatmap.  One way to do so is like this:
map <- pheatmap(as.matrix(m_minus_wt_diff[,2:lastcol]), color = colorRampPalette(rev(brewer.pal(n = 7, name ="PiYG")))(40),border_color = "grey60",
          breaks=myBreaks, cellwidth = NA, cellheight = NA,
          scale = "none", cluster_rows = TRUE,
          cluster_cols = FALSE, show_rownames = F)
save_pheatmap_pdf(map, paste0(output,"_mC_diff_value_hmap.pdf"))
# Capture the row order from the clustering dendrogram of the first heatmap
row_order <- map$tree_row$order

# Reorder the diff dataframe based on the clustering order
m_minus_wt_diff_reordered <- m_minus_wt_diff[row_order, ]


# make a plot of just wt methylation
med <- function(df1, df2, df3) {
  # Combine the data frames
  group1 <- rbind(df1, df2, df3)

  # Compute medians for group1
  group1_med <- group1 %>%
    group_by(ID) %>%
    mutate(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
    distinct(ID, .keep_all = TRUE) %>%
    ungroup()
  
  # Return the final dataframe with the medians difference
  return(group1_med)
}

#make wt median dataframe to plot
wt_med <- med(wt1, wt2, wt3) 
wt_med[,2:lastcol] <- sapply(wt_med[,2:lastcol], as.numeric)

# Reorder the wt_med dataframe based on the clustering order
wt_med_order <- wt_med[match(m_minus_wt_diff_reordered$ID, wt_med$ID), ]


map <- pheatmap(as.matrix(wt_med_order[,2:lastcol]), color = colorRampPalette(brewer.pal(n = 7, name ="RdPu"))(20),
          border_color = "grey60",
          cellwidth = NA, cellheight = NA,
          scale = "none", cluster_rows = FALSE,
          cluster_cols = FALSE, show_rownames = F)
save_pheatmap_pdf(map, paste0(output,"_mC_wt_median_value_hmap.pdf"))

#make mutant median dataframe to plot
r3_med <- med(m1, m2, m3) 
r3_med[,2:lastcol] <- sapply(r3_med[,2:lastcol], as.numeric)

r3_med_order <- r3_med[match(m_minus_wt_diff_reordered$ID, r3_med$ID), ]

map <- pheatmap(as.matrix(r3_med_order[,2:lastcol]), color = colorRampPalette(brewer.pal(n = 7, name ="RdPu"))(20),
          border_color = "grey60",
          cellwidth = NA, cellheight = NA,
          scale = "none", cluster_rows = FALSE,
          cluster_cols = FALSE, show_rownames = F)
save_pheatmap_pdf(map, paste0(output,"_mC_r3_median_value_hmap.pdf"))

