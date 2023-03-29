#!Rscript

library(argparse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(patchwork)

parser <- ArgumentParser(description='')
parser$add_argument('pdf',help='name for output pdf file.')
parser$add_argument('--bin',dest='bin',help='bin annotation file')
parser$add_argument('--cell',dest='cell',help='cell label')
parser$add_argument('--hic',dest='hic',help='hic pair file')
parser$add_argument('--edge',dest='edge',help='edge file (SPIN input)')
parser$add_argument('--hic_merge_factor',dest='hic_merge_factor',type='integer',help='whether to merge hic bin, eg. 1 means do not merge, 10 means merge 10 bins into a big bin')
parser$add_argument('--chrom_1',dest='chrom_1',help='first chromosome')
parser$add_argument('--start_1',type='integer',dest='start_1',help='first start position')
parser$add_argument('--end_1',type='integer',dest='end_1',help='first end position')
parser$add_argument('--chrom_2',dest='chrom_2',help='first chromosome')
parser$add_argument('--start_2',type='integer',dest='start_2',help='first start position')
parser$add_argument('--end_2',type='integer',dest='end_2',help='first end position')
parser$add_argument('--pdf_height',dest='height',type='integer',help='pdf figure height in inch')
parser$add_argument('--pdf_width',dest='width',type='integer',help='pdf figure width in inch')


# parse argument
args <- parser$parse_args()
bin_file = args$bin
cell_label = args$cell
hic_pair_file = args$hic
hic_merge_factor = args$hic_merge_factor
edge_file = args$edge
chrom_1 = args$chrom_1
start_1 = as.numeric(args$start_1)
end_1 = as.numeric(args$end_1)
chrom_2 = args$chrom_2
start_2 = as.numeric(args$start_2)
end_2 = as.numeric(args$end_2)
pdf_file = args$pdf
pdf_width = args$width
pdf_height = args$height

## funtions
load_bin_anno <- function(filename, cell_label) {
  anno = read.table(file = filename, sep = '\t', header = F, stringsAsFactors = F)
  colnames(anno) <- c("chrom", "start", "end", "idx", "cell")
  anno <- subset(anno, cell == cell_label)
  # fix index starting from 0
  anno$idx = anno$idx + 1
  #idx_list = as.character(sort(as.numeric(anno$idx), decreasing = F))
  #anno$idx <- factor(anno$idx, levels = idx_list)
  return(anno)
}
parse_edge_by_region <- function(filename, bin_anno, query_chrom_1, query_start_1, query_end_1, query_chrom_2, query_start_2, query_end_2) {
  # load all edge
  data = read.table(file = filename, sep = '\t', header = F, stringsAsFactors = F)
  colnames(data) <- c("idx_1", "idx_2", "score")
  # fix index starting from 0
  data$idx_1 = data$idx_1 + 1
  data$idx_2 = data$idx_2 + 1
  if (max(data$idx_1) > max(as.numeric(bin_anno$idx)) | max(data$idx_2) > max(as.numeric(bin_anno$idx))) {
    error(paste0("bin id of Hi-C edge cannot be found in bin annotation"))
  }
  # filter
  idx_list_1 <- subset(bin_anno, (chrom == query_chrom_1) & (start >= query_start_1) & (end < query_end_1))[, 'idx']
  idx_list_2 <- subset(bin_anno, (chrom == query_chrom_2) & (start >= query_start_2) & (end < query_end_2))[, 'idx']
  print(paste0(length(idx_list_1), " bins are selected for the first chromosome"))
  print(paste0(length(idx_list_2), " bins are selected for the second chromosome"))  
  data_subset <- subset(data, idx_1 %in% idx_list_1 & idx_2 %in% idx_list_2)
  data_subset <- data_subset %>% distinct()
  print(paste0("After filtering ", nrow(data_subset), " edges left from ", nrow(data), " edges"))
  # build full matrix
  data_merge <- data_subset
  data_merge$count <- ifelse(data_merge$score > 0, 1, 0)
  # if (query_chrom_1 == query_chrom_2 & query_start_1 == query_start_2 & query_end_1 == query_end_2) {
  #   data_subset_rev <- data.frame(idx_1 = data_subset$idx_2, idx_2 = data_subset$idx_1, score = data_subset$score)
  #   data_merge <- rbind(data_subset, data_subset_rev)
  #   data_merge$count <- ifelse(data_merge$score > 0, 1, 0)
  # } else {
  #   data_merge <- data_subset
  #   data_merge$count <- ifelse(data_merge$score > 0, 1, 0)
  # }
  return(data_merge) 
}
parse_pair_by_region <- function(filename, query_chrom_1, query_start_1, query_end_1, query_chrom_2, query_start_2, query_end_2) {
  # load all pairs
  data = read.table(file = filename, sep = '\t', header = F, stringsAsFactors = F)
  colnames(data) <- c("start_1", "start_2", "score")
  # filter
  data_subset <- data %>%
    filter(start_1 >= query_start_1 & start_1 < query_end_1 & start_2 >= query_start_2 & start_2 < query_end_2)
  print(paste0("After filtering ", nrow(data_subset), " edges left from ", nrow(data), " edges"))
  #
  if (query_chrom_1 == query_chrom_2 & query_start_1 == query_start_2 & query_end_1 == query_end_2) {
    data_subset_rev <- data.frame(start_1 = data_subset$start_2, start_2 = data_subset$start_1, score = data_subset$score)
    data_merge <- rbind(data_subset, data_subset_rev) 
  } else {
    data_merge <- data_subset
  }
  return(data_merge)
}
set_axis_label <- function(axis_min, axis_max, spin_res) {
  if ( (axis_max - axis_min) >= 50000/spin_res*1000 - 0.1) {
    break_list = seq(axis_min, axis_max, 25000/spin_res*1000)
    label_list = paste0(seq(axis_min * spin_res / 1000000, axis_max * spin_res / 1000000, 25), "Mb")
  } else if ( (axis_max - axis_min) >= 10000/spin_res*1000 - 0.1) {
    break_list = seq(axis_min, axis_max, 10000/spin_res*1000)
    label_list = paste0(seq(axis_min * spin_res / 1000000, axis_max * spin_res / 1000000, 10), "Mb")
  } else if ((axis_max - axis_min) >= 5000/spin_res*1000 - 0.1) {
    break_list = seq(axis_min, axis_max, 5000/spin_res*1000)
    label_list = paste0(seq(axis_min * spin_res / 1000000, axis_max * spin_res / 1000000, 5), "Mb")
  } else {
    break_list = seq(axis_min, axis_max, 2500/spin_res*1000 - 0.1)
    label_list = paste0(seq(axis_min * spin_res / 1000000, axis_max * spin_res / 1000000, 2.5), "Mb")
  }
  return(list("break" = break_list, "label" = label_list))
}
plot_all_data <- function(bin_anno, hic_pair, edge_file, hic_merg_factor, spin_res, query_chrom_1, query_start_1, query_end_1, query_chrom_2, query_start_2, query_end_2) {
  # test parameters
  #hic_merge_factor = 1
  #query_chrom_1 = "chr1"
  #query_start_1 = 40000000
  #query_end_1 = 50000000
  #query_chrom_2 = "chr1"
  #query_start_2 = 40000000
  #query_end_2 = 50000000
  #hic_pair_file = "pre-process/K562_chr1_chr1_oe_VC_SQRT_25000.pair"
  #edge_file = "pre-process/K562.edge"
  #spin_res = 25000
  # extract edges
  spin_edge <- parse_edge_by_region(edge_file, bin_anno, query_chrom_1, query_start_1, query_end_1, query_chrom_2, query_start_2, query_end_2)
  hic_edge <- parse_pair_by_region(hic_pair_file, query_chrom_1, query_start_1, query_end_1, query_chrom_2, query_start_2, query_end_2)
  # add chromosome position
  spin_edge$start_1 <- bin_anno[match(spin_edge$idx_1, bin_anno$idx), "start"] 
  spin_edge$start_2 <- bin_anno[match(spin_edge$idx_2, bin_anno$idx), "start"] 
  # down sample edge
  spin_edge_ds <- spin_edge %>%
    mutate(start_1= hic_merge_factor * round(start_1 / spin_res / hic_merge_factor),
           start_2 = hic_merge_factor * round(start_2 / spin_res / hic_merge_factor)) %>%
    group_by(start_1, start_2) %>%
    summarise(count = n())
  # down-sample hi-c pair
  hic_edge_ds <- hic_edge %>%
    mutate(start_1= hic_merge_factor * round(start_1 / spin_res / hic_merge_factor),
           start_2 = hic_merge_factor * round(start_2 / spin_res / hic_merge_factor)) %>%
    group_by(start_1, start_2) %>%
    summarise(score = mean(score, na.rm = T), count = n())
  # set axis label
  axis_1 <- set_axis_label(min(hic_edge_ds$start_1), max(hic_edge_ds$start_1), spin_res)
  axis_2 <- set_axis_label(min(hic_edge_ds$start_2), max(hic_edge_ds$start_2), spin_res)
  # 
  plot_hic <- ggplot(hic_edge_ds, aes(x = start_2, y = start_1, fill = score )) +
    geom_raster(width = hic_merge_factor, height = hic_merge_factor) +
    #coord_fixed(ratio = 1) +
    xlab("") +
    ylab("") +
    scale_x_continuous(expand = c(0,0), limits = c(min(hic_edge_ds$start_2), max(hic_edge_ds$start_2)), breaks = as.numeric(unlist(axis_2["break"])), labels = as.list(unlist(axis_2["label"]))) +
    scale_y_continuous(expand = c(0,0), limits = c(min(hic_edge_ds$start_1), max(hic_edge_ds$start_1)), breaks = as.numeric(unlist(axis_1["break"])), labels = as.list(unlist(axis_1["label"]))) +
    #scale_fill_distiller(name = "count", type= "seq", palette = 'Purples', direction = 1) +
    scale_fill_viridis_c(name = "Average Hi-C score", option = "A", direction = -1, begin = 0, end = 0.9, limits = c(0, 40), oob = scales::squish) +
    theme_cowplot() +
    theme(legend.position = 'bottom',
          legend.margin  = margin(t = -1, unit = 'cm'), 
          legend.text = element_text(size = 6),
        legend.title = element_text(size = 8))
  plot_edge <- ggplot(spin_edge_ds, aes(x = start_2, y = start_1, fill = count / hic_merge_factor / hic_merge_factor)) +
    geom_tile(width = hic_merge_factor, height = hic_merge_factor) +
    #coord_fixed(ratio = 1) +
    xlab("") +
    ylab("") +
    scale_x_continuous(expand = c(0,0), limits = c(min(hic_edge_ds$start_2), max(hic_edge_ds$start_2)), breaks = as.numeric(unlist(axis_2["break"])), labels = as.list(unlist(axis_2["label"]))) +
    scale_y_continuous(expand = c(0,0), limits = c(min(hic_edge_ds$start_1), max(hic_edge_ds$start_1)), breaks = as.numeric(unlist(axis_1["break"])), labels = as.list(unlist(axis_1["label"]))) +
    #scale_fill_distiller(name = "count", type= "seq", palette = 'Purples', direction = 1) +
    scale_fill_viridis_c(name = "Average edge count per 25kb bin", option = "A", direction = -1, begin = 0, end = 0.5) +
    theme_cowplot() +
    theme(legend.position = 'bottom',
          legend.margin  = margin(t = -1, unit = 'cm'), 
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8))
  # combine plots
  plot_final <- plot_hic + plot_edge 
  return(plot_final)
}

## process
bin_anno <- load_bin_anno(bin_file, cell_label)
spin_res = 25000
pdf(file = pdf_file, width = pdf_width, height = pdf_height)
print(plot_all_data(bin_anno, hic_file, edge_file, hic_merge_factor, spin_res, chrom_1, start_1, end_1, chrom_2, start_2, end_2))
dev.off()
