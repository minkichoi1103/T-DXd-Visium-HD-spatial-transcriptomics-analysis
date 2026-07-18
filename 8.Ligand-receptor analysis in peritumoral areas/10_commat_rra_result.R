
# ================= 0. libraries ============================================
suppressPackageStartupMessages({
  library(RobustRankAggreg)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

plot_slope_final <- function(df_all, sel_names,
                                top_n           = 5,          
                                label_indices   = NULL,       
                                title_txt       = "",
                                palette         = NULL,
                                add_labels      = FALSE,
                                base_font_size  = 14,
                                title_font_size = 18,
                                label_font_size = 4
) {
  
  highlight_names <- head(sel_names, top_n)
  
  
  if (!is.null(label_indices)) {
    valid_indices <- label_indices[label_indices <= length(sel_names)]
    label_names <- sel_names[valid_indices]
  } else {
    label_names <- highlight_names
  }
  
  long <- df_all %>%
    tidyr::pivot_longer(cols = c(Rank_sens, Rank_res),
                        names_to = "Cond", values_to = "Rank") %>%
    dplyr::mutate(
      Cond  = dplyr::recode(Cond,
                            "Rank_sens" = "Sensitive",
                            "Rank_res"  = "Resistant"),
      Cond  = factor(Cond, levels = c("Sensitive","Resistant")),

      in_sel = Name %in% highlight_names
    )
  
  bg_df  <- subset(long, !in_sel) # for the gray background
  sel_df <- subset(long,  in_sel) # for the colored highlight
  
  if (is.null(palette)) palette <- scales::hue_pal()(length(highlight_names))
  if(length(palette) < length(highlight_names)) {
    palette <- rep(palette, length.out = length(highlight_names))
  }
  
  # Draw the plot
  p <- ggplot() +
    # 1) Background (gray)
    geom_line(data = bg_df, aes(Cond, Rank, group = Name),
              color = "grey80", linewidth = 0.5) +
    geom_point(data = bg_df, aes(Cond, Rank, group = Name),
               color = "grey80", size = 1) +
    
    # 2) Highlight (color) - plot all data in top_n
    geom_line(data = sel_df,
              aes(Cond, Rank, group = Name, color = Name),
              linewidth = 0.9) +
    geom_point(data = sel_df,
               aes(Cond, Rank, color = Name),
               size = 2) +
    
    scale_color_manual(values = palette, name = "LR interaction") +
    scale_y_reverse() +
    labs(title = title_txt, x = 'Sensitivity of T-Dxd', y = "Rank") +
    
    theme_minimal(base_size = base_font_size) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = title_font_size, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = base_font_size + 2, face = "plain"),
      axis.text = element_text(size = base_font_size),
      legend.text = element_text(size = base_font_size - 2)
    )
  
  if (add_labels && length(label_names) > 0) {

    label_df <- subset(sel_df, (Cond == "Resistant") & (Name %in% label_names))
    
    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(label = Name, color = Name, x = Cond, y = Rank),
        max.overlaps  = Inf,
        force         = 3.0,
        box.padding   = 1.0,
        point.padding = 0.6,
        nudge_x       = 0.25,
        size          = label_font_size, 
        show.legend   = FALSE
      )
  }
  
  return(p)
}


# ================= 1. first-stage RRA (batch → sample_id) ==========
setwd("../data/")

## 1. Load batch–level ranked lists -----------------------------------------
df_batches <- read.csv("commot_ranked_result.csv",
                       row.names = 1, check.names = FALSE)

result_list <- lapply(seq_len(nrow(df_batches)), function(i) {
  v <- unlist(df_batches[i, ], use.names = FALSE)
  v[v != ""]
})
names(result_list) <- rownames(df_batches)

## 2. Load metadata ----------------------------------------------------------
meta <- read.csv(
  "../data/results_df_spatial_vessel_linkerenz.csv",
  row.names   = 1,
  check.names = FALSE
)
stopifnot(all(c("batch", "sample_id", "clinical_type", "HER2_type")
              %in% colnames(meta)))

## 3. First-stage RRA (aggregate per sample_id) ------------------------------
vec_from_rra <- function(agg_df) agg_df$Name[order(agg_df$Score)]

rank_per_list <- lapply(split(meta, meta$sample_id), function(sub) {
  lst <- result_list[sub$batch]
  lst <- lst[!sapply(lst, is.null)]
  if (length(lst) == 0) return(NULL)
  if (length(lst) == 1) return(lst[[1]])
  vec_from_rra(aggregateRanks(lst))
})
rank_per_list <- rank_per_list[!sapply(rank_per_list, is.null)]


# ================= 2. union-mode comparison ===============================
compare_rra <- function(ids_res, ## condition of interest
                        ids_sens,  
                        rank_list,
                        ties_method = "average",
                        cutoff      = 0.05) {
  # ----- RRA -----------------------------------------------------
  sens_rra <- aggregateRanks(rank_list[ids_sens])
  res_rra  <- aggregateRanks(rank_list[ids_res])
  
  sens_rra <- sens_rra %>% mutate(Rank_sens = rank(Score, ties.method = ties_method))
  res_rra  <- res_rra  %>% mutate(Rank_res  = rank(Score, ties.method = ties_method))
  n_sens <- nrow(sens_rra); n_res <- nrow(res_rra)
  
  # ----- union & NA-score handling ----------------------------------------
  df_uni <- full_join(
    sens_rra %>% select(Name, Score_sens = Score, Rank_sens),
    res_rra  %>% select(Name, Score_res  = Score, Rank_res ),
    by = "Name"
  ) %>%
    # fill missing scores with 1.5 (non-significant high value)
    mutate(
      Score_sens = ifelse(is.na(Score_sens), 1.5, Score_sens),
      Score_res  = ifelse(is.na(Score_res ), 1.5, Score_res ),
      Rank_sens  = ifelse(is.na(Rank_sens),  n_sens + 1, Rank_sens),
      Rank_res   = ifelse(is.na(Rank_res ),  n_res  + 1, Rank_res ),
      Delta      = Rank_sens - Rank_res           # + : better in Resistant
    )
  
  # ----- Resistant-specific set (UP) --------------------------------------
  up_df <- df_uni %>% 
    filter(Score_res < cutoff,           # significant only in Res
           Score_sens >= cutoff,
           Delta > 0) %>%               # rank improved
    arrange(desc(Delta))
  
  list(
    res_specific_df  = up_df,
    full_union       = df_uni
  )
}



# ================= 3. Comparison between Sensitive and Resistant ===========================================
sample_info <- meta %>% distinct(sample_id, clinical_type, HER2_type)

# HER2 positive
ids_sens_pos <- sample_info$sample_id[
  sample_info$clinical_type == "Sensitive" & sample_info$HER2_type == "Pos"]
ids_res_pos  <- sample_info$sample_id[
  sample_info$clinical_type == "Resistant" & sample_info$HER2_type == "Pos"]

out_pos <- compare_rra(ids_res_pos, ids_sens_pos, rank_per_list,
                       cutoff = 0.05)

# HER2 low
ids_sens_low <- sample_info$sample_id[
  sample_info$clinical_type == "Sensitive" & sample_info$HER2_type == "Low"]
ids_res_low  <- sample_info$sample_id[
  sample_info$clinical_type == "Resistant" & sample_info$HER2_type == "Low"]

out_low <- compare_rra(ids_res_low, ids_sens_low, rank_per_list,
                       cutoff = 0.05)

################### Result plot ################### 
##### HER2 positive ######
# plot of HER2 positive 

set.seed(400)
pos_highlight_names <- head(out_pos$res_specific_df$Name, 8)
palette <- sample(ggsci::pal_npg("nrc")(8)) 
names(palette) <- pos_highlight_names # assign colors in the order of the LR names


plt_up_positive <- plot_slope_final(
  out_pos$full_union, out_pos$res_specific_df$Name,
  top_n     = 8,
  label_indices = c(1, 7, 8),
  title_txt = paste0(""),
  palette   = palette,   
  add_labels = TRUE,
  base_font_size  = 15, 
  title_font_size = 15, 
  label_font_size = 4   
)

print(plt_up_positive)
ggsave(filename = "../figures/HER2_pos_res_lr_plot_pos.png", 
       plot = plt_up_positive, 
       dpi = 400, 
       width = 8, height = 6, units = "in", bg = "white") 


# Compare in HER2 low
compare_lr <- c("AREG-EGFR", "PSAP", "PSAP-GPR37")
plt_up_negative <- plot_slope_final(
  out_low$full_union, compare_lr,
  top_n     = length(compare_lr),
  title_txt = paste0(""),
  #palette   = scales::hue_pal()(length(compare_lr)), 
  palette   = palette,
  add_labels = TRUE,
  base_font_size  = 15, 
  title_font_size = 15, 
  label_font_size = 4 
)

print(plt_up_negative)
ggsave(filename = "../figures/HER2_pos_res_lr_plot_low.png", 
       plot = plt_up_negative, 
       dpi = 400, 
       width = 8, height = 6, units = "in", bg = "white") 






##### HER2 low ######
out_low_backup <- out_low
out_low$full_union <- out_low$full_union[diff != -304, ]

# plot of HER2 low (selected)
set.seed(123)
# 1. Extract colors from each high-quality palette
pal_npg <- ggsci::pal_npg("nrc")(10)      # Nature (10 colors)
pal_lancet <- ggsci::pal_lancet("lanonc")(9) # Lancet (9 colors)
pal_jama <- ggsci::pal_jama("default")(6)   # JAMA (8 colors)

# 2. Combine all colors (25 total)
palette_low <- sample(c(pal_npg, pal_lancet, pal_jama))

low_highlight_names <- out_low$res_specific_df$Name
names(palette_low) <- low_highlight_names # assign colors in the order of the LR names


#1. rank1-13
compare_lr <- out_low$res_specific_df$Name[1:13]
plt_low_total_1 <- plot_slope_final(
  out_low$full_union, compare_lr,
  top_n     = length(compare_lr),
  title_txt = paste0(""),
  palette   = palette_low[1:13],
  add_labels = FALSE,
  base_font_size  = 15, 
  title_font_size = 15, 
  label_font_size = 4   
)

print(plt_low_total_1)
ggsave(filename = "../figures/HER2_low_res_lr_plot_1.png", 
       plot = plt_low_total_1, 
       dpi = 400, 
       width = 8, height = 6, units = "in", bg = "white") 

#1. rank14-25
compare_lr <- out_low$res_specific_df$Name[14:25]
plt_low_total_2 <- plot_slope_final(
  out_low$full_union, compare_lr,
  top_n     = length(compare_lr),
  title_txt = paste0(""),
  palette   = palette_low[14:25],
  add_labels = FALSE,
  base_font_size  = 15, 
  title_font_size = 15, 
  label_font_size = 4   
)

print(plt_low_total_2)
ggsave(filename = "../figures/HER2_low_res_lr_plot_2.png", 
       plot = plt_low_total_2, 
       dpi = 400, 
       width = 8, height = 6, units = "in", bg = "white") 



# Compare in HER2 positive : HER2 low resistant specific, VEGF related
compare_lr <- c("VEGFB-FLT1", "VEGFA-KDR", "VEGFA-FLT1")
plt_pos_vegf <- plot_slope_final(
  out_pos$full_union, compare_lr,
  top_n     = length(compare_lr),
  title_txt = paste0(""),
  palette   = palette_low,
  add_labels = TRUE,
  base_font_size  = 15, 
  title_font_size = 15, 
  label_font_size = 4.5   
)
print(plt_pos_vegf)
ggsave(filename = "../figures/HER2_pos_lr_plot_vegf.png", 
       plot = plt_pos_vegf, 
       dpi = 400, 
       width = 8, height = 6, units = "in", bg = "white") 


compare_lr <- c("VEGFB-FLT1", "VEGFA-KDR", "VEGFA-FLT1")
plt_low_vegf <- plot_slope_final(
  out_low$full_union, compare_lr,
  top_n     = length(compare_lr),
  title_txt = paste0(""),
  palette   = palette_low,
  add_labels = TRUE,
  base_font_size  = 15, 
  title_font_size = 15, 
  label_font_size = 4.5   
)

print(plt_low_vegf)
ggsave(filename = "../figures/HER2_low_lr_plot_vegf.png", 
       plot = plt_low_vegf, 
       dpi = 400, 
       width = 8, height = 6, units = "in", bg = "white") 


# Compare in HER2 positive : HER2 low resistant specific, PDGFR related
compare_lr <- c("PDGFC-PDGFRA", "PDGFB-PDGFRA", "PDGFA-PDGFRB")
plt_pos_pdgf <- plot_slope_final(
  out_pos$full_union, compare_lr,
  top_n     = length(compare_lr),
  title_txt = paste0(""),
  palette   = palette_low,
  add_labels = TRUE,
  base_font_size  = 15, 
  title_font_size = 15, 
  label_font_size = 4.5   
)
print(plt_pos_pdgf)  # HER2 pos
ggsave(filename = "../figures/HER2_pos_lr_plot_pdgf.png", 
       plot = plt_pos_pdgf, 
       dpi = 400, 
       width = 8, height = 6, units = "in", bg = "white") 


compare_lr <- c("PDGFC-PDGFRA", "PDGFB-PDGFRA", "PDGFA-PDGFRB")
plt_low_pdgf <- plot_slope_final(
  out_low$full_union, compare_lr,
  top_n     = length(compare_lr),
  title_txt = paste0(""),
  palette   = palette_low,
  add_labels = TRUE,
  base_font_size  = 15, 
  title_font_size = 15,
  label_font_size = 4.5   
)
print(plt_low_pdgf)# HER2 low
ggsave(filename = "../figures/HER2_low_lr_plot_pdgf.png", 
       plot = plt_low_pdgf, 
       dpi = 400, 
       width = 8, height = 6, units = "in", bg = "white") 
