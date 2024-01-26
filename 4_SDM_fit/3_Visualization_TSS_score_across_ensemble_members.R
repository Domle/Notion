## Boxplot with species on the x-axis and the y-axis shows TSS score. Each boxplot summarizes the TSS score across the 5 predictor ensembles.

# Input files: 1. CSV files containing information on species and tss scores
# Output files: 1. A png plot with species on x-axis and y-axis contains the tss scores

## Author:  Dominic Eriksson
##          Environmental Physics Group, UP
##          ETHZ, Zürich
##          Switzerland

# deriksson@ethz.ch, 19 of September 2023 ------------------------------------------------------------------------------------




# Functions
source("/home/deriksson/Projects/Notion_DE/Code/Functions/Dominic/call_libraries.R")

# Libraries
lib_vec <- c("ggplot2", "patchwork")
call_libraries(lib_vec)

# Directories
wd_eval <- "/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/"
wd_out <- "/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/4_Output/"

# Vector of filenames
vec_folders <- c("Total_dataset", "MicroscopyBased_dataset", "SequenceBased_dataset")
f <- 1
#
vec_fname <- list.files(paste0(wd_eval, vec_folders[f]))
vec_fname_cr_bg_nonov <- grep("cr_bg_nonov", vec_fname, value = TRUE)
vec_fname_cr_bg_overl <- grep("cr_bg_overl", vec_fname, value = TRUE)
vec_fname_gr_bg_nonov <- grep("gr_bg_nonov", vec_fname, value = TRUE)
vec_fname_gr_bg_overl <- grep("gr_bg_overl", vec_fname, value = TRUE)
vec_fname_tot_bg_nonov <- grep("tot_bg_nonov", vec_fname, value = TRUE)
vec_fname_tot_bg_overl <- grep("tot_bg_overl", vec_fname, value = TRUE)
# Save in
l_vec_fname <- list()
l_vec_fname[[1]] <- vec_fname_tot_bg_overl
l_vec_fname[[2]] <- vec_fname_tot_bg_nonov
l_vec_fname[[3]] <- vec_fname_gr_bg_overl
l_vec_fname[[4]] <- vec_fname_gr_bg_nonov
l_vec_fname[[5]] <- vec_fname_cr_bg_overl
l_vec_fname[[6]] <- vec_fname_cr_bg_nonov

## Format data
l_d <- list()
for(d in seq_along(l_vec_fname)){
    print(d)
    df_gam <- read.csv(paste0(wd_eval, vec_folders[f], "/", l_vec_fname[[d]][1]))
    df_glm <- read.csv(paste0(wd_eval, vec_folders[f], "/", l_vec_fname[[d]][2]))
    df_rf <- read.csv(paste0(wd_eval, vec_folders[f], "/", l_vec_fname[[d]][3]))
    #
    df_gam <- df_gam[, c("taxon", "tss.xval4")]
    df_glm <- df_glm[, c("taxon", "tss.xval4")]
    df_rf <- df_rf[, c("taxon", "tss.xval4")]
    #
    coln_gam <- gsub("\\.csv", "", l_vec_fname[[d]][1])
    coln_glm <- gsub("\\.csv", "", l_vec_fname[[d]][2])
    coln_rf <- gsub("\\.csv", "", l_vec_fname[[d]][3])
    #
    df_gam$pseudoAbs <- coln_gam
    df_glm$pseudoAbs <- coln_glm
    df_rf$pseudoAbs <- coln_rf
    #
    df <- rbind(df_glm, df_gam, df_rf)
    l_d[[d]] <- df
}

# Loop
l_p <- list()
for(v in seq_along(l_d)){

    # Print progress
    print(v)

    # Get data
    d <- l_d[[v]]
    d <- d[which(!is.na(d$tss.xval4)), ]

    # Plot
    gg <- ggplot(d, aes(x = taxon, y = tss.xval4)) + 
        geom_boxplot(show.legend = FALSE, width = .5, lwd = 0.1, outlier.size = 0.1) + xlab("") + ylab("") + 
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2.5 ),
            axis.text.y = element_text(size = 3 ),
            axis.title.y = element_text(size = 5)) + 
            facet_wrap(~pseudoAbs, ncol = 3, scales ="free") +
            theme(strip.text.x = element_text(size = 3)) + ylab("TSS score")
            # ggtitle(vec_fname[v]) 

    l_p[[v]] <- gg

    # Store plot     
    # ggsave(paste0(wd_out, f, ".png"), plot = gg)

} # end of loop over vec_fname

# Plot with patchwork
pp <- l_p[[1]] +
l_p[[2]] +
l_p[[3]] + 
l_p[[4]] + 
l_p[[5]] + 
l_p[[6]] + 
plot_layout(nrow = 3)


# Save plot
# ggsave(paste0(wd_out, vec_folders[f], "/TSS_scores_across_ensemble_members.pdf"), plot = pp, dpi = 300)
ggsave(paste0(wd_out, vec_folders[f], "/TSS_scores_across_ensemble_members.png"), plot = pp, dpi = 300)
