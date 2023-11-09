library(augsynth)
library(data.table)
library(dplyr)
library(fixest)
library(ggplot2)


# Getting the workign directory 
shell_root <- "/scratch/gpfs/ar8787/groupdata2/india_forest_land" 
dbox_root <- "~/Dropbox/india_forest_land" 
root <- shell_root
setwd( root )

### Geting State names
state_names <- fread( "A_MicroData/state_names_shrug.csv" )

# Getting files
main_path <- "E_Estimates/_vcf_10_augsyth/state"
all_sum <- list.files( main_path, pattern = "*_sum*")

# Getting the num of files
df1 <- fread("A_MicroData/data_sysdif.csv")



main_path <- "/scratch/gpfs/ar8787/groupdata2/india_forest_land/E_Estimates/_vcf_10_augsyth/state"
all_sum <- list.files( main_path, pattern = "*_sum*")

for (name in all_sum){
  
  state_id_code <- gsub(".*st_(.*?)_sum.*", "\\1", name )  
  st_name <- state_names[ state_names$state_id == state_id_code ]$state_name
  ppool_syn_time_summ <- readRDS( file.path(main_path, name ))
  
  
  # Get the plot
  est_df <- ppool_syn_time_summ$att
  est_df <- est_df[est_df$Level == "Average", ]
  est_df <- est_df[apply(!is.na(est_df), 1, all),]
  
  df2 <- df1[df1$pc11_state_id  == state_id_code, ]
  n_df <- df2[t2ev_ror_data_entry %in% est_df$Time][, .N, by = t2ev_ror_data_entry]
  n_df$Time <- n_df$t2ev_ror_data_entry
  setorder(n_df,Time )
  setorder(est_df,Time )
  
  
  
  # Getting factors
  rect.length <- (max(est_df[,"Estimate"], na.rm = TRUE) - 
                    min(est_df[,"Estimate"], na.rm = TRUE))/2
  scale_fac <- 0.8 * rect.length / ( max(n_df[,"N"]) )
  min_y_lim <- round(min(est_df$lower_bound), 2) * 1.05
  max_y_lim <- round( max(est_df$upper_bound), 2) * 0.95
  
  est_df[,"xmin"] <- est_df[,"Time"] - 0.2
  est_df[,"xmax"] <- est_df[,"Time"] + 0.2
  est_df[,"ymin"] <- min_y_lim
  est_df[,"ymax"] <- est_df[,"ymin"] + ( n_df[,"N"] * scale_fac )
  
  
  # Getting the plot
  p <- ggplot(est_df, aes(x = Time, y = Estimate)) + 
    geom_line(color = "black") +
    geom_point(color = "blue") +
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), 
                  width = 0.2, color = "red") +
    geom_rect(data = est_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
              fill = "grey70", colour = "grey69", alpha = 0.4, size = 0.2) + 
    coord_cartesian(ylim = c( min_y_lim, max_y_lim ) ) +
    scale_y_continuous(name = "Estimate", 
                       sec.axis = sec_axis(~(.+ (-1*min_y_lim)) * (1/scale_fac), 
                                           name = "Number of Observations")) + 
    labs(y = "ATT", x = "Time to Event", title = "Ror Data Entry at Block Level if more than 50% Villages are Treated") + 
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +  # Dashed vertical line at x = 5
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "gray") + 
    ggtitle( st_name ) + 
    theme_minimal()
  
  filename <- paste0( "espec_st_" , state_id_code, "_plot.png")
  path <- file.path("F_Figures/_vcf_10_augsyth/state", filename)
  ggsave(filename = path, plot = p, width = 6, height = 4, dpi = 300)
  
  
  
}