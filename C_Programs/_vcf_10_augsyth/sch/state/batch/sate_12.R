

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


# Getting the num of files
df1 <- fread("A_MicroData/data_sysdif.csv")


# Use tstrsplit to split the 'shrid' column and assign the output to new columns:
df1[, c("pc11_s_id", "pc11_d_id", "pc11_sd_id", "pc11_tv_id") 
    := tstrsplit(shrid, "-", fixed=TRUE)]
cols_to_convert <- c("pc11_s_id", "pc11_d_id", "pc11_sd_id", "pc11_tv_id")
df1[, (cols_to_convert) := lapply(.SD, as.numeric), .SDcols = cols_to_convert]

# Import village
vill_sch <- fread("A_MicroData/vill_sch.csv")
cols_to_convert <- c("pc11_s_id", "pc11_d_id", "pc11_sd_id", "pc11_tv_id")
vill_sch[, (cols_to_convert) := lapply(.SD, as.numeric), .SDcols = cols_to_convert]


df2_inner <- df1 %>% inner_join( vill_sch, 
                           by = c("pc11_s_id", "pc11_d_id", "pc11_sd_id", "pc11_tv_id"))
list_ids = unique(df2_inner$pc11_state_id)
### Geting State names
state_names <- fread( "A_MicroData/state_names_shrug.csv" )
state_id_code <- list_ids[11]
st_name <- state_names[ state_names$state_id == state_id_code ]$state_name


df2 = df2_inner[df2_inner$pc11_state_id == state_id_code]
df2$vill_id <- as.numeric( df2$vill_id )

df2_sch <- df2[df2$sch == 1]
df2_nonsch <- df2[df2$sch == 0]

n_sch = dim(df2_sch)[1]
n_nonsch = dim(df2_sch)[1]

#-------------------------------------------------------------------------------
# sch the median
#-------------------------------------------------------------------------------



result <- tryCatch({
  
  # Estimatio Result
  ppool_syn <- multisynth(per_treecover ~ post_ror_data_entry, 
                          vill_id, year, df2_sch )
  
  filename <- paste0( "sch", state_id_code, ".RDS")
  path <- file.path("E_Estimates/_vcf_10_augsyth/sch/by_state", filename)
  saveRDS( ppool_syn, path )
  
  
  # Get Summary result
  ppool_syn_time_summ <- summary( ppool_syn )
  filename <- paste0( "sch", state_id_code, "_sum.RDS")
  path <- file.path("E_Estimates/_vcf_10_augsyth/sch/by_state", filename)
  saveRDS( ppool_syn_time_summ,  path )
  
  
  
  
  # Get the plot
  est_df <- ppool_syn_time_summ$att
  est_df <- est_df[est_df$Level == "Average", ]
  est_df <- est_df[apply(!is.na(est_df), 1, all),]
  
  
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
  
  
  filename <- paste0( "sch_", state_id_code, "_plot.png")
  path <- file.path("F_Figures/_vcf_10_augsyth/sch/by_state", filename)
  ggsave(filename = path, plot = p, width = 6, height = 4, dpi = 300)
  
}, error = function(err) {
  
  # Code to handle the error
  # You can access information about the error in 'err'
  cat("An error occurred:", conditionMessage(err), "
")
  
  return( 'NA')
})
print(result)





#-------------------------------------------------------------------------------
# nonsch the median
#-------------------------------------------------------------------------------



result <- tryCatch({
  
  
  # Estimatio Result
  ppool_syn <- multisynth(per_treecover ~ post_ror_data_entry, 
                          vill_id, year, df2_nonsch, time_cohort = TRUE, n_leads = 9 )
  
  filename <- paste0( "nonsch", state_id_code, ".RDS")
  path <- file.path("E_Estimates/_vcf_10_augsyth/sch/by_state", filename)
  saveRDS( ppool_syn, path )
  
  
  # Get Summary result
  ppool_syn_time_summ <- summary( ppool_syn )
  filename <- paste0( "nonsch", state_id_code, "_sum.RDS")
  path <- file.path("E_Estimates/_vcf_10_augsyth/sch/by_state", filename)
  saveRDS( ppool_syn_time_summ,  path )
  
  
  
  
  # Get the plot
  est_df <- ppool_syn_time_summ$att
  est_df <- est_df[est_df$Level == "Average", ]
  est_df <- est_df[apply(!is.na(est_df), 1, all),]
  
  
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
  
  
  filename <- paste0( "nonsch_", state_id_code, "_plot.png")
  path <- file.path("F_Figures/_vcf_10_augsyth/sch/by_state", filename)
  ggsave(filename = path, plot = p, width = 6, height = 4, dpi = 300)
  
  
  
}, error = function(err) {
  # Code to handle the error
  # You can access information about the error in 'err'
  cat("An error occurred:", conditionMessage(err), "
")
  
  # Optionally, return a default or custom value
  return( "NA")
})


print(result)


