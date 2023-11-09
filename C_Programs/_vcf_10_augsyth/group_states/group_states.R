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

### Geting State names
state_names <- fread( "A_MicroData/state_names_shrug.csv" )


# Union territories except delhi
# Ladakh
# Jammu & Kashmir
# Puducherry
# Lakshadweep
# Chandigarh
# Dadra and Nagar Haveli and Daman & Diu
# Andaman and Nicobar Islands

unioncat <- c(1, 34, 31, 4, 26, 35 )
state_names %>% 
  filter( state_id %in% unioncat )

df2_union <- df1 %>% 
  filter(pc11_state_id %in% unioncat )


# "northeastern states (except Assam)
# Arunachal Pradesh, Manipur, Meghalaya, Mizoram, Nagaland, Sikkim, Tripura.

north_st <- c(12, 14, 17, 15, 13, 11, 16)
state_names %>% 
  filter( state_id %in% north_st )


df2_north <- df1 %>% 
  filter(pc11_state_id %in% north_st )

# rest of states
`%nin%` = Negate(`%in%`)
df2_rest <- df1 %>% 
  filter(pc11_state_id %nin% unioncat ) %>% 
  filter(pc11_state_id %nin% north_st )

#-------------------------------------------------------------------------------
# Union territories except delhi
#-------------------------------------------------------------------------------

result <- tryCatch({
  
  # Estimatio Result
  ppool_syn <- multisynth(per_treecover ~ post_ror_data_entry, 
                          vill_id, year, df2_union )
  
  filename <- paste0( "union_states.RDS")
  path <- file.path("E_Estimates/_vcf_10_augsyth/group_states", filename)
  saveRDS( ppool_syn, path )
  
  
  # Get Summary result
  ppool_syn_time_summ <- summary( ppool_syn )
  filename <- paste0( "union_states_sum.RDS" )
  path <- file.path("E_Estimates/_vcf_10_augsyth/group_states", filename)
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
    labs(y = "ATT", x = "Time to Event", title = "Ror Data Entry at Village Level") + 
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +  # Dashed vertical line at x = 5
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "gray") + 
    ggtitle( "Union Territories (Except Delhi)" ) +
    theme_minimal()
  
  
  filename <- paste0( "union_states_plot.png")
  path <- file.path("F_Figures/_vcf_10_augsyth/group_states", filename)
  ggsave(filename = path, plot = p, width = 6, height = 4, dpi = 300)

}, error = function(err) {
  
  # Code to handle the error
  # You can access information about the error in 'err'
  cat("An error occurred:", conditionMessage(err), "\n")
  
  return( 'NA')
})
print(result)





#-------------------------------------------------------------------------------
# north the median
#-------------------------------------------------------------------------------



result <- tryCatch({
  
  
  # Estimatio Result
  ppool_syn <- multisynth(per_treecover ~ post_ror_data_entry, 
                          vill_id, year, df2_north, time_cohort = TRUE, n_leads = 9 )
  
  filename <- paste0( "north_states.RDS")
  path <- file.path("E_Estimates/_vcf_10_augsyth/group_states", filename)
  saveRDS( ppool_syn, path )
  
  
  # Get Summary result
  ppool_syn_time_summ <- summary( ppool_syn )
  filename <- paste0( "north_states_sum.RDS")
  path <- file.path("E_Estimates/_vcf_10_augsyth/group_states", filename)
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
    labs(y = "ATT", x = "Time to Event", title = "Ror Data Entry at Village Level") + 
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +  # Dashed vertical line at x = 5
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "gray") + 
    ggtitle( "Northeastern States (Except Assam)" ) +
    theme_minimal()
  
  
  filename <- paste0( "north_states_plot.png")
  path <- file.path("F_Figures/_vcf_10_augsyth/group_states", filename)
  ggsave(filename = path, plot = p, width = 6, height = 4, dpi = 300)
  
  
  
}, error = function(err) {
  # Code to handle the error
  # You can access information about the error in 'err'
  cat("An error occurred:", conditionMessage(err), "\n")
  
  # Optionally, return a default or custom value
  return( "NA")
})


print(result)




#-------------------------------------------------------------------------------
#  The rest of the States
#-------------------------------------------------------------------------------

result <- tryCatch({
  
  
  # Estimatio Result
  ppool_syn <- multisynth(per_treecover ~ post_ror_data_entry, 
                          vill_id, year, df2_rest, time_cohort = TRUE, n_leads = 9 )
  
  filename <- paste0( "rest_states.RDS")
  path <- file.path("E_Estimates/_vcf_10_augsyth/group_states", filename)
  saveRDS( ppool_syn, path )
  
  
  # Get Summary result
  ppool_syn_time_summ <- summary( ppool_syn )
  filename <- paste0( "rest_states_sum.RDS")
  path <- file.path("E_Estimates/_vcf_10_augsyth/group_states", filename)
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
    labs(y = "ATT", x = "Time to Event", title = "Ror Data Entry at Village Level") + 
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +  # Dashed vertical line at x = 5
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "gray") + 
    ggtitle( "Rest of States" ) +
    theme_minimal()
  
  
  filename <- paste0( "rest_states_plot.png")
  path <- file.path("F_Figures/_vcf_10_augsyth/group_states", filename)
  ggsave(filename = path, plot = p, width = 6, height = 4, dpi = 300)
  
  
  
}, error = function(err) {
  # Code to handle the error
  # You can access information about the error in 'err'
  cat("An error occurred:", conditionMessage(err), "\n")
  
  # Optionally, return a default or custom value
  return( "NA")
})


print(result)

