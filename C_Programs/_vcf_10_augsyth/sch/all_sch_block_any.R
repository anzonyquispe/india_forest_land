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


# Import village
vill_sch <- fread("A_MicroData/vill_sch.csv")
cols_to_convert <- c("pc11_s_id", "pc11_d_id", "pc11_sd_id", "pc11_tv_id")
vill_sch[, (cols_to_convert) := lapply(.SD, as.numeric), .SDcols = cols_to_convert]

# Calculate the mean of 'value' by 'group'
block_sch <- vill_sch[, .(sch = max(sch)), by = .(pc11_s_id, pc11_d_id, pc11_sd_id)]
# Rename columns
setnames(block_sch, old = c('pc11_s_id', 'pc11_d_id', 'pc11_sd_id') , 
         new = c("state_id", "district_id", "subdist_id" ))


# Import data
df1 <- fread( "A_MicroData/data_sysdif_block_lvl.csv" )
cols_to_convert <- c("state_id", "district_id", "subdist_id" )
df1[, (cols_to_convert) := lapply(.SD, as.numeric), .SDcols = cols_to_convert]


# Generation of t2event
# Create a data.table with minimum year by block_id
min_years_data <- df1[
  post_ror_data_entry_block_any == 1,
  .(min_year = min(year)), # Replace 'year_column' with your actual year column name
  by = .(block_id)
]

# Merge the min_years_data back to the original dataset
df1 <- merge(df1, min_years_data, by = "block_id", all.x = TRUE)
df1$t2ev_ror_data_entry_block_any <- df1$year - df1$min_year


df2 <- df1 %>% inner_join( block_sch, 
                           by = c("state_id", "district_id", "subdist_id" ) )





df2_sch <- df2[df2$sch == 1]
df2_nonsch <- df2[df2$sch == 0]




#-------------------------------------------------------------------------------
# sch the median
#-------------------------------------------------------------------------------

# Estimatio Result
ppool_syn <- multisynth(per_treecover ~ post_ror_data_entry_block_any, 
                        block_id, year, df2_sch, time_cohort = TRUE, n_leads = 9 )

filename <- paste0( "sch_block_any.RDS")
path <- file.path("E_Estimates/_vcf_10_augsyth/sch", filename)
saveRDS( ppool_syn, path )


# Get Summary result
ppool_syn_time_summ <- summary( ppool_syn )
filename <- paste0( "sch_sum_block_any.RDS")
path <- file.path("E_Estimates/_vcf_10_augsyth/sch", filename)
saveRDS( ppool_syn_time_summ,  path )

# Import File
filename <- paste0( "sch_sum_block_any.RDS")
path <- file.path("E_Estimates/_vcf_10_augsyth/sch", filename)
ppool_syn_time_summ <- readRDS( path )

# Get the plot
est_df <- ppool_syn_time_summ$att
est_df <- est_df[est_df$Level == "Average", ]
est_df <- est_df[apply(!is.na(est_df), 1, all),]


n_df <- df2_sch[t2ev_ror_data_entry_block_any %in% est_df$Time][, .N, by = t2ev_ror_data_entry_block_any]
n_df$Time <- n_df$t2ev_ror_data_entry_block_any
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
  labs(y = "ATT", x = "Time to Event", title = "Sch" ) + 
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +  # Dashed vertical line at x = 5
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "gray") + 
  theme_minimal()

filename <- paste0( "sch_block_any_plot.png")
path <- file.path("F_Figures/_vcf_10_augsyth/sch", filename)
ggsave(filename = path, plot = p, width = 6, height = 4, dpi = 300)




#-------------------------------------------------------------------------------
# nonsch the median
#-------------------------------------------------------------------------------


# Estimatio Result
ppool_syn <- multisynth(per_treecover ~ post_ror_data_entry_block_any, 
                        block_id, year, df2_nonsch, time_cohort = TRUE, n_leads = 9 )

filename <- paste0( "nonsch_all_block_any.RDS")
path <- file.path("E_Estimates/_vcf_10_augsyth/sch", filename)
saveRDS( ppool_syn, path )


# Get Summary result
ppool_syn_time_summ <- summary( ppool_syn )
filename <- paste0( "nonsch_all_sum_block_any.RDS")
path <- file.path("E_Estimates/_vcf_10_augsyth/sch", filename)
saveRDS( ppool_syn_time_summ,  path )


# Import File
filename <- paste0( "nonsch_all_sum_block_any.RDS")
path <- file.path("E_Estimates/_vcf_10_augsyth/sch", filename)
ppool_syn_time_summ <- readRDS( path )


# Get the plot
est_df <- ppool_syn_time_summ$att
est_df <- est_df[est_df$Level == "Average", ]
est_df <- est_df[apply(!is.na(est_df), 1, all),]


n_df <- df2_nonsch[t2ev_ror_data_entry_block_any %in% est_df$Time][, .N, by = t2ev_ror_data_entry_block_any]
n_df$Time <- n_df$t2ev_ror_data_entry_block_any
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
  labs(y = "ATT", x = "Time to Event", title = "Nonsch" ) + 
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +  # Dashed vertical line at x = 5
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "gray") + 
  theme_minimal()

filename <- paste0( "nonsch_all_block_any_plot.png")
path <- file.path("F_Figures/_vcf_10_augsyth/sch", filename)
ggsave(filename = path, plot = p, width = 6, height = 4, dpi = 300)

