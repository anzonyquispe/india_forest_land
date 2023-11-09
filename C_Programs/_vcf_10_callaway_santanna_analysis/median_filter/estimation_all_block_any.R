library(augsynth)
library(data.table)
library(dplyr)
library(fixest)
library(ggplot2)
library(did)

# Getting the workign directory 
shell_root <- "/scratch/gpfs/ar8787/groupdata2/india_forest_land" 
dbox_root <- "~/Dropbox/india_forest_land" 
root <- shell_root
setwd( root )


# Getting the num of files
# Import data
df1 <- fread( "A_MicroData/data_sysdif_block_lvl.csv" )
names(df1)

df_year_2010 <- df1[df1$year == 2010]


# Calculate the median of the 'value' column
median_value <- median(df_year_2010$per_treecover)
# Generate a new variable 'above_median'; it will be TRUE if 'value' is greater than the median
df_year_2010[, above_median_all := (per_treecover > median_value)*1 ]
df_year_2010 <- df_year_2010[, .(block_id, above_median_all )]

df2 <- df1 %>% left_join( df_year_2010 )

# Generation of t2event
# Create a data.table with minimum year by block_id
min_years_data <- df2[
  post_ror_data_entry_block_any == 1,
  .(min_year = min(year)), # Replace 'year_column' with your actual year column name
  by = .(block_id)
]

# Merge the min_years_data back to the original dataset
df2 <- merge(df2, min_years_data, by = "block_id", all.x = TRUE)
df2$t2ev_ror_data_entry_block_any <- df2$year - df2$min_year



df2_avobe <- df2[df2$above_median_all == 1]
df2_below <- df2[df2$above_median_all == 0]



#-------------------------------------------------------------------------------
# Below the median
#-------------------------------------------------------------------------------

# Start estimation
out <- att_gt(yname = "per_treecover",
              gname = "min_year",
              idname = "block_id",
              tname = "year",
              xformla = ~ 1 ,
              data = df2_below,
              est_method = "reg", 
              control_group = c( "notyettreated")
)

filename <- paste0( "below_median_block_any.RDS")
path <- file.path("E_Estimates/_vcf_10_callaway_santanna_analysis/median_filter", 
                  filename)
saveRDS( out, path )



# aggregate the group-time average treatment effects
dynamic_sum <- aggte( out, type = "dynamic")
res <- ggdid(dynamic_sum)
data_sum <- res$data
data_sum$ymin=(data_sum$att-data_sum$c*data_sum$att.se)
data_sum$ymax = ( data_sum$att + data_sum$c * data_sum$att.se )
est_df <- data.table( data_sum )

# Vectors of old and new column names
oldnames <- c( "year", "att", "ymax", "ymin" )
newnames <- c( "Time", "Estimate", 'upper_bound', 'lower_bound' )
setnames(est_df, oldnames, newnames)

# Getting the number of obseravtions
n_df <- df2_below[t2ev_ror_data_entry_block_any %in% est_df$Time][, .N, by = t2ev_ror_data_entry_block_any]
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
  geom_line(color = "blue") +
  geom_point(color = "red") +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), 
                width = 0.2, color = "red") +
  geom_rect(data = est_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, 
                               ymax = ymax), 
            fill = "grey70", colour = "grey69", alpha = 0.4, size = 0.2) + 
  coord_cartesian(ylim = c( min_y_lim, max_y_lim ) ) +
  scale_y_continuous(name = "Estimate", 
                     sec.axis = sec_axis(~(.+ (-1*min_y_lim)) * (1/scale_fac), 
                                         name = "Number of Observations")) + 
  labs(y = "ATT", x = "Time to First Year of Treatment", 
       title = "Ror Data Entry at Level") + 
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "gray") + 
  ggtitle( "Ror Data Entry at Block Level when Any \n Village is Treated (Below Median)" ) +
  theme_minimal()

filename <- paste0( "below_median_block_any.png")
path <- file.path("F_Figures/_vcf_10_callaway_santanna_analysis/group_states", 
                  filename)
ggsave(filename = path, plot = p, width = 6, height = 4, dpi = 300)

#-------------------------------------------------------------------------------








