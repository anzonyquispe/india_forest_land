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
df1 <- fread("A_MicroData/data_sysdif.csv")

# Add a new column 'is_duplicate' which is TRUE for duplicated rows based on 'vill_id' and 'year'
df1[, is_duplicate := duplicated(.SD) | duplicated(.SD, fromLast = TRUE), .SDcols = c('vill_id', 'year')]
df1[ df1$is_duplicate == 0 ]


# Get year of treatment
# Get the minimum year where ps_ror_data_entry equals 1 for each vill_id
min_yr_tr <- df1[post_ror_data_entry == 1, .(min_year = min(year)), by = .(vill_id)]
# Merge with the original data.table by vill_id
df1 <- merge(df1, min_yr_tr, by = "vill_id", all.x = TRUE)
df1[is.na(min_year), min_year := 0]


# Calculate the median of the 'value' column
df_year_2010 <- df1[df1$year == 2010]
median_value <- median(df_year_2010$per_treecover)
# Generate a new variable 'above_median'; it will be TRUE if 'value' is greater than the median
df_year_2010[, above_median_all := (per_treecover > median_value)*1 ]
df_year_2010 <- df_year_2010[, .(shrid, above_median_all )]

df2 <- df1 %>% left_join( df_year_2010 )
df2_avobe <- df2[df2$above_median_all == 1]
df2_below <- df2[df2$above_median_all == 0]



#-------------------------------------------------------------------------------
# Avobe the median
#-------------------------------------------------------------------------------

# Start estimation
out <- att_gt(yname = "per_treecover",
              gname = "min_year",
              idname = "vill_id",
              tname = "year",
              xformla = ~ 1 ,
              data = df2_avobe,
              est_method = "reg", 
              control_group = c("nevertreated", "notyettreated")
)

filename <- paste0( "above_median.RDS")
path <- file.path("E_Estimates/_vcf_10_callaway_santanna_analysis/median_filter", 
                  filename)
saveRDS( out, path )



#-------------------------------------------------------------------------------
# Below the median
#-------------------------------------------------------------------------------


# Start estimation
out <- att_gt(yname = "per_treecover",
              gname = "min_year",
              idname = "vill_id",
              tname = "year",
              xformla = ~ 1 ,
              data = df2_below,
              est_method = "reg", 
              control_group = c("nevertreated", "notyettreated")
)

filename <- paste0( "bellow_median.RDS")
path <- file.path("E_Estimates/_vcf_10_callaway_santanna_analysis/median_filter", 
                  filename)
saveRDS( out, path )


