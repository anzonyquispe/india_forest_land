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

# Generation of t2event
# Create a data.table with minimum year by block_id
min_years_data <- df1[
  post_ror_data_entry == 1,
  .(min_year = min(year)), # Replace 'year_column' with your actual year column name
  by = .(vill_id)
]

# Merge the min_years_data back to the original dataset
df1 <- merge(df1, min_years_data, by = "vill_id", all.x = TRUE)



# Use tstrsplit to split the 'shrid' column and assign the output to new columns:
df1[, c("pc11_s_id", "pc11_d_id", "pc11_sd_id", "pc11_tv_id") 
    := tstrsplit(shrid, "-", fixed=TRUE)]
cols_to_convert <- c("pc11_s_id", "pc11_d_id", "pc11_sd_id", "pc11_tv_id")
df1[, (cols_to_convert) := lapply(.SD, as.numeric), .SDcols = cols_to_convert]

# Import village
vill_sch <- fread("A_MicroData/vill_sch.csv")
cols_to_convert <- c("pc11_s_id", "pc11_d_id", "pc11_sd_id", "pc11_tv_id")
vill_sch[, (cols_to_convert) := lapply(.SD, as.numeric), .SDcols = cols_to_convert]


df2 <- df1 %>% inner_join( vill_sch, 
                           by = c("pc11_s_id", "pc11_d_id", "pc11_sd_id", "pc11_tv_id"))

df2_sch <- df2[df2$sch == 1]
df2_nonsch <- df2[df2$sch == 0]




#-------------------------------------------------------------------------------
# sch the median
#-------------------------------------------------------------------------------



# Start estimation
out <- att_gt(yname = "per_treecover",
              gname = "min_year",
              idname = "vill_id",
              tname = "year",
              xformla = ~ 1 ,
              data = df2_sch,
              est_method = "reg", 
              control_group = c(  "notyettreated")
)

filename <- paste0( "sch_median_vill.RDS")
path <- file.path("E_Estimates/_vcf_10_callaway_santanna_analysis/sch", 
                  filename)
saveRDS( out, path )


# aggregate the group-time average treatment effects
dynamic_sum <- aggte( out, type = "dynamic")
filename <- paste0( "sch_median_vill_sum.RDS")
path <- file.path("E_Estimates/_vcf_10_callaway_santanna_analysis/sch", 
                  filename)
saveRDS( out, path )


#-------------------------------------------------------------------------------
# nonsch the median
#-------------------------------------------------------------------------------
# Start estimation
out <- att_gt(yname = "per_treecover",
              gname = "min_year",
              idname = "vill_id",
              tname = "year",
              xformla = ~ 1 ,
              data = df2_nonsch,
              est_method = "reg", 
              control_group = c(  "notyettreated")
)

filename <- paste0( "nonsch_median_vill.RDS")
path <- file.path("E_Estimates/_vcf_10_callaway_santanna_analysis/sch", 
                  filename)
saveRDS( out, path )


# aggregate the group-time average treatment effects
dynamic_sum <- aggte( out, type = "dynamic")
filename <- paste0( "nonsch_median_vill_sum.RDS")
path <- file.path("E_Estimates/_vcf_10_callaway_santanna_analysis/sch", 
                  filename)
saveRDS( out, path )

