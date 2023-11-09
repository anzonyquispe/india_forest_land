library(did)
library(data.table)
library(dplyr)
library(boot)
library(ssynthdid)
library( tidyr )

setwd("/scratch/gpfs/ar8787/groupdata2/india_forest_land")

# Input estimate
data <- fread("A_MicroData/data_sysdif.csv")

# Generation of year of treatment
data$tyear <- data$ev_ror_data_entry
data[is.na(data$tyear), 'tyear'] <-  0 

# Generation of deciles
# Convert your data frame to data.table (if it's not already)
data <- as.data.table(data)

# 1. Generate deciles for per_treecover for the year 2000
data[year == 2010, treedec := ntile(per_treecover, 10)]
# 2. For each village, calculate the maximum decile value
data[, treedec1 := max(treedec, na.rm = TRUE), by = vill_id]
# 3. Tabulate the max decile value and generate indicator variables
# Since data.table doesn't have a direct tabulate and generate functionality,
# we'll achieve this using a join
unique_treedec1 <- unique(data$treedec1)
for (val in unique_treedec1) {
  col_name <- paste0("itreedec", val)
  data[, (col_name) := 0L]
  data[treedec1 == val, (col_name) := 1L]
}
# 4. Drop the temporary decile variables
data[, c("treedec", "treedec1") := NULL]

# Add the per_treecover of 2010
data[year == 2010, treecover2010yr_aux := per_treecover]
data[, treecover2010yr := max( treecover2010yr_aux, na.rm = TRUE), 
     by = vill_id]
data[, c( "treecover2010yr_aux" ) := NULL]

# Order the dataset
setorder(data, vill_id, year)

# Start estimation
out <- att_gt(yname = "per_treecover",
              gname = "tyear",
              idname = "vill_id",
              tname = "year",
              xformla = ~ 1 + treecover2010yr,
              data = data,
              est_method = "reg", 
              control_group = c("nevertreated", "notyettreated")
)


saveRDS( out, 'E_Estimates/_vcf_10_callaway_santanna_analysis/3_espec_cs.RDS')
