library(did)
library(data.table)
library(dplyr)
library(boot)
library(ssynthdid)
library( tidyr )
library(Matrix)


setwd("/scratch/gpfs/ar8787/groupdata2/india_forest_land")


# # Load estimation
espec1 <- readRDS( 'E_Estimates/_vcf_10_callaway_santanna_analysis/output_1.RDS')
espec1_plot <- aggte(espec1, type = "dynamic")
saveRDS(espec1_plot, 'E_Estimates/_vcf_10_callaway_santanna_analysis/output_1.RDS' )

# Load estimation
espec2 <- readRDS( 'E_Estimates/_vcf_10_callaway_santanna_analysis/2_espec_cs.RDS')
espec2_plot <- aggte(espec2, type = "dynamic")
saveRDS(espec2_plot, 'E_Estimates/_vcf_10_callaway_santanna_analysis/2_espec_cs_ag.RDS' )

# Load estimation
espec3 <- readRDS( 'E_Estimates/_vcf_10_callaway_santanna_analysis/3_espec_cs.RDS')
espec3_plot <- aggte(espec3, type = "dynamic")
saveRDS(espec3_plot, 'E_Estimates/_vcf_10_callaway_santanna_analysis/3_espec_cs_ag.RDS' )

# # Load estimation
espec4 <- readRDS( 'E_Estimates/_vcf_10_callaway_santanna_analysis/4_espec_cs.RDS')
espec4_plot <- aggte(espec4, type = "dynamic")
saveRDS(espec4_plot, 'E_Estimates/_vcf_10_callaway_santanna_analysis/4_espec_cs_ag.RDS' )

# Load estimation
espec5 <- readRDS( 'E_Estimates/_vcf_10_callaway_santanna_analysis/5_espec_cs.RDS')
espec5_plot <- aggte(espec5, type = "dynamic")
saveRDS(espec5_plot, 'E_Estimates/_vcf_10_callaway_santanna_analysis/5_espec_cs_ag.RDS' )


# # Load estimation
espec6 <- readRDS( 'E_Estimates/_vcf_10_callaway_santanna_analysis/6_espec_cs.RDS')
espec6_plot <- aggte(espec6, type = "dynamic")
saveRDS(espec6_plot, 'E_Estimates/_vcf_10_callaway_santanna_analysis/6_espec_cs_ag.RDS' )




