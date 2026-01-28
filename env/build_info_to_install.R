#!/usr/bin/env RScript

# CREATION : 13 dec. 2021
# UPDATE : 13 janv. 2026
# GOAL : This script is used to create the file info_to_install.txt for Docker

rm(list = ls())

library(dplyr)

## What are the packages of interest ?
package_names = c("base", "rmdformats", "rmarkdown", "knitr", "plotly", "gridtext",
                  "ggplot2", "dplyr", "stringr", "RColorBrewer", "reshape2",
                  "SummarizedExperiment", "Matrix", "circlize", "kableExtra",
                  "ComplexHeatmap", "ggdendro", "DESeq2", "ggrepel", "FactoMineR") %>% unique()

## Installation order
info_to_install = aquarius::repro_installation_order(wanted_packages = package_names)

# Which packages have no url ?
info_to_install %>% dplyr::filter(url == "another_exception")

# Add exception
info_to_install["GenomeInfoDbData", "url"] = "https://bioconductor.org/packages/3.21/data/annotation/src/contrib/GenomeInfoDbData_1.2.14.tar.gz"

# Which packages have no url ?
info_to_install %>% dplyr::filter(url == "another_exception")

# Remove base packages
info_to_install = info_to_install %>% dplyr::filter(url != "this_is_base")

# Update order
info_to_install$order = c(1:nrow(info_to_install))

# Save directory
save_dir = "/import/oceans06/analyses/DESimple_A2025/4_env_Docker/"

# Save
today_date = format(Sys.time(), "%Y_%m_%d")
save_name = paste0(save_dir, "/info_to_install_", today_date, ".txt")

utils::write.table(x = info_to_install,
                   file = save_name,
                   quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Download all packages
out_dir = paste0(save_dir, "/", today_date)
dir.create(out_dir)

apply(info_to_install, 1, FUN = function(one_row) {
  package_name = one_row["package_name"]
  dl_url = one_row["url"]
  
  message(package_name)
  download.file(url = dl_url,
                destfile = paste0(out_dir, "/", package_name, ".tar.gz"),
                quiet = TRUE)
})


