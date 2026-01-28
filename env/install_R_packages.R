#!/usr/bin/env Rscript

file_path = commandArgs(trailingOnly = TRUE)
info_to_install = utils::read.table(file_path, header = TRUE,
                                    colClasses = c("integer", "character", "character", "character", "character"))

for (i in c(1:nrow(info_to_install))) {
  package = info_to_install[i, "package_name"]

  install.packages(paste0("/packages_dir/", package, ".tar.gz"),
                   repos = NULL, type = "source",
                   clean = TRUE,
                   dependencies = c("Depends", "Imports", "LinkingTo"),
                   extra = "--no-check-certificate")
  
  if (!(package %in% rownames(utils::installed.packages()))) {
    stop(paste0("PACKAGE ", package, " NOT INSTALLED. Add the missing OS-level dependency in the definition file and build again."))
  }
}
