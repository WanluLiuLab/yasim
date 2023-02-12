#!/usr/bin/env Rscript
ip <- installed.packages()
write.csv(
    ip[, c("Package", "Version")],
    "installed.csv",
    row.names=FALSE
)
