#!/usr/bin/env Rscript 

argv <- commandArgs(TRUE)
from_file = file("stdin")
to_file = stdout()
correction_method = argv[1]

dataT <- read.table(from_file)
values <- data.frame(dataT)
newValues <- p.adjust(values$V1, method = correction_method)
write.table(newValues, to_file , FALSE, FALSE, "\t", "\n", "NA", ".", FALSE, FALSE)

