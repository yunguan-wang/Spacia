# Backward-compatible standalone job entry point, now with a persistent build cache.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1L) stop("Missing Spacia code directory")
source(file.path(args[1], "spacia_runtime.R"))
spacia_run_job(args)
