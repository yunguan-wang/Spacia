# One R process, one model load, multiple independent receiver-feature jobs.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("Usage: Rscript spacia_batch.R manifest.json")
jobs <- rjson::fromJSON(file=args[1])
if (!length(jobs)) quit(status=0)
source(file.path(jobs[[1]][[1]], "spacia_runtime.R"))
state <- new.env(parent=emptyenv())
for (job in jobs) spacia_run_job(unlist(job, use.names=FALSE), state)
