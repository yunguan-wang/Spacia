# Shared runtime for standalone jobs and bounded batches of jobs in one R process.
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(rjson))

spacia_cache_identity <- function(spacia_path) {
  existing_hashes <- function(paths) {
    paths <- unique(paths[nzchar(paths) & file.exists(paths)])
    tools::md5sum(paths)
  }
  identity <- list(
    source=existing_hashes(file.path(spacia_path, "Fun_MICProB_C2Cinter.cpp")),
    R=R.version, platform=Sys.info()[c("sysname", "release", "machine")],
    packages=vapply(c("Rcpp", "RcppArmadillo"), function(p) as.character(packageVersion(p)), ""),
    build_environment=Sys.getenv(c("R_MAKEVARS_USER", "R_MAKEVARS_SITE", "CC", "CXX",
                                  "CXXFLAGS", "PKG_CPPFLAGS", "PKG_CXXFLAGS", "PKG_LIBS")),
    makefiles=existing_hashes(c(file.path(R.home("etc"), "Makeconf"),
                               file.path(R.home("etc"), "Makevars.site"),
                               path.expand("~/.R/Makevars"),
                               path.expand(paste0("~/.R/Makevars-", R.version$platform)),
                               Sys.getenv("R_MAKEVARS_USER"), Sys.getenv("R_MAKEVARS_SITE"))))
  tmp <- tempfile()
  on.exit(unlink(tmp))
  saveRDS(identity, tmp, version=2)
  list(key=unname(tools::md5sum(tmp)), identity=identity)
}

spacia_load_model <- function(spacia_path) {
  spacia_path <- normalizePath(spacia_path, mustWork=TRUE)
  cache_root <- Sys.getenv("SPACIA_CACHE_DIR", unset=file.path(dirname(spacia_path), ".spacia_cache"))
  dir.create(cache_root, recursive=TRUE, showWarnings=FALSE)
  cache_root <- normalizePath(cache_root, mustWork=TRUE)
  identity <- spacia_cache_identity(spacia_path)
  cache_dir <- file.path(cache_root, identity$key)
  lock_dir <- paste0(cache_dir, ".lock")
  started <- Sys.time()
  repeat {
    if (dir.create(lock_dir, showWarnings=FALSE)) break
    if (as.numeric(difftime(Sys.time(), started, units="secs")) > 180)
      stop("Timed out waiting for compiler-cache lock: ", lock_dir,
           ". If the owning process has exited, remove that stale lock directory.")
    Sys.sleep(0.1)
  }
  on.exit(unlink(lock_dir, recursive=TRUE), add=TRUE)
  writeLines(c(paste("pid", Sys.getpid()), paste("host", Sys.info()[["nodename"]]),
               paste("started", Sys.time())), file.path(lock_dir, "owner.txt"))
  dir.create(cache_dir, recursive=TRUE, showWarnings=FALSE)
  model <- new.env(parent=globalenv())
  Rcpp::sourceCpp(file.path(spacia_path, "Fun_MICProB_C2Cinter.cpp"),
                  env=model, cacheDir=cache_dir, rebuild=FALSE, showOutput=TRUE)
  saveRDS(identity$identity, file.path(cache_dir, "build_identity.rds"))
  source(file.path(spacia_path, "MICProB_MIL_C2Cinter.R"), local=model)
  source(file.path(spacia_path, "MIL_wrapper.R"), local=model)
  model$spacia_path <- spacia_path
  cat(sprintf("SPACIA_MODEL_READY pid=%d cache=%s elapsed_seconds=%.3f\n",
              Sys.getpid(), cache_dir, as.numeric(difftime(Sys.time(), started, units="secs"))))
  model
}

spacia_run_job <- function(args, state=new.env(parent=emptyenv())) {
  if (!length(args) %in% c(12L,13L,14L)) stop("Expected 12, 13 or 14 Spacia job arguments")
  spacia_path <- normalizePath(args[1], mustWork=TRUE)
  output_path <- args[10]
  job_id <- args[5]
  plot_mcmc <- as.logical(args[11])
  if (is.na(plot_mcmc)) stop("Invalid plot_mcmc value")
  dir.create(output_path, recursive=TRUE, showWarnings=FALSE)
  depth <- sink.number()
  sink(file.path(output_path, paste0(job_id, "_log.txt")))
  on.exit(while (sink.number() > depth) sink(), add=TRUE)
  cat(sprintf("SPACIA_WORKER_PID %d\n", Sys.getpid()))
  if (is.null(state$model) || !identical(state$model$spacia_path, spacia_path))
    state$model <- spacia_load_model(spacia_path)
  model <- state$model
  if (plot_mcmc && !exists("BetaB2MCMCPlots", envir=model, inherits=FALSE))
    source(file.path(spacia_path, "BetaB2MCMCPlots.R"), local=model)

  # Cache only the most recent common inputs, bounding memory per worker.
  # Content hashes prevent stale reuse if a caller replaces a file in place.
  distance_scale <- if (length(args)==14L) as.numeric(args[14]) else NULL
  if (!is.null(distance_scale) && (!is.finite(distance_scale) || distance_scale <= 0))
    stop("Explicit distance normalizer must be finite and positive")
  input_key <- list(files=tools::md5sum(normalizePath(args[2:3], mustWork=TRUE)),
                    distance_scale=distance_scale)
  if (is.null(state$input_key) || !identical(state$input_key, input_key)) {
    tmp <- rjson::fromJSON(file=args[2])
    state$sender <- lapply(tmp, function(x) do.call(rbind, as.list(x)))
    distance <- rjson::fromJSON(file=args[3])
    max_dist <- if (is.null(distance_scale)) max(sapply(distance, function(x) x[which.max(abs(x))])) else distance_scale
    state$distance <- lapply(distance, function(x) x/max_dist)
    state$input_key <- input_key
    cat("SPACIA_INPUTS_LOADED\n")
  } else cat("SPACIA_INPUTS_REUSED\n")
  receiver <- read.csv(args[4], header=FALSE, row.names=NULL, stringsAsFactors=FALSE)$V1 == 1
  ntotal <- as.integer(args[6]); nwarm <- as.integer(args[7])
  nthin <- as.integer(args[8]); nchain <- as.integer(args[9])
  prior <- if (length(args)>=13L) as.numeric(args[13]) else 1
  thetas <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  # Match independent Rscript jobs, regardless of which worker handles a job.
  set.seed(0)
  t0 <- Sys.time()
  res <- model$MIL_C2Cinter(receiver, state$distance, state$sender,
                           ntotal,nwarm,nthin,nchain,thetas,prior)
  print(Sys.time()-t0)
  print(gc())
  for (name in names(res)) {
    path <- file.path(output_path,paste0(job_id,"_",name,".txt"))
    if (name=="FDRs") {
      fdr <- res$FDRs; fdr[is.na(fdr)] <- 1
      write.table(fdr,path,sep="\t")
    } else write.table(res[name],path,sep="\t")
  }
  if (plot_mcmc) {
    beta <- as.matrix(res$beta); b <- as.matrix(res$b)
    colnames(beta) <- paste0("beta.",seq_len(ncol(beta)))
    colnames(b) <- c("b.1","b.2")
    model$BetaB2MCMCPlots(beta,b,nwarm,ntotal,nthin,nchain,job_id,output_path,args[12])
  }
  invisible(NULL)
}
