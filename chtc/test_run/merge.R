path <- "chtc/test_run/results"
files <- list.files(path=path)
setwd(path)
results <- sapply(files, function(x) mget(load(x)[1]), simplify = TRUE)
all <- do.call("rbind", results)
