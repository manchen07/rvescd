path <- "chtc/main_job/results"
files <- list.files(path=path)
setwd(path)
all <- lapply(files, function(x) mget(load(x)))
results <- sapply(files, function(x) mget(load(x)[1]), simplify = TRUE)
res_dat <- do.call("rbind", results)
rownames(res_dat) <- c()

save(all, file = "../sim-res-all.Rdata")
save(res_dat, file = "../simulation-results.Rdata")
