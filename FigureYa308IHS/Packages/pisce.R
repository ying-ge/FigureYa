parameter <- commandArgs(T)

pisce.path <- parameter[1]
out.path <- parameter[2]
spec <- ifelse(is.na(parameter[3]), yes = "human", no = parameter[3])
key <- ifelse(is.na(parameter[4]), yes = "symbol", no = parameter[4])
P_threshold <- ifelse(is.na(parameter[5]), yes =  1e-8, no = as.numeric(parameter[5]))
nBootstrap <- ifelse(is.na(parameter[6]), yes =  200, no = as.numeric(parameter[6]))

message(quote(pisce.path), " = ", pisce.path)
message(quote(out.path), " = ", out.path)
message(quote(spec), " = ", spec)
message(quote(key), " = ", key)
message(quote(P_threshold), " = ", P_threshold)
message(quote(nBootstrap), " = ", nBootstrap)


pisce.script.path <- file.path(pisce.path, "script")
pisce.data.path <- file.path(pisce.path, "data")
collections = list.files(file.path(pisce.data.path, paste0(spec, "_", key)), 
                         pattern = "txt$|dat$", full.names = T)
message("collections:", paste(gsub("\\.txt|\\.dat", "", basename(collections)), collapse = ", "))
message("")

message(Sys.time())

system(paste("Rscript", file.path(pisce.script.path, "aracne_data-prep.r"),
             "--rds", file.path(out.path, "expr.rds"),
             "--out", file.path(out.path, "expr.tsv")))

for (collection in collections){
  name = gsub("(cotfs|sig|tfs|surface)(.+)", "\\1", basename(collection))
  message("Calculating ", name)
  if (!dir.exists(file.path(out.path, name))) dir.create(file.path(out.path, name))
  
  if (length(list.files(path = file.path(out.path, name), pattern = "^miThreshold"))==0){
    system(paste("java -Xmx5G -jar", file.path(pisce.script.path, "aracne-ap.jar"),
                 "-e", file.path(out.path, "expr.tsv"),
                 "-o", file.path(out.path, name),
                 "--tfs", collection,
                 "--pvalue", P_threshold,
                 "--seed", 666,
                 "--calculateThreshold"), ignore.stdout = T, ignore.stderr = T)
  }
  
  for (i in 1:nBootstrap){
    if (file.exists(file.path(out.path, name, paste0("bootstrapNetwork_seed", i, ".txt")))) next
    
    system(paste("java -Xmx5G -jar", file.path(pisce.script.path, "aracne-ap.jar"),
                 "-e", file.path(out.path, "expr.tsv"),
                 "-o", file.path(out.path, name),
                 "--tfs", collection,
                 "--pvalue", P_threshold,
                 "--seed", i), ignore.stdout = T, ignore.stderr = T)
    
    filename <- list.files(file.path(out.path, name), "^bootstrapNetwork")
    filename <- filename[!grepl("seed", filename)]
    
    if (length(filename) == 0) {
      warning("No bootstrapNetwork file found in ", file.path(out.path, name))
      next
    }
    
    if (length(filename) > 1) {
      warning("Multiple bootstrapNetwork files found, using the first one.")
      filename <- filename[1]
    }
    
    file.rename(from = file.path(out.path, name, filename),
                to = file.path(out.path, name, paste0("bootstrapNetwork_seed", i, ".txt")))
  }
  
  if(!file.exists(file.path(out.path, paste0(name, "_finalNetwork_4col.tsv")))){
    system(paste("Rscript", file.path(pisce.script.path, "aracne_consolidate.r"),
                 file.path(out.path, name),
                 file.path(out.path, "expr.tsv"),
                 collection,
                 "bonferroni 0.01"), ignore.stdout = T, ignore.stderr = T)
    invisible(file.copy(from = file.path(out.path, name, "finalNetwork_4col.tsv"),
                        to = file.path(out.path, paste0(name, "_finalNetwork_4col.tsv"))))
  }
}

message("Merging collections")
system(paste("bash", file.path(pisce.script.path, "aracne_merge.sh"),
             paste(list.files(file.path(out.path), pattern = "finalNetwork_4col.tsv$", full.names = T), collapse = " "),
             file.path(out.path, "finalNet_merged.tsv")), ignore.stdout = T, ignore.stderr = T)

system(paste("Rscript", file.path(pisce.script.path, "aracne_regProc.r"),
             "--a_file", file.path(out.path, "finalNet_merged.tsv"),
             "--exp_file", file.path(out.path, "expr.tsv"),
             "--out_dir", paste0(file.path(out.path), "/"),
             "--out_name", "out"), ignore.stdout = T, ignore.stderr = T)

