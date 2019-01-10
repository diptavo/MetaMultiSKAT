save.Meta.Files <- function(obj,pre,study,work.dir = NULL){
  if(is.null(work.dir)){
    work.dir <- getwd()
  }
  filename <- paste(work.dir,"/",pre,"_st.",study,".RDS",sep = "")
  if(class(obj) != "MultiSKAT.Test.Item")
    stop("Please provide an object of class MultiSKAT.Test.Item")
  saveRDS(object = obj,file = filename)
  print(paste("MultiSKAT.Test.Item saved in ",filename, " successfully",sep = ""))
}

read.Meta.Files <- function(file.vec){
  obj.list <- list()
  l1 <- length(file.vec)
  for(i in 1:l1){
    obj.list[[i]] <- readRDS(normalizePath(file.vec[i]))
  }
  return(obj.list)
}
