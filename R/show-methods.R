setMethod(f="show",signature="Trajectories", function(object){
  cat("class: Trajectories")
  cat("\n\tSource: ");cat(object@path)
  cat("\n\tNames: ");cat(object@filenames)
  cat("\n\tTrajectories: ");cat(length(object@trajlist))
  cat("\n\tPoints per trajectory: ");for(i in 1:length(object@trajlist)) cat(nrow(object@trajlist[[i]]),"")
  cat("\n\tadf p-values: ");cat(object@avd)
  cat("\n\nAvailable slots: "); cat(slotNames(object))
  cat("\n")
})

setMethod(f="show",signature="TransTrajectories",function(object){
  cat("class: TransTrajectories")
  if(length(object@tmethod)!=0){
    cat("\n\tMethod: ");cat(object@tmethod)
    cat("\n\tNames: "); cat(object@tfilenames)
    cat("\n\tTrajectories: ");cat(length(object@tfilenames))
    cat("\n\tPoints per trajectory: ");for(i in 1:length(object@ttrajlist)) cat(nrow(object@ttrajlist[[i]]),"")
    cat("\n\tadf p-values: "); cat(object@tavd)
    if(object@tmethod=='differentiation'){
      cat("\n\tDifferentiated trajectories: "); for(i in length(object@difftraj)) cat(names(object@difftraj)[i])
    }
    else {
      cat("\n\tSegment splits per trajectory: ")
      if(length(object@breakpoints)==0) cat("No segments defined.") 
      else if(length(object@breakpoints)==1) cat(unlist(object@breakpoints)) 
      else {
        for(i in 1:length(object@breakpoints)) cat(length(object@breakpoints[[i]]),"")
      }
    }    
  }
  else{
    cat("\n\tNo transformation was done.")
  }
  cat("\n\n")
  callNextMethod()
})

setMethod(f="show",signature="SegTrajectories",function(object){
  cat("class: SegTrajectories")
  paramnames<-c("method","maxQ","filter.number","factor","threshold.level","minobs")
  params<-object@sparams
  cat("\n\tParameters: ");for(i in 1:length(params)) cat(paste(paramnames[i],"=",params[i],sep=""),"")
  cat("\n\t# Segments: ");cat(length(unique(object@smatrix[,5])))
  cat("\n\n")
  callNextMethod()
})

setMethod(f="show",signature="SegSeriesTrajectories",function(object){
  cat("class: SegSeriesTrajectories")
  paramnames<-c("intervention")
  params<-object@ssparams
  cat("\n\tParameters: ");for(i in 1:length(params)) cat(paste(paramnames[i],"=",params[i],sep=""),"")
  cat("\n\t# identified subpopulations: ");cat(length(unique(object@ssmatrix$subpopulation)))
  cat("\n\n")
  callNextMethod()  
})
