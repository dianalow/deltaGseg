#########################################################
# Non S4                                                #
#########################################################

parseTraj<-function(path=getwd(),files=NULL,fromfile=TRUE){
  # gets trajectories from path folder, if files are not specified, assume all files in dir will be used.  
  # computes adf.test
  # define default break as NULL, so that if you plot without splitting, breaks would have been defined
  # returns: 1) path (for tracking), 2) nTraj, 3) trajlist, 4) adf.test 5) status of breakpoints (for plotting)
  options(warn=-1)
  if(is.null(files)) files<-list.files(path)
  ntraj<-length(files)
  trajlist<-list()
  avd<-rep(0,ntraj)
  if(Sys.info()["sysname"]=="Windows") fullpath=paste(path,"\\",sep="") else fullpath=paste(path,"/",sep="")
  for(i in 1:ntraj){
    if(fromfile) dat<-as.matrix(read.table(paste(fullpath,files[i],sep=""),header=F))
    else dat<-files[[i]]
    if(ncol(dat)!=2 | !is.numeric(dat)) stop(paste("Please check if",files[i],"is in correct format."))
    avd[i]<-adf.test(dat[,2])$p.value
    trajlist[[i]]<-dat
  }
  if(fromfile) {
    names(trajlist)<-files
    return(new(Class="Trajectories",path=fullpath,filenames=files,trajlist=trajlist,avd=avd))
  }
  else{
    return(new(Class="Trajectories",path=fullpath,filenames=as.character(seq(1:ntraj)),trajlist=trajlist,avd=avd))
  }
}

######################
# internal functions #
######################

segden1<-function(data,cpoint_low,cpoint_high,fn=1,factor,thresh_level,segment_number,series_number){
  xnew<-c(data,data[(length(data)-1):1])[1:2^ceiling(log(length(data),2))]
  wdy1<-wd(xnew, filter.number=fn, family="DaubExPhase",type="wavelet",bc="symmetric", verbose=F, min.scale=0,precond=T)
  th<-threshold(wdy1,return.threshold=TRUE,by.level=thresh_level)
  filty1<-threshold(wdy1, value = factor*th, policy="manual")
  filtered <-wr(filty1)
  b<-filtered[1:length(data)]
  c<-data-b
  dlow<-rep(cpoint_low,length(c))
  dhigh<-rep(cpoint_high,length(c))
  e<-rep(segment_number,length(c))
  f<-t(matrix(rep(quantile(b,seq(0,1,0.05)),length(c)),nrow=21))
  g<-rep(series_number,length(c))
  all<-matrix(cbind(data,b,c,dlow,dhigh,e,f,g),ncol=28)
  return(all)
}

join.segmented.series<-function(data){
  data1<-data[[1]]
  if(length(data)>1){
    for(i in 2:length(data)){
      dd<-data[[i]]
      dd[,4]<-dd[,4]+data1[nrow(data1),5]
      dd[,5]<-dd[,5]+data1[nrow(data1),5]
      data1<-matrix(rbind(data1,dd),ncol=ncol(data1))
    }
  }
  return(data1)
}

mydist<-function(data){
  metrics<-c(data)
  dist1<-as.matrix(abs(dist(metrics)),ncol=length(metrics))
  start_seg<-seq(1,nrow(dist1),nrow(data))
  end_seg<-seq(nrow(data),nrow(dist1),nrow(data))
  dist2<-matrix(0,ncol(data),ncol(data))
  
  for(i in 1:(length(start_seg)-1)){
    for(j in (i+1):length(end_seg)){
      dist2[i,j]<-mean(dist1[start_seg[i]:end_seg[i],start_seg[j]:end_seg[j]])
      dist2[j,i]<-dist2[i,j]
    }
  }
  distances<-as.dist(dist2)
  return(distances)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
