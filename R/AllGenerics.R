#############
# Accessors #
#############

setGeneric("getAVD",function(object){standardGeneric ("getAVD")})

setMethod("getAVD","Trajectories",function(object){
  return(object@avd)
})

setMethod("getAVD","TransTrajectories",function(object){
  return(object@tavd)
})

setGeneric("getBreaks",function(object){standardGeneric ("getBreaks")})

setMethod("getBreaks","TransTrajectories",function(object){
  return(object@breakpoints)
})

setGeneric("getSegments",function(object){standardGeneric ("getSegments")})

setMethod("getSegments","SegTrajectories",function(object){
  return(object@smatrix)
})

setMethod("getSegments","SegSeriesTrajectories",function(object){
  return(object@ssmatrix)
})

setGeneric("getTNames",function(object){standardGeneric ("getTNames")})

setMethod("getTNames","Trajectories",function(object){
  return(object@filenames)
})

setMethod("getTNames","TransTrajectories",function(object){
  if(length(object@tfilenames)==0) return(object@filenames) else return(object@tfilenames)
})

setGeneric("getTraj",function(object){standardGeneric ("getTraj")})

setMethod("getTraj","Trajectories",function(object){
  return(object@trajlist)
})

setMethod("getTraj","TransTrajectories",function(object){
  return(object@ttrajlist)
})

#############
# deltaGseg #
#############

setGeneric("chooseBreaks",
           function(breakpoints,numbreaks){
             standardGeneric("chooseBreaks")
           })

setMethod("chooseBreaks",
          signature=c(breakpoints="list",numbreaks="numeric"),
          definition=function(breakpoints,numbreaks){
            # returns evenly spaced breakpoints given desired number of breaks
            # not recommended for accuracy, use only if you do not know how to choose breakpoints
            bl<-list()
            for(n in 1:length(breakpoints)){
              width<-length(breakpoints[[n]])/numbreaks
              x<-rep(0,0,numbreaks)
              for(i in 1:numbreaks) {
                x[i]<-breakpoints[[n]][(i-1)*width+1]
              }
              bl[[n]]<-x
            }
            names(bl)<-names(breakpoints)
            return(bl)  
          })

setGeneric("clusterPV",
           function(object,bootstrap=500){
             standardGeneric("clusterPV")
           })

setMethod("clusterPV",
          signature=c(object="SegTrajectories"),
          definition=function(object,bootstrap=500){
            data<-object@smatrix
            pv<-pvclust_new(data=t(unique(data[,7:27])), nboot=bootstrap)
            return(pv)
          })
setGeneric("clusterSegments",
           function(object,intervention="groups",pv=NULL,graphics=NULL){
             standardGeneric("clusterSegments")
           })

setMethod("clusterSegments",
          signature=c(object="SegTrajectories"),
          definition=function(object,intervention="groups",pv=NULL,graphics=NULL){
            #check arguments
            if(!is.element(intervention,c('groups','pvclust'))) stop(paste(intervention,'is not a valid input for argument "intervention".'))
            
            data<-object@smatrix
            changepoints<-c(c(1,which(diff(data[,28])!=0)+1))
            ser<-unique(data[,28])
            
            if(is.null(graphics)) graphics<-hue_pal()(length(unique(data[,6]))) else graph<-graphics
            graph<-c()
            for(i in 1:length(ser)){
              ser_data<-unique(data[which(data[,28]==ser[i]),5])
              graph<-c(graph,graphics[1:length(ser_data)])
            }
            ll<-paste("s",1:length(unique(data[,5])),sep="")
            cp<-c(0,unique(data[,5]))
            
            if(intervention=="pvclust"){
              hc<-hclust(mydist(t(unique(data[,7:27]))),"average")
              cc<-c()
              message("Segment grouping. Click on the root of the groups you want clustered.")
              message("Please ensure that ALL segments are grouped (boxed).\nOtherwise, function will not exit.")
              message("To exit, click Esc (Windows/Linux) or Ctrl-click (Mac)")
              layout(matrix(c(1,2), 2, 1, byrow = TRUE))
              
              while(length(as.numeric(unlist(cc)))!=(length(cp)-1)){
                plot(1:cp[2],data[1:cp[2],1],xlab="Time",ylab="FreeEnergy",xlim=c(0,nrow(data)),ylim=c(min(data[,1]),(max(data[,1])+5)),col=graph[1],pch=20,cex=0.5,main="")
                lines(1:cp[2],data[1:cp[2],2],col=graphics[1],lwd=2)
                for(i in 2:(length(cp)-1)){
                  points((cp[i]+1):cp[(i+1)],data[(cp[i]+1):cp[(i+1)],1],col=graph[i],pch=20,cex=0.5)
                  lines((cp[i]+1):cp[(i+1)],data[(cp[i]+1):cp[(i+1)],2],col=graph[i],lwd=2)
                }
                abline(v=changepoints[2:length(changepoints)],col="gray",lty=2)
                posy<-max(data[,1])+3
                for(i in 1:length(ll)){
                  text(mean(c(cp[i],cp[(i+1)])),posy,ll[i],cex=0.7)
                }
                plot(hc,labels=ll,xlab="",sub="",main="Hierarchical clustering with P-values")
                if(!is.null(pv)) text.pvclust(pv,cex=0.7)
                cc<-identify(hc,N=length(hc$order),MAXCLUSTER=length(hc$order))
                
                #case where nothing is selected, return 1 subpopulation
                if(length(unlist(cc))==0) cc<-list(hc$order)
              }
              
              rescc<-matrix(0,2,1)
              for(i in 1:length(cc)){
                rescc<-matrix(cbind(rescc,matrix(rbind(cc[[i]],(rep(i,length(cc[[i]])))),nrow=2)),nrow=2)
              }
              rescc<-rescc[,-1]
              rescc<-rescc[,sort.list(rescc[1,])]
              ct<-rescc[2,]
            }
            
            else if(intervention=='groups'){
              par(mfrow=c(2,1))
              plot(1:cp[2],data[1:cp[2],1],xlab="Time",ylab="FreeEnergy",xlim=c(0,nrow(data)),ylim=c(min(data[,1]),(max(data[,1])+5)),col=graph[1],pch=20,cex=0.5,main="")
              lines(1:cp[2],data[1:cp[2],2],col=graphics[1],lwd=2)
              for(i in 2:(length(cp)-1)){
                points((cp[i]+1):cp[(i+1)],data[(cp[i]+1):cp[(i+1)],1],col=graph[i],pch=20,cex=0.5)
                lines((cp[i]+1):cp[(i+1)],data[(cp[i]+1):cp[(i+1)],2],col=graph[i],lwd=2)
              }
              abline(v=changepoints[2:length(changepoints)],col="gray",lty=2)
              posy<-max(data[,1])+3
              for(i in 1:length(ll)){
                text(mean(c(cp[i],cp[(i+1)])),posy,ll[i],cex=0.7)
              }
              
              dd<-unique(data[,7:27])
              distances<-mydist(t(dd))
              hc<-hclust(distances,"average")
              plot(hc,labels=ll,main="Hierarchical Clustering (Euclidean)",xlab="",sub="")
              if(intervention=="heights"){
                options(warn=-1)
                switch1<-0
                while(switch1==0){
                  threshold_loc<-locator()
                  if(max(threshold_loc$y)<max(distances) | length(threshold_loc$y)==0){
                    switch1<-1
                  }
                }
                if(length(threshold_loc$y)>0){
                  threshold<-threshold_loc$y[length(threshold_loc$y)]
                  ct<-cutree(hc,h=threshold)
                  abline(h=threshold,lt=2,col=2)
                  mtext(paste("Manually identified groups = ",length(unique(ct))," at threshold ",threshold,sep=""), side=1, line=3)
                }
              }
              else if(intervention=="groups"){
                gg<-0
                while(gg<1){
                  gg<- readline("enter the number of groups (a positive integer): ")
                }
                ct<-cutree(hc,k=gg)                  
                mtext(paste("User defined identified groups = ",gg,sep=""), side=1, line=3)
              }
              
            }
            
            uct<-unique(ct)
            for(i in 1:length(uct)){
              wh<-which(ct==uct[i])
              mm<-c()
              for(j in 1:length(wh)){
                mm<-c(mm,which(data[,4]==cp[wh[j]] & data[,5]==cp[(wh[j]+1)]))
              }
              data[mm,6]<-i
            }
            
            gg<-unique(data[,6])
            for(i in 1:length(gg)){
              wh<-which(data[,6]==gg[i])
              mm<-median(data[wh,7])
              data[wh,7]<-mm
            }
            
            data<-data.frame(cbind(data[,1:3],data[,6],paste(getTNames(object)[data[,ncol(data)]]))) 
            colnames(data)<-c("observed","estimated","residuals","subpopulation","seriesID")
            tempTraj<-new("SegSeriesTrajectories",ct=ct,ssmatrix=data,ssparams=c(intervention))
            as(tempTraj,"Trajectories") <- as(object,"Trajectories")
            as(tempTraj,"TransTrajectories") <- as(object,"TransTrajectories")
            as(tempTraj,"SegTrajectories") <- as(object,"SegTrajectories")
            plot(tempTraj)
            return(tempTraj)  
          })

setGeneric("denoiseSegments",
           function(object,seg_method="BinSeg",maxQ=15,fn=1,factor=0.8,thresh_level=TRUE,minobs=200){
             standardGeneric("denoiseSegments")
           })

setMethod("denoiseSegments",
          signature=c(object="Trajectories"),
          definition=function(object,seg_method="BinSeg",maxQ=15,fn=1,factor=0.8,thresh_level=TRUE,minobs=200){
            #check arguments
            if(!is.element(seg_method,c('BinSeg','SegNeigh'))) stop(paste(seg_method,'is not a valid input for argument "seg_method".'))
            
            #check series length
            files<-getTraj(object)
            message("Checking series length against minobs.....")
            for(i in 1:length(files)) 
              if(nrow(files[[i]])<2*minobs){
                message(paste("Series length less than 2*minobs!\n",getTNames(object)[i],"=",nrow(files[[i]]) ))
                stop("Either define a longer series (minimum 2*minobs) or reduce the minobs parameter (length(series)/2)")
              }
            message("Start denoising.....")
            totalavd<-getAVD(object)
            filenames<-getTNames(object)
            if(length(totalavd[totalavd>0.05])>0) {
              message('Non-stationary segments detected!');cat(filenames[totalavd>0.05])
              stop("Splitting or differentiation is required.")
            }
                    
            results<-as.list(rep(0,length(files)))
            numfiles<-length(files)
            pb <- txtProgressBar(style=3)
            for(ff in 1:numfiles){
              setTxtProgressBar(pb, ff/numfiles)
              results[[ff]]<-matrix(0,1,28)
              dat<-files[[ff]][,2]
              options(warn=-1)
              adf<-totalavd[ff]
              
              if(class(object)=="TransTrajectories") {
                
                if(filenames[ff] %in% names(object@difftraj)) {
                  dat0<-dat
                  dat<-object@difftraj[[filenames[ff]]][,2]
                }
              }
              #changed due to changepoint function wrap!
              res<-cpt.mean(dat,method=seg_method,penalty="Asymptotic",pen.value=0.01,Q=maxQ,class=FALSE)
              
              if(seg_method=="BinSeg") {
                mm_est_binseg<-res$cps[1,1:res$op.cpts]
                cp<-sort(c(0,mm_est_binseg,length(dat)))
              }
              
              else if(seg_method=="SegNeigh") {
                mm_est_segneigh<-res$cps[(res$op.cpts+1),1:res$op.cpts]
                cp<-sort(c(0,mm_est_segneigh,length(dat)))                
              }

              dd<-diff(cp)
              for(i in 1:length(dd)) if(dd[i]<15) cp[(i+1)]<- -999
              cp<-cp[cp!=-999]

              if(length(cp)==1) cp<-c(0,length(dat))
              else cp[length(cp)]<-length(dat)

              dd<-diff(cp)
              if(class(object)=="TransTrajectories") {
                if(filenames[ff] %in% names(object@difftraj)) dat<-dat0
              }

              for(i in 1:(length(cp)-1)){
                a<-dat[(cp[i]+1):cp[(i+1)]]
                all<-segden1(data=a,cpoint_low=cp[i],cpoint_high=cp[(i+1)],fn,factor,thresh_level,segment_number=i,series_number=ff)
                results[[ff]]<-matrix(rbind(results[[ff]],all),ncol=28)
              }
              
              results[[ff]]<-results[[ff]][-1,]

              k<-1
              while(length(dd[dd<minobs])>0){
                if(k==1 & dd[k]<minobs){
                  cp[(k+1)]<-cp[(k+2)]
                  cp<-unique(cp)
                  a<-dat[(cp[k]+1):cp[(k+1)]]
                  all<-segden1(data=a,cpoint_low=cp[k],cpoint_high=cp[(k+1)],fn,factor,thresh_level,segment_number=k,series_number=ff)
                  results[[ff]][(cp[k]+1):cp[(k+1)],]<-all
                  dd<-diff(cp)
                }
                else if(k==length(dd) & dd[k]<minobs){
                  cp[k]<-cp[(k+1)]
                  cp<-unique(cp)
                  a<-dat[(cp[(k-1)]+1):cp[k]]
                  sn<-results[[ff]][which(results[[ff]][,4]==results[[ff]][(cp[(k-1)]+1),4] & results[[ff]][,5]==(results[[ff]][cp[k],5]-dd[k])),6][1]
                  all<-segden1(data=a,cpoint_low=cp[(k-1)],cpoint_high=cp[k],fn,factor,thresh_level,segment_number=sn,series_number=ff)
                  results[[ff]][(cp[(k-1)]+1):cp[k],]<-all
                  dd<-diff(cp)  
                }
                
                else if(k>1 & k<length(dd) & dd[k]<minobs){
                  pre_b<-results[[ff]][(cp[(k-1)]+1),17]
                  centerb<-results[[ff]][(cp[k]+1),17]
                  post_b<-results[[ff]][(cp[(k+1)]+1),17]
                  if(abs(centerb-pre_b)>abs(centerb-post_b)){
                    cp[(k+1)]<-cp[(k+2)]
                    cp<-unique(cp)
                    a<-dat[(cp[k]+1):cp[(k+1)]]
                    all<-segden1(data=a,cpoint_low=cp[k],cpoint_high=cp[(k+1)],fn,factor,thresh_level,segment_number=(results[[ff]][(cp[(k-1)]+1),6]+1),series_number=ff)
                    results[[ff]][(cp[k]+1):cp[(k+1)],]<-all
                    if(cp[(k+1)]<length(dat)){
                      results[[ff]][(cp[(k+1)]+1):nrow(results[[ff]]),6]<-results[[ff]][(cp[(k+1)]+1):nrow(results[[ff]]),6]-1
                    }
                    dd<-diff(cp)
                  }
                  
                  else if(abs(centerb-pre_b)<abs(centerb-post_b)){
                    cp[k]<-cp[(k+1)]
                    cp<-unique(cp)
                    a<-dat[(cp[(k-1)]+1):cp[k]]
                    sn<-results[[ff]][which(results[[ff]][,4]==results[[ff]][(cp[(k-1)]+1),4] & results[[ff]][,5]==(results[[ff]][cp[k],5]-dd[k])),6][1]
                    all<-segden1(data=a,cpoint_low=cp[(k-1)],cpoint_high=cp[k],fn,factor,thresh_level,segment_number=sn,series_number=ff)
                    results[[ff]][(cp[(k-1)]+1):cp[k],]<-all
                    results[[ff]][(cp[k]+1):nrow(results[[ff]]),6]<-results[[ff]][(cp[k]+1):nrow(results[[ff]]),6]-1
                    dd<-diff(cp)
                  }
                  
                }
                
                else if(k<length(dd) & dd[k]>minobs){
                  results[[ff]][,6]<-as.numeric(factor(results[[ff]][,6]))
                  k<-k+1
                }
              }
            }
            
            if(length(results)>1) results<-join.segmented.series(results)
            tempTraj<-new("SegTrajectories",sparams=c(seg_method,maxQ,fn,factor,thresh_level,minobs),smatrix=data.frame(results))
            as(tempTraj,"Trajectories") <- as(object,"Trajectories")
            if(class(object)=="TransTrajectories") as(tempTraj,"TransTrajectories") <- as(object,"TransTrajectories")
            message("Complete.")
            return(tempTraj)
          })

setGeneric("diagnosticPlots",
           function(object,norm.test="KS",single.series=FALSE){
             standardGeneric("diagnosticPlots")
           })

setMethod("diagnosticPlots",
          signature=c(object="SegSeriesTrajectories"),
          definition=function(object,norm.test="KS",single.series=FALSE){
            #check arguments
            if(!is.element(norm.test,c('KS','Shapiro','Agost'))) stop(paste(norm.test,'is not a valid input for argument "norm.test".'))
            
            files<-getTNames(object)
            data<-object@ssmatrix
            uu<-as.character(unique(data$seriesID))
            dd<-data.frame(data[,c(3,5)])
            colnames(dd)<-c("residuals","series")
            
            if(single.series==FALSE){
              for(i in 1:length(uu)){
                dd0<-dd[which(dd[,2]==uu[i]),]
                if(norm.test=="KS"){
                  pp<-lillieTest(as.numeric(as.vector(data$residuals))[which(data$seriesID==uu[i])])@test$p.value[[1]]
                  ttext<-paste("KS p-value = ",round(pp,3))
                  textplot<-ggplot(data.frame(ttext)) +
                    geom_text(data=data.frame(ttext),hjust=0,aes(x=0.4,y=0.8,label=ttext))
                }
                else if(norm.test=="Shapiro"){
                  pp<-shapiroTest(as.numeric(as.vector(data$residuals))[which(data$seriesID==uu[i])])@test$p.value[[1]]
                  mtext<-paste("Shapiro p-value = ",round(pp,3))
                }
                else if(norm.test=="Agost"){
                  pp<-as.numeric(dagoTest(as.numeric(as.vector(data$residuals))[which(data$seriesID==uu[i])])@test$p.value)
                  ttext<-c()
                  ttext.1<-c("D'Agostino p-value = ","Skewness p-value = ","Kurtosis p-value = ")
                  for(i in 1:3){
                    ttext<-rbind(ttext,paste(ttext.1[i],round(pp[[i]],3)))
                  }
                  textplot<-ggplot(data.frame(ttext)) +
                    geom_text(data=data.frame(ttext),hjust=0,aes(x=c(0.3,0.3,0.3),y=c(0.4,0.6,0.8),label=ttext))
                }
                textplot<-textplot+ylim(0,1)+xlim(0,1)+
                  theme(axis.title.x = element_blank(),axis.title.y = element_blank(), panel.border = element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank(),axis.text.y = element_blank())+
                  theme(panel.grid.minor=element_line(colour=NA), panel.grid.major=element_line(colour=NA),panel.background=element_rect(fill=NA, colour=NA))
                
                p1<-ggplot(dd0,aes(x=residuals)) + 
                  geom_histogram(binwidth=.5, colour="gray25", fill="cornflowerblue") +
                  ggtitle(paste("Residuals of",files[i])) 
                
                actvals<-acf(dd0[,1],ci=0,ci.type="ma",plot=FALSE)
                significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(dd0[,1])))
                #ggplot prep
                ff<-data.frame(actvals$acf)
                ff<-cbind(1:dim(ff)[1],ff)
                names(ff)[1]<-c("lag")
                ff<-melt(ff,id="lag")
                
                p2<-qplot(lag,value,data=ff,geom="bar",stat="identity",position="identity",xlab="Lag",ylab="ACF",main="Autocorrelations") +
                  geom_hline(yintercept=0,color="red",alpha=0.9,size=0.2) +
                  geom_hline(yintercept=-significance_level,color="blue",alpha=0.9,size=0.3,linetype="dashed") +
                  geom_hline(yintercept=significance_level,color="blue",alpha=0.9,size=0.3,linetype="dashed")    	  
                grid.newpage()
                pushViewport(viewport(layout=grid.layout(4,2)))
                print(p1, vp = viewport(layout.pos.row=1:3, layout.pos.col=1))
                print(p2, vp = viewport(layout.pos.row=1:3, layout.pos.col=2))
                print(textplot, vp = viewport(layout.pos.row=4, layout.pos.col=1))
              }
            }
            
            if(single.series==TRUE){
              if(norm.test=="KS"){
                pp<-lillieTest(as.numeric(as.vector(data$residuals)))@test$p.value[[1]]
                ttext<-paste("KS p-value = ",round(pp,3))
                textplot<-ggplot(data.frame(ttext)) +
                  geom_text(data=data.frame(ttext),hjust=0,aes(x=0.4,y=0.8,label=ttext))
              }
              else if(norm.test=="Shapiro"){
                pp<-shapiroTest(as.numeric(as.vector(data$residuals)))@test$p.value[[1]]
                ttext<-paste("Shapiro p-value = ",round(pp,3))
                textplot<-ggplot(data.frame(ttext)) +
                  geom_text(data=data.frame(ttext),hjust=0,aes(x=0.4,y=0.8,label=ttext))
              }
              else if(norm.test=="Agost"){
                pp<-as.numeric(dagoTest(as.numeric(as.vector(data$residuals)))@test$p.value)
                ttext<-c()
                ttext.1<-c("D'Agostino p-value = ","Skewness p-value = ","Kurtosis p-value = ")
                for(i in 1:3){
                  ttext<-rbind(ttext,paste(ttext.1[i],round(pp[[i]],3)))
                }
                textplot<-ggplot(data.frame(ttext)) +
                  geom_text(data=data.frame(ttext),hjust=0,aes(x=c(0.3,0.3,0.3),y=c(0.4,0.6,0.8),label=ttext))
              }
              textplot<-textplot+ylim(0,1)+xlim(0,1)+
                theme(axis.title.x = element_blank(),axis.title.y = element_blank(), panel.border = element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank(),axis.text.y = element_blank())+
                theme(panel.grid.minor=element_line(colour=NA), panel.grid.major=element_line(colour=NA),panel.background=element_rect(fill=NA, colour=NA))
              p1<-ggplot(dd,aes(x=residuals)) + 
                geom_histogram(binwidth=.5, colour="gray25", fill="cornflowerblue") +
                ggtitle("Distribution") 
              actvals<-acf(data[,3],ci=0,ci.type="ma",plot=FALSE)
              significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(data[,3])))
              ##ggplot prep
              ff<-data.frame(actvals$acf)
              ff<-cbind(1:dim(ff)[1],ff)
              names(ff)[1]<-c("lag")
              ff<-melt(ff,id="lag")
              p2<-qplot(lag,value,width=0.2,data=ff,geom="bar",stat="identity",position="identity",xlab="Lag",ylab="ACF",main="Autocorrelations")+
                geom_hline(yintercept=0,color="red",alpha=0.9,size=0.2) +
                geom_hline(yintercept=-significance_level,color="blue",alpha=0.9,size=0.3,linetype="dashed") +
                geom_hline(yintercept=significance_level,color="blue",alpha=0.9,size=0.3,linetype="dashed")
              grid.newpage()
              pushViewport(viewport(layout=grid.layout(4,2)))
              print(p1, vp = viewport(layout.pos.row=1:3, layout.pos.col=1))
              print(p2, vp = viewport(layout.pos.row=1:3, layout.pos.col=2))
              print(textplot, vp = viewport(layout.pos.row=4, layout.pos.col=1))
            }
            
          })

setGeneric("getIntervals",
           function(object){
             standardGeneric("getIntervals")
           })

setMethod("getIntervals",
           signature=c(object="SegSeriesTrajectories"),
           definition=function(object){
             dd<-getSegments(object)
             changepoints<-c(c(1,which(diff(dd$subpopulation)!=0)+1))
             subpop<-dd[changepoints,4]
             sp<-list(subpop,changepoints)
             names(sp)<-c('subpopulation','interval.start')
             return(sp)
           }
           )

setGeneric("plotDiff",
           function(object,name=NULL){
             standardGeneric("plotDiff")
           })

setMethod("plotDiff",
           signature=c(object="TransTrajectories"),
           definition=function(object,name=NULL){
             if(is.null(name)) stop('Please provide name of (sub)trajectory to plot.')
             data<-data.frame(object@ttrajlist[[name]]);colnames(data)<-c('V1','V2')
             ddata<-data.frame(object@difftraj[[name]]);colnames(ddata)<-c('V1','V2')
             
             p1<-ggplot(data,aes_string(x='V1',y='V2')) + 
                geom_point(alpha=0.3,size=3) + 
                ggtitle(paste(name,'(before differentiation)',sep=""))+
                scale_y_continuous("Free energy")+
                scale_x_continuous("Time")
             p2<-ggplot(ddata,aes_string(x='V1',y='V2')) + 
               geom_point(alpha=0.3,size=3) + 
               ggtitle(paste(name,'(after differentiation)',sep=""))+
               scale_y_continuous("Free energy")+
               scale_x_continuous("Time")
             multiplot(p1,p2,cols=1)
           })

setGeneric("splitTraj",
           function(object,segsplits=rep(5,length(object@filenames))){
             standardGeneric("splitTraj")
           })

setMethod("splitTraj",
          signature=c(object="Trajectories"),
          definition=function(object,segsplits=rep(5,length(object@filenames))){
            ntraj<-length(object@trajlist)
            tnames<-object@filenames
            if(class(segsplits)!="numeric") stop(paste("Please provide a numeric vector of points! Current input:",class(segsplits)))
            else if(length(segsplits)!=ntraj) stop(paste("Please provide number of splits equal to number of trajectories (",ntraj,")",sep=""))
            mm<-as.list(rep(0,ntraj))
            for(i in 1:ntraj){
              dat<-object@trajlist[[i]]
              mm_est_binseg<-cpt.mean(dat[,2],method="BinSeg",penalty="Asymptotic",pen.value=0.01,Q=segsplits[i],class=FALSE)
              mm[[i]]<-mm_est_binseg#sort(mm_est_binseg$cps[1,1:mm_est_binseg$op.cpts])

              names(mm)[i]<-tnames[i]
            }
            return(mm)
          })

setGeneric("transformSeries",
           function(object,method="splitting",breakpoints=1){
             standardGeneric("transformSeries")
           })

setMethod("transformSeries",
          signature=c(object="Trajectories",breakpoints="ANY"),
          definition=function(object,method='splitting',breakpoints=1){
            options(warn=-1)
            #check arguments
            if(!is.element(method,c('splitting','override_splitting','differentiation'))) stop(paste(method,'is not a valid input for argument "method".'))
            
            files<-getTraj(object)
            totalavd<-getAVD(object)
            filenames<-getTNames(object)
            
            if(method=="override_splitting" & length(files)!=length(breakpoints)) stop("parameter segsplits must be a list of breakpoints equal to # of trajectories.")
            else if(method=="splitting" & length(breakpoints)!=1) stop("parameter segsplits must be an integer specifying number of splits required.") 
            
            avd<-leg_files<-c()
            results<-difftraj<-list()
            
            for(i in 1:length(files)){
              dat<-files[[i]]
              
              if(method=="override_splitting"){
                
                splits1<-sort(c(breakpoints[[i]],nrow(dat)))
                k<-1
                
                for(j in 1:length(splits1)){
                  dat1<-dat[k:splits1[j],]#dat1<-dat[k:splits1[j],]
                  results[[length(results)+1]]<-dat1
                  avd<-c(avd,adf.test(dat1[,2])$p.value)
                  leg_files<-c(leg_files,paste(filenames[i],"_",j,sep=""))
                  k<-splits1[j]+1
                }
                
              }
              # no transformation done
              else if(totalavd[i]<=0.05){
                leg_files<-c(leg_files,filenames[i])
                results[[length(results)+1]]<-dat
                avd<-c(avd,totalavd[i])
              }
              
              else if(method=="differentiation"){
                results[[length(results)+1]]<-dat
                avd<-c(avd,totalavd[i])
                dat[2:nrow(dat),2]<-diff(dat[,2])
                dat1<-dat[-1,]
                difftraj<-c(difftraj,list(dat1))
                names(difftraj)[length(difftraj)]<-filenames[i]
                avd[i]<-adf.test(dat1[,2])$p.value
                leg_files<-c(leg_files,filenames[i])
              }
              
              else if(method=="splitting"){
                mm_est_binseg<-cpt.mean(dat[,2],method="BinSeg",penalty="Asymptotic",pen.value=0.01,Q=breakpoints,class=FALSE)
                mm<-mm_est_binseg#c(sort(mm_est_binseg$cps[1,1:mm_est_binseg$op.cpts]),nrow(dat))
                k<-1
                
                for(j in 1:length(mm)){
                  dat1<-dat[k:mm[j],]  	
                  results[[length(results)+1]]<-dat1
                  avd<-c(avd,adf.test(dat1[,2])$p.value)
                  leg_files<-c(leg_files,paste(filenames[i],"_",j,sep=""))
                  k<-mm[j]+1
                }
              }
            }
            names(results)<-leg_files
            tempTraj<-new("TransTrajectories",breakpoints=as.list(breakpoints),tmethod=method,ttrajlist=results,tfilenames=leg_files,tavd=avd,difftraj=difftraj)
            as(tempTraj,"Trajectories") <- as(object,"Trajectories")
            return(tempTraj)
          })
