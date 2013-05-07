setMethod(f="plot",signature="Trajectories",
          function(x,y,name='all',breakpoints=NULL,...){
            # make as generic plot function
            # takes object of class Trajectories with or without breakpoints
            filenames<-getTNames(x)
            avd<-getAVD(x)
            trajs<-getTraj(x)
            ntraj<-length(trajs)

            ll<-0
            all<-mm<-posx<-labelp<-c()
            if(!is.null(breakpoints)) mm_all<-c()
            
            if(name=='all'){
              for(i in 1:ntraj){
                dat<-data.frame(trajs[[i]])
                colnames(dat)<-c("Time","FreeEnergy")
                ll<-c(ll,(nrow(dat)+ll[i]))
                all<-rbind(all,dat)
                posx<-c(posx,(ll[(i+1)]+ll[i])/2)
                labelp<-c(labelp,paste("adf.test p-value = ",round(avd[i],3)))
                if(!is.null(breakpoints)){
                  mm<-breakpoints[[i]]
                  mm_all<-c(mm_all,(ll[(length(ll)-1)]+mm))
                }                
              }
              all<-data.frame(all)
              colnames(all)<-c("Time","FreeEnergy")
              all[,1]<-1:nrow(all)
              allavd<-adf.test(all[,2])$p.value
              p1<-ggplot(all,aes_string(x='Time',y='FreeEnergy')) + 
                geom_point(size = 3,alpha=0.3)+
                theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))+
                theme(axis.title.x=element_text(size=16))+theme(axis.title.y=element_text(size=16))+
                scale_x_continuous("snapshots extracted from replicates")
              
              option1<-option2<-NULL
              if(!is.null(breakpoints)) option2<- geom_vline(xintercept=mm_all,color="#3399FF")
              if(ntraj>1) option1<-geom_vline(xintercept=ll[2:(length(ll)-1)],color="#CC0000",linetype="dashed",size=1.2)
              
              posy<-c(rep(max(all[,2])+2,ntraj),rep(min(all[,2])-2,ntraj))
              posx<-rep(posx,2)
              plot(p1+option1+option2+annotate("text",x=posx,y=posy,label=c(filenames,labelp)))              
            }
            else{
              dat<-data.frame(trajs[[name]])
              colnames(dat)<-c("Time","FreeEnergy")
              if(!is.null(breakpoints)) mm<-breakpoints[[name]]
              p1<-ggplot(dat,aes_string(x='Time',y='FreeEnergy')) + 
                geom_point(alpha=0.3,size=3) + 
                theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))+
                theme(axis.title.x=element_text(size=16))+theme(axis.title.y=element_text(size=16))+
                ggtitle(paste(name,'\n',paste("adf.test p-value = ",round(avd[which(x@filenames==name)],3),sep="")))
              if(length(mm)>0) option1<-geom_vline(xintercept=mm,color="#3399FF")
              plot(p1+option1)
            }        
          })

setMethod(f="plot",signature="TransTrajectories",
          function(x,y,labelling=TRUE,...){
            filenames<-getTNames(x)
            avd<-getAVD(x)
            trajs<-getTraj(x)
            ntraj<-length(trajs)
            
            ll<-0
            all<-posx<-labelp<-c()
            
            for(i in 1:ntraj){
              dat<-trajs[[i]]
              all<-rbind(all,dat)
              ll<-c(ll,(nrow(dat)+ll[length(ll)]))
              posx<-c(posx,(ll[(i+1)]+ll[i])/2)
              labelp<-c(labelp,paste("adf.test p-value = ",round(avd[i],3)))
            }
            allavd<-adf.test(all[,2])$p.value
            all<-data.frame(all)
            all[,1]<-1:nrow(all)
            colnames(all)<-c("Time","FreeEnergy")
            
            option1<-option2<-NULL
            
            p1<-ggplot(all,aes_string(x='Time',y='FreeEnergy')) + 
              geom_point(size = 3,alpha=0.3) +
              theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))+
              theme(axis.title.x=element_text(size=16))+theme(axis.title.y=element_text(size=16))
            if(ntraj>1) option1<-geom_vline(xintercept=ll[2:(length(ll)-1)],color="#CC0000")
            
            if(labelling){
              posy<-c(rep(max(all[,2])+2,ntraj),rep(min(all[,2])-2,ntraj))
              posx<-rep(posx,2)       
              option2<-annotate("text",x=posx,y=posy,label=c(filenames,labelp),size=4)
            }
            plot(p1+option1+option2)
          })

setMethod(f="plot",signature="SegTrajectories",
          function(x,y,...){
            data<-getSegments(x)
            changepoints<-c(c(1,which(diff(data[,28])!=0)+1))
            ser<-unique(data[,28])
            aa<-data.frame(cbind(1:nrow(data),data[,c(1,6)]))
            colnames(aa)<-c("Time","FreeEnergy","Groups")
            aa[,3]<-as.character(aa[,3])
            ll<-paste("s",1:length(unique(data[,5])),sep="")
            cp<-c(0,unique(data[,5]))
            posx<-c()
            posy<-rep(max(aa[,2]+3),length(ll))
            for(i in 1:length(ll)) posx<-c(posx,mean(c(cp[i],cp[(i+1)])))
      
            p1<-ggplot(aa,aes_string(x='Time',y='FreeEnergy',color='Groups'))+geom_point()+scale_color_hue()+theme(legend.position="none")+
              theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))+
              theme(axis.title.x=element_text(size=16))+theme(axis.title.y=element_text(size=16))+
              scale_y_continuous("Free energy")+
              scale_x_continuous("snapshots extracted from replicates")+
              geom_vline(xintercept=changepoints[2:length(changepoints)],linetype = "longdash",color="darkgray")
            
            plot(p1+annotate("text",x=posx,y=posy,label=ll,size=4))
          })

setMethod(f="plot",signature="SegSeriesTrajectories",
          function(x,y,...){
            data<-getSegments(x)
            changepoints<-c(c(1,which(diff(as.numeric(as.factor(data$seriesID)))!=0)+1))
            aa<-data.frame(cbind(1:nrow(data),data$observed,data$subpopulation))
            colnames(aa)<-c("Time","FreeEnergy","Subpopulations")
            aa[,3]<-as.character(aa[,3])
            p1<-ggplot(aa,aes_string(x='Time',y='FreeEnergy',color='Subpopulations')) +
              geom_point(aes_string(shape='Subpopulations'),size=3)+scale_color_grey(end=0.7)+
              ggtitle("Identified subpopulations")+
              theme(legend.title=element_text(size=15),legend.text=element_text(size=12),plot.title=element_text(size=20))+
              theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))+
              theme(axis.title.x=element_text(size=16))+theme(axis.title.y=element_text(size=16))+
              scale_x_continuous("snapshots extracted from replicates")+
              geom_vline(xintercept=changepoints[2:length(changepoints)],linetype = "longdash",color="darkgray")
            plot(p1)
          })
