##First read in the arguments listed at the command line
#R CMD BATCH --no-save --no-restore '--args [model dir] [estimate dir + timestamp] [outdir] ' plot_recal.R /dev/stdout

args=(commandArgs(TRUE))



scoremax<-20;
qualmax<-50;

if(length(args)==0){
  print("No arguments supplied.")
}else{
  print(args);
  modelsdir<-paste(args[1],"/",sep="");
  estdir<-paste(args[2],"_",sep="");

                                        #counting cycles
  print(modelsdir);
  maxcycle<--1;
  for(cycle in 1:50000) {
    if(file.exists(paste(modelsdir,"SVMlight_recal_cycle_",cycle,"_1.par",sep=""))){
    }else{
      maxcycle<-cycle-1;
      break;
    }
  }

  if(maxcycle == -1){
    print("Unable to determine max cycles");
    quit();
  }
  print(maxcycle);

  cycleinfile<-50;

  for(cycle in 1:maxcycle) {
    if( (cycle %% cycleinfile) == 1){
      if( !(cycle == 1) ){
        dev.off();
      }
      pdf(paste(args[3],"/recalcurve",ceiling(cycle/cycleinfile),".pdf",sep=""),height=cycleinfile*2,width=4);
      mat<-matrix(c(seq(1,cycleinfile*4)), cycleinfile, 4, byrow = TRUE);
      layout(mat);
      par("mar"=c(4.5, 1, 4.5, 1)); #bottom, left, top and right margins
    }
    
    for(nuc in 1:4) {

      if(!file.exists(paste(modelsdir,"SVMlight_recal_cycle_",cycle,"_",nuc,".par",sep=""))){
        print(paste("Error file ",modelsdir,"SVMlight_recal_cycle_",cycle,"_",nuc,".par does not exist",sep=""));
        quit();
      }else{
        if(!file.exists(paste(estdir,"SVMlight_recal_cycle_",cycle,"_",nuc,".estp",sep=""))){
          print(paste("Error file ",estdir,"SVMlight_recal_cycle_",cycle,"_",nuc,".estp does not exist",sep=""));
          quit();          
        }else{
          
          datapar<-read.table(paste(modelsdir,"SVMlight_recal_cycle_",cycle,"_",nuc,".par",sep=""));
          dataest<-read.table(paste(estdir,"SVMlight_recal_cycle_",cycle,"_",nuc,".estp",   sep=""));

          par(xpd=FALSE);
          #plot(cycle,type="n",xlab="",ylab="",xlim=c(1,scoremax),ylim=c(1,qualmax));
          plot(dataest$V1,dataest$V4,xlab="",ylab="",xlim=c(1,scoremax),ylim=c(1,qualmax));
          abline(datapar$V1[2],datapar$V1[3],col = 'red');
          abline(datapar$V1[4],datapar$V1[5],col = 'red');
          segments(datapar$V1[1], -1.5, datapar$V1[1], datapar$V1[2] + datapar$V1[3]*datapar$V1[1],lty=2);
          par(xpd=NA);
          text(datapar$V1[1], -3,"k");
          if(nuc ==1 ){
            text(7,qualmax+5,paste("cycle:",cycle));
          }
        }
      }
    }
    

  }
  dev.off();

}


