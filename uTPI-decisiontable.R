get.decision.table<-function(u1,u2,targetE,targetT,cutoff.eliE,cutoff.eliT,nstar=9){

tox.l<-seq(0.0,0.9,by=0.1)
tox.u<-seq(0.1,1,by=0.1)
u.l<-seq(0.0,1-0.1,by=0.1)
u.u<-seq(0.1,1,by=0.1)
 
temp.tab<-NULL

for(n in c(0,3,6,9)){
  for(yT in 0:n){
    if(pbeta(targetT,yT+1,n-yT+1)<(1-cutoff.eliT)){
      tox.ind<-which.max(pbeta(tox.u,1+yT,1+n-yT)-pbeta(tox.l,1+yT,1+n-yT))
      temp.tab=rbind(temp.tab,c(n,paste(">=",yT),">=0",tox.ind,0));
      break;
    }
    
    for(yE in 0:n){
      if(pbeta(targetE,yE+1,n-yE+1)>cutoff.eliE)
      {
        tox.ind<-which.max(pbeta(tox.u,1+yT,1+n-yT)-pbeta(tox.l,1+yT,1+n-yT))
        temp.tab=rbind(temp.tab,c(n,yT,paste("<=",yE),tox.ind,0));next}
      if(n>=nstar){
        uhat = (u1*yE/n+(1-yT/n)*u2)/100*n 
      } else {
        uhat = (u1*yE/n+u2)/100*n
      }
      posH<-which.max(pbeta(u.u,1+uhat,1+n-uhat)-pbeta(u.l,1+uhat,1+n-uhat))
      posH<-posH+(1-pbeta(u.u[posH],1+uhat,1+n-uhat))
      tox.ind<-which.max(pbeta(tox.u,1+yT,1+n-yT)-pbeta(tox.l,1+yT,1+n-yT))
      if(n==0){posH=(u1*0.5+u2)/100*length(u.l); tox.ind=0}
      temp.tab=rbind(temp.tab,c(n,yT,yE,tox.ind,posH));
    }
  }
}





aa<-as.numeric(noquote(temp.tab[,5]))
temp.tab<-cbind(temp.tab,pmax(rank(aa)-sum(aa==0),0))
temp.tab<-temp.tab[,-5]
colnames(temp.tab)<-c("Num. Patients","Num. Tox","Num. Eff","Tox Interval","Desirability Score")
temp.tab<-noquote(temp.tab)
temp.tab

return(temp.tab)
}


u1<-70
u2<-30
targetE = 0.25
targetT =0.30
cutoff.eliE = 0.90
cutoff.eliT = 0.95
get.decision.table(u1,u2,targetE,targetT,cutoff.eliE,cutoff.eliT)