
rm(list=ls())
library(Iso)

get.oc <- function(targetE,targetT, pE.true,pT.true,u1,u2, ncohort, cohortsize, cutoff.eliT=0.95, cutoff.eliE=0.90, nstar=9,ntrial=10)
{

safe = 1
stop=150
startdose=1

  peestimate<-function(yE,n){
    ndose<-length(yE)
    lik<-rep(0,ndose)
    pe<-(yE+0.05)/(n+0.1)
    p.e<-matrix(NA,ndose,ndose)
    for (i in 1:ndose){
      if (i==1) {x<-seq(ndose,1,by=-1)} else {x<-c(1:(i-1),seq(ndose,i))}
      #x<-x
      p.e[i,]<-ufit(pe,lmode=i,x=x,w=n+0.5)[[2]]
      lik[i]<-prod(dbinom(yE,n,p.e[i,]))		
    }
    lik<-lik/sum(lik)
    pe<-t(p.e)%*%lik+0.01*seq(1,ndose)
    return(pe)}
  
  ## to make get.oc as self-contained function, we copied functions get.boundary() and select.mtd() here.
  get.boundary <- function(target, ncohort, cohortsize=3, p.saf=NA, p.tox=NA,  cutoff.eli=0.95)
  {
    
    # if the user does not provide p.saf and p.tox, use the default values
    if(is.na(p.saf)) p.saf=0.6*target;
    if(is.na(p.tox)) p.tox=1.4*target;
    
    ### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
    npts = ncohort*cohortsize;
    ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;
    for(n in (1:ncohort)*cohortsize)
    {
      error.min=3;
      for(m1 in 0:(n-1))
      {
        for(m2 in (m1+1):n)
        {
          
          error1 = pbinom(m1, n, target)+1-pbinom(m2-1, n, target);
          error2 = 1-pbinom(m1, n, p.saf);
          error3 = pbinom(m2-1, n, p.tox);
          
          error=error1+error2+error3;
          if(error<error.min) {error.min=error; cutoff1=m1; cutoff2=m2;}
        }
      }
      ntrt = c(ntrt, n);
      b.e = c(b.e, cutoff1);
      b.d = c(b.d, cutoff2);
      
      elimineed=0; # indicating whether elimination is needed
      if(n<3) { elim = c(elim, NA); }  # require treating at least 3 patients before eliminating a dose
      else
      {
        for(ntox in 3:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
        {
          if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli) {elimineed=1; break;}
        }
        if(elimineed==1) { elim = c(elim, ntox); }
        else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs
      }
    }
    for(i in 1:length(b.d)) { if(!is.na(elim[i]) && (b.d[i]>elim[i])) b.d[i]=elim[i]; }
    boundaries = rbind(ntrt, elim, b.d, b.e);
    rownames(boundaries) = c("Number of patients treated", "Eliminate if # of DLT >=", 
                             "Deescalate if # of DLT >=",  "Escalate if # of DLT <=");
    colnames(boundaries) = rep("", ncohort);
    
    return(boundaries);
  }
  
  set.seed(32);
  ndose=length(pE.true)	
  npts = ncohort*cohortsize;
  YT=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
  YE=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
  N=matrix(rep(0, ndose*ntrial), ncol=ndose); # store the number of patients
  dselect = rep(0, ntrial); # store the selected dose level
  sel=rep(0,ndose);
  pts=rep(0,ndose);
  dlt=rep(0,ndose);
  eff=rep(0,ndose);
  ntox=0;
  neff=0;
  temp=get.boundary(targetT, ncohort, cohortsize,cutoff.eli=cutoff.eliT) 	
  b.e=temp[4,];   # escalation boundary
  b.d=temp[3,];   # deescalation boundary
  b.elim=temp[2,];  # elimination boundary
  u.true<-(u1*pE.true+(1-pT.true)*u2);
  bd=which.max(u.true*(pT.true<(targetT+0.1))*(pE.true>(targetE-0.05)));
  poorall=0;
  incoherent=0;
  tox.l<-seq(0.0,0.9,by=0.1)
  tox.u<-seq(0.1,1,by=0.1)
  u.l<-seq(0.0,1-0.1,by=0.1)
  u.u<-seq(0.1,1,by=0.1)
  ################## simulate trials ###################
  for(trial in 1:ntrial)
  {
    yT<-yE<-rep(0, ndose);    ## number of DLT at each dose level
    n<-rep(0, ndose);         ## number of patients treated at each dose level
    earlystop=0;              ## indiate whether the trial terminates early
    d=startdose;              ## starting dose level
    elimi = rep(0, ndose);    ## indicate whether doses are eliminated due to toxicity
    elimiE=  rep(0,ndose);    ## indicate whether doses are eliminated due to efficacy
    incoh=0;                  ## count incoherent movement
    
    posH<-rep((u1*0.5+u2)/100*length(u.l),ndose)
    #posH<-rep(10,ndose)
    safe<-0
    for(i in 1:ncohort)  
    {  			
      ### generate toxicity outcome
      #if(d<ndose&n[min(d+1,ndose)]>0){safe=1}
      wT = sum(runif(cohortsize)<pT.true[d])
      yT[d] = yT[d] + wT;
      wE = sum(runif(cohortsize)<pE.true[d])
      yE[d] = yE[d] + wE;
      n[d] = n[d] + cohortsize;
      nc = n[d]/cohortsize;
      if(n[d]>=nstar){
        uhat = (u1*yE[d]/n[d]+(1-yT[d]/n[d])*u2)/100*n[d] 
      } else {
        uhat = (u1*yE[d]/n[d]+u2)/100*n[d]
      }
      posH[d]<-which.max(pbeta(u.u,1+uhat,1+n[d]-uhat)-pbeta(u.l,1+uhat,1+n[d]-uhat))
      posH[d]<-posH[d]+(1-pbeta(u.u[posH[d]],1+uhat,1+n[d]-uhat))
      tox.ind<-which.max(pbeta(tox.u,1+yT[d],1+n[d]-yT[d])-pbeta(tox.l,1+yT[d],1+n[d]-yT[d]))
      # print(uhat)
      #print(posH)
      if(n[d]>=nstar){safe=1} else{safe=0}
      if(n[d]>=stop){break}
      if(!is.na(b.elim[nc]))
      {
        if(yT[d]>=b.elim[nc]) 
        {      
          elimi[d:ndose]=1;
          if(d==1) {earlystop=1; break;} 
        }
        
      }
      eff_cut<-pbeta(targetE,yE[d]+1,n[d]-yE[d]+1)
      if(eff_cut>cutoff.eliE) {elimi[d]=1;}
      posH <- posH*(1-elimi);
      
      if (tox.l[tox.ind]>(targetT+0.00001) && d!=1) {
        if(sum(elimi[1:(d-1)]==0)>0){d_opt=max(which(elimi[1:(d-1)]==0))} else {d_opt=d}
      } else if (tox.l[tox.ind]>(targetT+0.00001) && d==1) {if(elimi[d]==0){d_opt=d} else{earlystop=1;break}
      } else{
        admi_set=d;
        if(d>1){
          if(sum(elimi[1:(d-1)]==0)>0){admi_set<-c(admi_set,max(which(elimi[1:(d-1)]==0)))} 
        }
        if(d<ndose){
          if(safe==0){
            if(sum(elimi[(d+1):ndose]==0)>0){admi_set<-c(admi_set,d+min(which(elimi[(d+1):ndose]==0)))}
          } else {
            if(tox.l[tox.ind]<targetT & sum(elimi[(d+1):ndose]==0)>0){admi_set<-c(admi_set,d+min(which(elimi[(d+1):ndose]==0)))}
          }
        }
        #if(length(admi_set)>1 & eff_cut>0.5 & n[d]>=9){admi_set<-admi_set[-1]}
        temp.posH<-posH[admi_set]+runif(length(admi_set))*(10^-15)
        d_opt=admi_set[which.max(temp.posH)]
      }
      if (elimi[d_opt]==1) {earlystop=1; break} 
      if (sum(elimi)==ndose) {earlystop=1; break}
      if(((yT[d]/n[d])>targetT)&d_opt>d){incoh<-incoh+1}
      # print(c(d,d_opt))
      # print(yT)
      # print(yE)
      # print(n)
      d<-d_opt
    }
    incoherent = incoherent+(incoh/i)/ntrial*100
    YT[trial,]=yT;
    YE[trial,]=yE;
    N[trial,]=n;
    ntox=ntox+sum(yT)/ntrial
    neff=neff+sum(yE)/ntrial
    if (earlystop==0){
      pT<-(yT+0.05)/(n+0.1)
      pE<-(yE+0.05)/(n+0.1)
      pT<-pava(pT,n+0.1)+0.001*seq(1,ndose)
      pE<-peestimate(yE,n)
      u<-u1*pE+(1-pT)*u2
      u[elimi==1]<--100
      u[elimiE==1]<--100
      u[n==0]<--100
      #u[pT>(targetT+0.1)]<--100
      d_mtd<-which.min(abs(pT-targetT))
      d_opt<-which.max(u[1:d_mtd])	
      dselect[trial]=d_opt
      sel[dselect[trial]]=sel[dselect[trial]]+1/ntrial*100
    } else {dselect[trial]<-99}	
    pts<-pts+n/ntrial
    dlt<-dlt+yT/ntrial
    eff<-eff+yE/ntrial
    if(n[bd]<(npts/ndose)){poorall<-poorall+1/ntrial*100}
  }	
  pts<-round(pts,1)
  dlt<-round(dlt,1)
  bd.sel<-sel[bd]
  od.sel<-sum(sel[which(abs(u.true-max(u.true))<=3)])
  overdose<-sum(pts[pT.true>targetT])
  earlystop<-sum(dselect==99)/ntrial*100
  results=list(pT=pT.true,pE=pE.true,u=u.true,sel=sel,pts=pts,dlt=dlt,eff=eff,earlystop=earlystop,ntox=ntox,neff=neff,
               bd.sel=bd.sel,od.sel=od.sel,overdose=overdose,poorall=poorall)
  print(paste0("The selection percentages of all the doses are "))
  print(sel)
  print(paste0("The average number of patients allocated at all the doses are ",pts))
  print(pts)
  print(paste0("The average number of efficacy outcomes at all the doses are "))
  print(dlt)  
  print(paste0("The average number of toxicity outcomes at all the doses are "))
  print(eff)
  print(paste0("The percentage of early stopped trials is ",earlystop)  )
  print(paste0("The average total number of of toxicity outcomes is ",ntox))
  print(paste0("The average total number of of toxicity outcomes is ",neff))
  print(paste0("The selection percentage of the best dose is ",bd.sel))
  print(paste0("The selection percentage of the optimal dose is ",bd.sel))
  print(paste0("The number of overdoses is ",overdose))
  print(paste0("The percentage of poor allocation is ",poorall))
  
  return(results)	
}






cohortsize = 3
ncohort = 12

targetE=0.25
targetT=0.3
u1 = 70
u2 = 30


pT.true<-c(0.20,0.40,0.45,0.50,0.55)
pE.true<-c(0.40,0.50,0.60,0.70,0.80)
oc=get.oc(targetE,targetT, pE.true,pT.true,u1,u2, ncohort, cohortsize, cutoff.eliT=0.95, cutoff.eliE=0.90, nstar=9,ntrial=10)
print(oc)

  



