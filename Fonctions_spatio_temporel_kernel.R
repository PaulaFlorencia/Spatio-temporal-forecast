library(evd)
library(latex2exp)
library(mgcv)
library(gamlss)
##################################################################################
################## Main fonctoin used in file #####################################
DataStructure.func <- function(Tnext = 86, S = 8, dependence.type = "log",
                               dependence.param = 0.3, xi = -0.2, support = 50,
                               phi.mat= matrix(replicate(8,seq(20, 35, length.out = 86)),ncol=8)){
  # objective : simulate max-stable dependent trajectories with same xi and support
  # input: Distribution parameters, dependance stucture, number of sites and time steps
  # output : list of 2 elements: 1) matrix containing trajectories between t=1 and T = Tnext-1
  #                              2) vector of values at Tnext
  
  #############################################################################
  
  # code: 
  
  # Generate stationary maxstable data with defined dependance structure
  dat_maxsable <- evd::rmvevd(Tnext, dep= dependence.param, model = c(dependence.type), d = S, mar=c(0, 1,xi)) # matrix
  
  # empty matrix
  dat_gev <- matrix(NA, nrow=dim(dat_maxsable)[1], ncol=dim(dat_maxsable)[2]) 
  
  # Add trend to our times series
  if (xi < 0){
    sig.mat <- -xi*(support-phi.mat)
    for ( i in 1:dim(dat_maxsable)[2]){
      dat_gev[,i] <- dat_maxsable[,i]*sig.mat[,i] + phi.mat[,i]}
  }
  if (xi > 0){
    mu.mat <- support + phi.mat/xi
    for ( i in 1:dim(dat_maxsable)[2]){
      dat_gev[,i] <- dat_maxsable[,i]*phi.mat[,i] + mu.mat[,i]}
  }
  
  return(list("SimulatedObservations"=dat_gev[1:Tnext-1,],"SimulatedNext"=dat_gev[Tnext,]))
}

TruePnext <- function(sigma.mat, xi, alpha, dependence.type = "log",siteid = 0){
  # Objective : compute true records probability at Tnext, from theoretical parameters
  # input : sigma.mat (matrix of all scale parameters), xi(shape parameter) and type of dependence and it's parameter
  #         sitsid is a value or a vector with the "locations of interest" (where we want to know records probability)
  #         if siteid = 0, all locations are locations of interest
  # output : list of 5 arrays.
  #          for the 3-dimension arrays we have 1d = time (rows) , 2d=sites of comparison (cols), 2d=sites of interest
  
  #############################################################################
  # code: 
  
  Tobs <- (dim(sigma.mat)[1] -1)
  nbsites <- dim(sigma.mat)[2]
  sigma.obs <- sigma.mat[1:Tobs,]
  sigma.next <- sigma.mat[Tobs+1,]
  sigma.obs_s <- t(t(apply(sigma.obs^(1/xi), 2,sum)))^xi 
  
  if (siteid[1] == 0){ # site id = 0 means all sites
    
    #arrays
    interest_sites <- c(1:nbsites)
    P.next <- array(NA,dim = c(1,1,nbsites))
    lambda.next <- array(NA,dim = c(1,1,nbsites))
    PMarginal.next <- array(NA,dim = c(1,nbsites,nbsites))
    lambdaMarginal.next <- array(NA,dim = c(1,nbsites,nbsites))
    
    # reste
    for (i in 1:nbsites){
      lambda.next_s <- (sigma.obs_s/sigma.next[i])^(1/xi)
      lambdaMarginal.next[1,,i] <- lambda.next_s
      PMarginal.next[1,,i] <- 1/(1+lambda.next_s)
      lambda.next_s_inv_sansV <- lambda.next_s^(-1)
      if (dependence.type == "log"){
        lambda.next_s_inv_sansV <- t(as.numeric(lambda.next_s_inv_sansV))
        lambda.nextV <- V(lambda.next_s_inv_sansV,alpha)
        lambda.next[,,i] <- lambda.nextV 
        Precord <- 1/(1+lambda.nextV)
        P.next[,,i] <- Precord
      }
    }
  }
  
  if (siteid[1] != 0){
    sites.so_length <- length(siteid)
    
    #arrays
    interest_sites <- siteid
    P.next <- array(NA,dim = c(1,1,sites.so_length))
    lambda.next <- array(NA,dim = c(1,1,sites.so_length))
    PMarginal.next <- array(NA,dim = c(1,nbsites,sites.so_length))
    lambdaMarginal.next <- array(NA,dim = c(1,nbsites,sites.so_length))
    
    # reste
    for (i in 1:length(siteid)){
      site <- siteid[i]
      lambda.next_s <- (sigma.obs_s/sigma.next[site])^(1/xi)
      lambdaMarginal.next[1,,i] <- lambda.next_s
      PMarginal.next[1,,i] <- 1/(1+lambda.next_s)
      lambda.next_s_inv_sansV <- lambda.next_s^(-1)
      if (dependence.type == "log"){
        lambda.next_s_inv_sansV <- t(as.numeric(lambda.next_s_inv_sansV))
        lambda.nextV <- V(lambda.next_s_inv_sansV,alpha)
        lambda.next[,,i] <- lambda.nextV
        Precord <- 1/(1+lambda.nextV)
      }
      P.next[,,i] <- Precord
    }
  }
  return(list("index_SitesInteret" = interest_sites,
              "ProbabiliteMarginal.next_SitesInteret"=PMarginal.next,
              "lambdaMarginal.next_SitesInteret"=lambdaMarginal.next,
              "ProbabiliteTotal.next_SitesInteret"=P.next,
              "lambdaTotal.next_SitesInteret"=lambda.next))
} 

TrueP1toT <- function(sigma.mat, sigma.next, xi, alpha){
  TruePrecord <- matrix(NA, ncol=(dim(sigma.mat)[2]+1), nrow=dim(sigma.mat)[1])
  TrueLambdas <- matrix(NA, ncol=(dim(sigma.mat)[2]+1), nrow=dim(sigma.mat)[1])
  for (i in 1: dim(sigma.mat)[2]){
    sigma_maxst <- cumsum((sigma.mat[,i]^(1/xi)))^xi
    lambda_st <- (sigma_maxst/ c(sigma.mat[-1,1],sigma_next))^(1/xi)
    P_record <- 1/(1+lambda_st)
    TruePrecord[,i] <- P_record
    TrueLambdas[,i] <- lambda_st
  }
  lambdaMatMatinv <- TrueLambdas^(-1)
  lambda_Vt <- V(lambdaMatMatinv[,1:dim(sigma.mat)[2]],alpha)
  Truerecord_all <- 1/(1+lambda_Vt)
  TrueLambdas[,(dim(sigma.mat)[2]+1)] <- lambda_Vt
  TruePrecord[,(dim(sigma.mat)[2]+1)] <- Truerecord_all
  return(list(TruePrecord,TrueLambdas))
}

TrueP2toTnext <- function(sigma.mat, xi, alpha, dependence.type = "log", siteid = 0){
  # Objective : compute true records probability at from t=2 to t=Ttoday=Tnext-1 from theoretical parameters
  # input : sigma.mat (matrix of all scale parameters), xi(shape parameter) and type of dependence and it's parameter
  #         sitsid is a value or a vector with the "locations of interest" (where we want to know records probability)
  #         if siteid = 0, all locations are locations of interest
  # output : list of 5 arrays.
  #          for the 3-dimension arrays we have 1d = time (rows) , 2d=sites of comparison (cols), 2d=sites of interest
  
  #############################################################################
  # code: 
  
  nbsites <- dim(sigma.mat)[2]
  Tnext <- dim(sigma.mat)[1]
  Tobs <- Tnext - 1
  sigma.obs <- sigma.mat[1:Tobs,]
  sigma.next <- as.numeric(sigma.mat[Tnext,])
  sigmaMax.obs <- (apply(sigma.obs^(1/xi),2,cumsum))^xi # ncol = nbsites , nrow = Tobs
  
  if (siteid[1] == 0){
    
    # vecteurs a retourner
    interest_sites <- c(1:nbsites)
    # a rendre dans une liste finale
    lambda.2toTnext_marginal <- array(NA,dim = c(Tobs,nbsites,nbsites))
    # a rendre dans une liste finale
    P.2toTnext_marginal <- array(NA,dim = c(Tobs,nbsites,nbsites))
    # a rendre dans une liste finale
    lambda.2toTnext <- matrix(NA,nrow=Tobs, ncol=nbsites)
    # a rendre dans une liste finale
    P.2toTnext <- matrix(NA,nrow=Tobs, ncol=nbsites)
    
    for (i in 1:nbsites){
      lambda.2toTnext_marginal_indx <- (sigmaMax.obs/ sigma.mat[-1, i])^(1/xi) # matrix a sauvgarder
      P.2toTnext_marginal_indx <- 1/(1+lambda.2toTnext_marginal_indx) # matrix a sauvgarder
      
      lambda.2toTnext_marginal[,,i]<- lambda.2toTnext_marginal_indx
      P.2toTnext_marginal[,,i]<- P.2toTnext_marginal_indx 
      
      if (dependence.type == "log"){
        lambdaInv.2toTnext_marginal_indx <- lambda.2toTnext_marginal_indx^(-1) # <- faire un array de vecterys 
        lambda.2toTnext_indx <- V(lambdaInv.2toTnext_marginal_indx,alpha) # <- mettre dans la liste
        P.2toTnext_indx <- 1/(1+lambda.2toTnext_indx) # un vecteur de la pro de t=2 a t=Tnext <- faire un array de vecterys 
        
        lambda.2toTnext[,i] <- lambda.2toTnext_indx
        P.2toTnext[,i] <- P.2toTnext_indx 
      }
    }
  }
  
  if (siteid[1] != 0){
    sites.so_length <- length(siteid)
    
    interest_sites <- siteid
    # a rendre dans une liste finale
    lambda.2toTnext_marginal <- array(NA,dim = c(p_timesteps,nbsites,sites.so_length))
    # a rendre dans une liste finale
    P.2toTnext_marginal <- array(NA,dim = c(p_timesteps,nbsites,sites.so_length))
    # a rendre dans une liste finale
    lambda.2toTnext <- array(NA,dim = c(p_timesteps,1,sites.so_length))
    # a rendre dans une liste finale
    P.2toTnext <- array(NA,dim = c(p_timesteps,1,sites.so_length))
    
    for (i in 1:sites.so_length){
      site <- siteid[i]
      lambda.2toTnext_marginal_indx <- (sigmaMax.obs/ sigma.mat[-1, site])^(1/xi) # matrix a sauvgarder
      P.2toTnext_marginal_indx <- 1/(1+lambda.2toTnext_marginal_indx) # matrix a sauvgarder
      lambda.2toTnext_marginal[,,i]<- lambda.2toTnext_marginal_indx
      P.2toTnext_marginal[,,i]<- P.2toTnext_marginal_indx 
      if (dependence.type == "log"){
        lambdaInv.2toTnext_marginal_indx <- lambda.2toTnext_marginal_indx^(-1) # <- faire un array de vecterys 
        lambda.2toTnext_indx <- V(lambdaInv.2toTnext_marginal_indx,alpha) # <- mettre dans la liste
        P.2toTnext_indx <- 1/(1+lambda.2toTnext_indx) # un vecteur de la pro de t=2 a t=Tnext <- faire un array de vecterys 
        lambda.2toTnext[,1,i] <- lambda.2toTnext_indx
        P.2toTnext[,1,i] <- P.2toTnext_indx 
      }
    }
  }
  
  return(list("index_SitesInteret" = interest_sites,
              "lambdaMarginal.2toTnext_SitesInteret"=lambda.2toTnext_marginal,
              "ProbabiliteMarginal.2toTnext_SitesInteret"=P.2toTnext_marginal,
              "lambdaTotal.2toTnext_SitesInteret"=lambda.2toTnext,
              "ProbabiliteTotal.2toTnext_SitesInteret"=P.2toTnext))
} 


Plot.ProbLambda <- function(ArrayProbLambda, siteid=1, site.other = 0, plot.type = "total",
                            par1=1,par2=2,div.screen="TRUE"){
  # Objective : Plot true theoretical marginal and aggregated record probabilities
  # Input : Array coming from function TrueP2toTnext()
  # Output: plots of marginal Prob & lambda or Total Prob & lambda
  # !!! this 
  
  #############################################################################
  # code:
  
  # siteid must be 1 value
  if(div.screen=="TRUE"){par(mfrow=c(par1,par2), mar = c(5, 5.3, 2, 2))}
  
  if(plot.type == "total"){
    
    if(siteid[1]==0){
      siteid <- as.numeric(ArrayProbLambda$index_SitesInteret)
    }
    
    if(length(siteid) == 1){
      proba_total.vec <- as.numeric(ArrayProbLambda$ProbabiliteTotal.2toTnext_SitesInteret[,siteid])
      lambda_total.vec <- as.numeric(ArrayProbLambda$lambdaTotal.2toTnext_SitesInteret[,siteid])
      x_axis <- c(1:length(proba_total.vec))
      plot(x_axis, proba_total.vec, type = "l", col='black',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s_o}$)'),
           main=bquote("Record prob at site "~ .(siteid)),lwd=2)
      plot(x_axis, lambda_total.vec, type = "l", col='black',xlab="t from 2 to T+1",ylab= TeX(r'($\lambda_t^{s_o}$)'),
           main=bquote("Associated" ~ lambda ~ "at site "~ .(siteid)),lwd=2)
    }
    
    if(length(siteid) > 1){
      proba_total.mat <- ArrayProbLambda$ProbabiliteTotal.2toTnext_SitesInteret[,siteid]
      lambda_total.mat <- ArrayProbLambda$lambdaTotal.2toTnext_SitesInteret[,siteid]
      
      x_axis <- c(1:dim(proba_total.mat)[1])
      
      for( i in 1:length(siteid)){
        iplot <- siteid[i]
        plot(x_axis, proba_total.mat[,i], type = "l", col='black',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s_o}$)'),
             main=bquote("Record prob at site "~ .(iplot)),lwd=2)
        plot(x_axis, lambda_total.mat[,i], type = "l", col='black',xlab="t from 2 to T+1",ylab= TeX(r'($\lambda_t^{s_o}$)'),
             main=bquote("Associated" ~ lambda ~ "at site "~ .(iplot)),lwd=2)
      }
    }
  }
  
  
  if(plot.type == "marginal"){
    site.other <- as.numeric(site.other)
    if(site.other[1]==0){
      site.other <- as.numeric(ArrayProbLambda$index_SitesInteret)
    }
    
    if(siteid[1]==0){
      siteid <- as.numeric(ArrayProbLambda$index_SitesInteret)
    }
    
    if(length(site.other) == 1){
      
      if(length(siteid) == 1){
        proba_marginal.vec <- as.numeric(ArrayProbLambda$ProbabiliteMarginal.2toTnext_SitesInteret[,site.other,siteid])
        lambda_marginal.vec <- as.numeric(ArrayProbLambda$lambdaMarginal.2toTnext_SitesInteret[,site.other,siteid])
        x_axis <- c(1:length(proba_marginal.vec))
        plot(x_axis, proba_marginal.vec, type = "l", col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s,s_o}$)'),
             main=bquote("Marginal record prob at site "~ .(siteid)~" w.r.t site " ~ .(site.other)),lwd=2)
        plot(x_axis, lambda_marginal.vec, type = "l", col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($\lambda_t^{s,s_o}$)'),
             main=bquote("Associated marginal" ~ lambda ~ "at site "~ .(siteid)~" w.r.t site " ~ .(site.other)),lwd=2)
      }
      if(length(siteid) > 1){
        for( i in 1:length(siteid)){
          proba_marginal.mat <- ArrayProbLambda$ProbabiliteMarginal.2toTnext_SitesInteret[,site.other,siteid]
          lambda_marginal.mat <- ArrayProbLambda$lambdaMarginal.2toTnext_SitesInteret[,site.other,siteid]
          x_axis <- c(1:dim(proba_marginal.mat)[1])
          iplot <- siteid[i]
          plot(x_axis, proba_marginal.mat[,i], type = "l", col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s,s_o}$)'),
               main=bquote("Marginal record prob at site "~ .(iplot)~" w.r.t site " ~ .(site.other)),lwd=2)
          plot(x_axis, lambda_marginal.mat[,i], type = "l", col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($\lambda_t^{s,s_o}$)'),
               main=bquote("Associated marginal " ~ lambda ~ "at site "~ .(iplot)~" w.r.t site " ~ .(site.other)),lwd=2)
        }
      }
    }
    
    
    if(length(site.other) > 1){
      
      if(div.screen=="TRUE"){par(mfrow=c(par1,par2), mar = c(5, 5.3, 4, 2))}
      
      if(length(siteid) == 1){
        
        proba_marginal.mat <- ArrayProbLambda$ProbabiliteMarginal.2toTnext_SitesInteret[,site.other,siteid]
        lambda_marginal.mat <- ArrayProbLambda$lambdaMarginal.2toTnext_SitesInteret[,site.other,siteid]
        x_axis <- c(1:dim(ArrayProbLambda$ProbabiliteMarginal.2toTnext_SitesInteret)[1])
        for( j in 1:length(site.other)){
          jplot <- site.other[j]
          plot(x_axis,proba_marginal.mat[,j],type='l', col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s,s_o}$)'),
               main=bquote(atop("Marginal Record prob at site " ~ .(siteid), " w.r.t site " ~ .(jplot))))
          plot(x_axis,lambda_marginal.mat[,j],type='l', col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($\lambda_t^{s,s_o}$)'),
               main=bquote(atop("Associated marginal " ~ lambda ~ "at site "~ .(siteid)," w.r.t site " ~ .(jplot))))
        }
      }
      
      if(length(siteid) > 1){
        proba_marginal.array <- ArrayProbLambda$ProbabiliteMarginal.2toTnext_SitesInteret[,site.other,siteid]
        lambda_marginal.array <- ArrayProbLambda$lambdaMarginal.2toTnext_SitesInteret[,site.other,siteid]
        x_axis <- c(1:dim(ArrayProbLambda$ProbabiliteMarginal.2toTnext_SitesInteret)[1])
        for (i in 1:length(siteid)){
          iplot <- siteid[i]
          for (j in 1:length(site.other)){
            jplot <- site.other[j]
            plot(x_axis,proba_marginal.array[,j],type='l', col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s,s_o}$)'),
                 main=bquote(atop("Marginal Record prob at site "~ .(iplot)," w.r.t site " ~ .(jplot))))
            plot(x_axis,proba_marginal.array[,j],type='l', col='blue',xlab="t from 2 to T+1",ylab= TeX(r'($P_t^{s,s_o}$)'),
                 main=bquote(atop("Associated marginal "~ lambda ~ "at site "~ .(iplot)," w.r.t site " ~ .(jplot))))
          }
        }
      }
    }
  }
  
  
  
}

PlotData.separate <- function(Data, siteid = 0, par1 = 4, par2 = 2 , add.max=FALSE , max.type="all" ){ # if  siteid = 0 we plot data in all sites 
  par(mfrow=c(par1,par2))
  timesteps <- dim(Data)[1]
  nbsites <- dim(Data)[2]
  if (length(siteid)==1){
    if (siteid>0){
      par1 <- 1
      par2 <- 1 
      plot(c(1:timesteps),Data[,siteid],ylim=c(min(Data),max(Data)),
           xlab="time",ylab='observations', main=bquote("Observations at site "~ .(siteid)))
      if(add.max==TRUE){
        max.traj <- apply(X=Data, MARGIN=1, FUN=max)
        points(c(1:timesteps),max.traj,col='red',pch=4)
        legend('bottomright',legend="max all sites",col="red",pch=4, bg="transparent")
      }
      #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1, cex=1.3,col="blue",paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
    }
    
    
    
    
    
    if (siteid==0){
      ylim.max = max(Data)
      ylim.min = min(Data)
      if(add.max==FALSE){
        for(i in 1:nbsites){
          plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="time",ylab='observations',
               main=bquote("Observations at site "~ .(i)))
          #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1,cex=1.3,col="blue", paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
        }
      }
      
      if(add.max==TRUE){
        max.traj <- apply(X=Data, MARGIN=1, FUN=max)
        for(i in 1:nbsites){
          plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="time",ylab='observations',
               main=bquote("Observations at site "~ .(i)))
          points(c(1:timesteps),max.traj,col='red',pch=4)
          #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1,cex=1.3,col="blue", paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
          legend('bottomright',legend="max all sites",col="red",pch=4, bg="transparent")
        }
      }
    }
  }
  
  if (length(siteid)>1){
    
    nb.plots <- length(siteid)
    ylim.max = max(Data)
    ylim.min = min(Data)
    
    if(add.max==FALSE){
      for(i in siteid ){
        plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="time",ylab='observations',
             main=bquote("Observations at site "~ .(i)))
        #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1, cex=1.3,col='blue',paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
      }
    }
    
    if(add.max==TRUE){
      if(max.type=="all"){
        max.traj <- apply(X=Data, MARGIN=1, FUN=max)
        colmax <- "red"
        legend.text <- "max all sites" }
      
      if(max.type=="not.all"){
        max.traj <- apply(X=Data[,c(siteid)], MARGIN=1, FUN=max)
        colmax <- "green"
        legend.text <- "max showed sites"}
      
      for(i in siteid ){
        plot(c(1:timesteps),Data[,i], ylim = c(ylim.min,ylim.max),xlab="time",ylab='observations',
             main=bquote("Observations at site "~ .(i)))
        points(c(1:timesteps),max.traj,col=colmax,pch=4)
        #text(c(1:timesteps)[10],max(Data[1:(timesteps/3)])-1, cex=1.3,col='blue',paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
        legend('bottomright',legend=legend.text,col=colmax,pch=4, bg="transparent")
      }
    }
    
  }
}

PlotData.compare <- function(Data, site.interest = 1, site.others = 0,  max.plot = FALSE){
  
  if(length(site.interest) != 1){
    stop("We must observe only 1 site of interest")}
  
  nbtimesteps <- dim(Data)[1]
  timesteps <- c(1:nbtimesteps)
  nbsites <- dim(Data)[2]
  
  site.id <- site.interest
  others.id <- site.others
  col.points.max <- "green"
  if(length(site.others) == 1){
    if(site.others == 0){
      others.id <- c(1:nbsites)
      col.points.max <- "red"}
  }
  
  if(max.plot == FALSE){par(mfrow=c(1,1), mar = c(5, 5.3, 4, 2))}
  if(max.plot == TRUE){par(mfrow=c(1,2), mar = c(5, 5.3, 4, 2))}
  
  ylim.max = max(Data)
  ylim.min = min(Data)
  
  sitesother.names <- paste(others.id,collapse=" ")
  
  plot(timesteps,Data[,site.id],ylim=c(ylim.min,ylim.max),xlab="time",ylab='observations',
       main=bquote(atop("Observations at site of interest "~ .(site.id),"w.r.t to other sites " ~ .(sitesother.names))))
  for(i in others.id ){points(timesteps,Data[,i],col='gray')}
  points(timesteps,Data[,site.id],col="black")
  # text(c(1:nbtimesteps)[13],max(Data)-2, cex=1.3,col='blue',paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
  legend('bottomright',legend=c("site interest","other sites"),col=c("black",'gray'),pch=1, bg="transparent")
  
  if(max.plot == TRUE){
    plot(timesteps,Data[,site.id],ylim=c(ylim.min,ylim.max),xlab="time",ylab='observations',
         main=bquote(atop("Observations at site of interest "~ .(site.id),"w.r.t to other sites " ~ .(sitesother.names))))
    for(i in others.id ){points(timesteps,Data[,i],col='gray')}
    points(timesteps,Data[,site.id],col="black")
    max.traj <- apply(X=Data[,c(site.id,others.id)], MARGIN=1, FUN=max)
    points(timesteps,max.traj,pch=4,col=col.points.max)
    # text(c(1:nbtimesteps)[13],max(Data)-2, cex=1.3,col='blue',paste0("xi = ", xi,"\n"," alpha = ",dep.para,"\n"))
    legend('bottomright',legend=c("site interest","other sites","max plotted sites"),col=c("black",'gray',col.points.max),pch=c(1,1,4), bg="transparent")
  }
}

GmaxY.estimation.function.step1 <- function(Data, site.interest = 0, site.others = 0, bandwidth= 35, plot=FALSE,
                                            plot.site.interest=1,plot.site.others=1){
  nbtimesteps <- dim(Data)[1]
  timesteps <- seq(1,nbtimesteps,by=1)
  nbsites <- dim(Data)[2]
  
  if(length(site.interest) > nbsites | length(site.others) > nbsites ){
    stop("the number of sites of study exceeds the number of columns of data ")}
  
  if(site.interest[1] == 0){site.interest <- seq(1,dim(Data)[2],by=1)}
  
  if(site.others[1] == 0){site.others <- seq(1,dim(Data)[2],by=1)}
  
  nbsites.interest <- length(site.interest)
  nbsites.others <- length(site.others)
  
  if (length(bandwidth) == 1){bandwidth <- rep(bandwidth,nbsites.others)}
  
  if(length(bandwidth)!= nbsites.others){
    stop("Bandwight must be a value or a vector of length nbsites.others")}
  
  array.Gmaxhat.mat <- array(NA,c((dim(Data)[1]-1),nbsites.others,nbsites.interest))
  
  for (i.dim3 in 1:nbsites.interest){
    i.interest <- site.interest[i.dim3]
    
    Gmaxhat.mat <- matrix(NA,ncol=nbsites.others,nrow=dim(Data)[1])
    
    for (i.dim2 in 1:nbsites.others){
      i.others <- site.others[i.dim2]
      Ghat <- G_estimator(timesteps, Data[,i.others], bandwidth[i.dim2])
      
      for (i.dim1 in 2:nbtimesteps){
        Gmaxhat_ti.func <- Gmax_glissant_estimator(t = i.dim1,Ghat)
        Gmaxhat <- Gmaxhat_ti.func(Data[i.dim1,i.interest]) 
        Gmaxhat.mat[i.dim1,i.dim2] <- Gmaxhat
      }
    }
    Gmaxhat.mat <-as.matrix(Gmaxhat.mat[-1,])
    array.Gmaxhat.mat[,,i.dim3] <- Gmaxhat.mat
  }
  
  id.1.plot <- which(as.numeric(site.interest) == plot.site.interest)
  id.2.plot <- which(as.numeric(site.others) == plot.site.others)
  
  if (plot==TRUE){
    par(mfrow=c(1,1), mar = c(5, 5.3, 4, 2))
    plot(timesteps[-1],array.Gmaxhat.mat[,id.2.plot,id.1.plot],xlab='t from 2 to T',
         ylab=TeX(r'($\hat{G}_{\max,t -1}^{s}(Y_{t}^{s_o})$)'),
         main=bquote("GmaxHat_s(Y_so) with so =  "~ .(id.1.plot) ~ "and s = " ~ .(id.2.plot)),col='red')
  }
  
  return(list("Array.GmaxHatY" = array.Gmaxhat.mat))
}

plot.separate.function.step.1 <- function(Array.step1, third.dimension = 1, second.dimension=0, par1=1, par2=1){
  par(mfrow=c(par1,par2), mar = c(5, 5.3, 4, 2))
  timesteps <- dim(Array.step1)[1]
  nbsites.interest <- dim(Array.step1)[3]
  nbsites.others <- dim(Array.step1)[2]
  
  if(second.dimension[1] == 0){second.dimension <- seq(1,nbsites.others,by=1)}
  
  if(length(second.dimension)==1){
    plot(c(1:timesteps),Array.step1[,second.dimension,third.dimension],ylim=c(min(Array.step1[,,third.dimension]),max(Array.step1[,,third.dimension])),
         xlab="t from 2 to T",ylab=TeX(r'($\hat{G}_{\max,t -1}^{s}(Y_{t}^{s_o})$)'), 
         main=bquote("GmaxHatY at 3d :"~ .(third.dimension) ~ " and 2d :"~ .(second.dimension)),col='red')
  }
  
  if(length(second.dimension)>1){
    for(i in 1:length(second.dimension)){
      plot(c(1:timesteps),Array.step1[,second.dimension[i],third.dimension],ylim=c(min(Array.step1[,,third.dimension]),max(Array.step1[,,third.dimension])),
           xlab="t from 2 to T",ylab=TeX(r'($\hat{G}_{\max,t -1}^{s}(Y_{t}^{s_o})$)'), 
           main=bquote("GmaxHatY at 3d :"~ .(third.dimension) ~ " and 2d :"~ .(second.dimension[i])),col="red")
    }
  }
}

HatEGmaxs <- function(h2, t, Gmaxhat.mat, Gmax.mat){
  EGmaxhat.mat <- matrix(NA, ncol=dim(Gmaxhat.mat)[2],nrow=dim(Gmaxhat.mat)[1])
  EGmax.mat <- matrix(NA, ncol=dim(Gmax.mat)[2],nrow=dim(Gmax.mat)[1])
  for (i in 1:dim(Gmaxhat.mat)[2]){
    Kij <- outer(t[-1],t[-1],function(zz,z) dEpan((zz - z) / h2[i])) 
    W <- Kij / rowSums(Kij)
    EGmax.mat[,i] <- W %*% as.numeric(Gmax.mat[,i])
    EGmaxhat.mat[,i] <- W %*% as.numeric(Gmaxhat.mat[,i])
  }
  return(list(EGmaxhat.mat,EGmax.mat))
}

EGmaxY.estimation.function.step2 <- function(Array.GmaxY, bandwidth= 45){
  nbtimesteps <- dim(Array.GmaxY)[1]
  nbsites.others <- dim(Array.GmaxY)[2]
  nbsites.interest <- dim(Array.GmaxY)[3]
  
  tt <- c(1:nbtimesteps)
  # timesteps <- seq(1,nbtimesteps,by=1)
  
  if (length(bandwidth) == 1){bandwidth <- rep(bandwidth,nbsites.interest)}
  
  if(length(bandwidth)!= nbsites.interest){
    stop("Bandwight must be a value or a vector of length nbsites.interest, i.e dim(Array.GmaxY)[3] ")}
  
  array.EGmaxHatY.mat <- array(NA,c(nbtimesteps,nbsites.others,nbsites.interest))
  
  for(i.interest in 1:nbsites.interest){
    #Gmaxhat_i.mat <- as.matrix(Array.GmaxY[,,i.interest])
    
    EGmaxhat_i.mat <- matrix(NA, ncol=nbsites.others,nrow=nbtimesteps)
    
    for (i.others in 1:nbsites.others){
      Kij <- outer(tt,tt,function(zz,z) dEpan((zz - z) / bandwidth[i.interest]))
      W <- Kij / rowSums(Kij)
      EGmaxhat_i.mat[,i.others] <- W %*% as.numeric(Array.GmaxY[,i.others,i.interest])
    }
    
    array.EGmaxHatY.mat[,,i.interest] <- EGmaxhat_i.mat
    
  }
  return(list("Array.EGmaxHatY" = array.EGmaxHatY.mat))
}

Kernel.marginal.estimation.step1step2 <- function(Data, site.interest = 0, site.others = 0, bandwidth.step1 = 35,bandwidth.step2 = 45){
  Array.step1 <- GmaxY.estimation.function.step1(Data = Data, site.interest = site.interest, site.others = site.others, bandwidth= bandwidth.step1)$Array.GmaxHatY
  Array.E.kernel.estimation <- EGmaxY.estimation.function.step2(Array.GmaxY= Array.step1, bandwidth= bandwidth.step2)$Array.EGmaxHatY
  Array.lambda.kernel.estimation <- (1-Array.E.kernel.estimation)/Array.E.kernel.estimation
  return(list("Marginal.probability.estimation" = Array.E.kernel.estimation,"Marginal.lambda.estimation" = Array.lambda.kernel.estimation))
} 

Marginal.aggregation.function.step3 <- function(Marginal.Array, dependence.type="log", dep.param = 0.3){
  
  marginal.lambdaHat <- (1-Marginal.Array)/Marginal.Array # Array
  
  nbtimesteps <- dim(marginal.lambdaHat)[1]
  nbsites.others <- dim(marginal.lambdaHat)[2]
  nbsites.interest <- dim(marginal.lambdaHat)[3]
  
  mat.aggregated.probability <- matrix(NA,ncol=nbsites.interest,nrow=nbtimesteps)
  mat.aggregated.lambda <- matrix(NA,ncol=nbsites.interest,nrow=nbtimesteps)
  
  if (dependence.type=="log"){
    for (i in 1:nbsites.interest){
      marginal.lambda.site.i <- as.matrix(marginal.lambdaHat[,,i])
      marginal.lambda.site.i.invSansV <- marginal.lambda.site.i^(-1)
      aggregated.lambda.site.i <- V(marginal.lambda.site.i.invSansV,dep.param)
      mat.aggregated.lambda[,i] <- aggregated.lambda.site.i
      aggregated.probability.site.i <- 1/(1+aggregated.lambda.site.i)
      mat.aggregated.probability[,i] <- aggregated.probability.site.i
    }
  }
  
  return(list("Kernel.aggregated.probability" = mat.aggregated.probability,"Kernel.aggregated.lambda" = mat.aggregated.lambda))
  
}

Estimation.RecordProbability <- function(Data, site.interest = 0, site.others = 0, 
                                         bandwidth.step1 = 35,bandwidth.step2 = 45,
                                         dependence.type="log", dep.param = 0.3){
  
  Array.step1 <- GmaxY.estimation.function.step1(Data = Data, site.interest = site.interest, site.others = site.others, bandwidth= bandwidth.step1)$Array.GmaxHatY
  Array.E.kernel.estimation <- EGmaxY.estimation.function.step2(Array.GmaxY= Array.step1, bandwidth= bandwidth.step2)$Array.EGmaxHatY
  
  Array.marginal.lambdaHat <- (1-Array.E.kernel.estimation)/Array.E.kernel.estimation 
  
  Aggregated.Kernel.estimation <- Marginal.aggregation.function.step3(Marginal.Array = Array.E.kernel.estimation, dependence.typ = dependence.type, dep.param = dep.param)
  
  Aggregated.probability <- Aggregated.Kernel.estimation$Kernel.aggregated.probability
  Aggregated.lambda <- Aggregated.Kernel.estimation$Kernel.aggregated.lambda
  
  return(list("Aggregated.probability.estimation" = Aggregated.probability,
              "Aggregated.lambda.estimation" = Aggregated.lambda,
              "Marginal.probability.estimation" = Array.E.kernel.estimation,
              "Marginal.lambda.estimation" = Array.marginal.lambdaHat))
} 

Compare.estimation.and.forecast.plot.func <- function(List.Forecast.object, List.Kernel.Estimation, Listes.True.quantities, type.plot="total", siteid = 0 , site.others = 1 , par1 = 3,par2 = 3 ){
  
  
  indexes.sites <- Listes.True.quantities$index_SitesInteret
  
  if(type.plot=="total"){
    
    par(mfrow=c(par1,par2), mar = c(5, 5.3, 2, 2))
    
    # lambda <- Listes.True.quantities$lambdaTotal.2toTnext_SitesInteret
    probability <- Listes.True.quantities$ProbabiliteTotal.2toTnext_SitesInteret
    nbtimesteps <- dim(probability)[1]
    HW.probability.evolution <- List.Forecast.object$HW.fit
    HW.probability.Tnext <- List.Forecast.object$HW.forecast
    probability.kernel <- List.Kernel.Estimation$Aggregated.probability.estimation
    # lambda.kernel <- List.Kernel.Estimation$Aggregated.lambda.estimation
    if(is.matrix(HW.probability.evolution) == FALSE){
      stop("Please use the aggregated forecast object as imput")
    }
    
    if(siteid[1]==0){
      siteid <- as.numeric(indexes.sites)
    }
    
    ylim.max <- max(probability,HW.probability.evolution)
    
    
    for (i in 1:length(siteid)){
      site.interest <- siteid[i]
      plot(c(1:nbtimesteps), probability[,site.interest], col='black', lwd=2, ylab=TeX(r'($P_{t}^{s_o}$)'), xlab="t from 2 to T+1",
           main=bquote("Record prob at site "~ .(site.interest)),ylim=c(0,ylim.max),xlim=c(1,nbtimesteps),type="l")
      points(c(4:nbtimesteps),HW.probability.evolution[,site.interest],col="green",pch=19)
      points(c(2:nbtimesteps),probability.kernel[,site.interest],pch=4,col="red")
      points(nbtimesteps,HW.probability.Tnext[site.interest],col="blue")
    }
    
  }
  
  
  if(type.plot=="marginal"){
    
    par(mfrow=c(par1,par2), mar = c(5, 5.3, 4, 2))
    
    # lambda <- Listes.True.quantities$lambdaMarginal.2toTnext_SitesInteret
    probability <- Listes.True.quantities$ProbabiliteMarginal.2toTnext_SitesInteret
    nbtimesteps <- dim(probability)[1]
    HW.probability.evolution <- List.Forecast.object$HW.fit
    HW.probability.Tnext <- List.Forecast.object$HW.forecast
    probability.kernel <- List.Kernel.Estimation$Marginal.probability.estimation
    # lambda.kernel <- List.Kernel.Estimation$Marginal.lambda.estimation
    
    if(is.array(HW.probability.evolution) == FALSE){
      stop("Please use the marginal forecast object as imput")
    }
    
    if(siteid[1]==0){
      siteid <- as.numeric(indexes.sites)
    }
    
    if(site.others[1]==0){
      site.others <- as.numeric(indexes.sites)
    }
    
    ylim.max <- max(probability,HW.probability.evolution)
    
    for(i in 1:length(siteid)){
      for (j in 1:length(site.others)){
        site.interest <- siteid[i]
        jplot <- site.others[j]
        plot(c(1:nbtimesteps), probability[,jplot,site.interest], col='blue', lwd=2, ylab=TeX(r'($P_{t}^{s,s_o}$)'), xlab="t from 2 to T+1",
             main=bquote(atop("Marginal Record prob at site " ~ .(site.interest), " w.r.t site " ~ .(jplot))),ylim=c(0,ylim.max),xlim=c(1,nbtimesteps),type="l")
        points(c(4:nbtimesteps),HW.probability.evolution[,jplot,site.interest],col="green",pch=19)
        points(c(2:nbtimesteps),probability.kernel[,jplot,site.interest],pch=4,col="red")
        points(nbtimesteps,HW.probability.Tnext[jplot,site.interest],col="blue")
      }
    }
    
  }
  
}

Compare.estimation.plot.func <- function(List.Kernel.Estimation, Listes.True.quantities, type.plot="total",
                                         plot.lambda = FALSE, siteid = 0, site.other = 0, par1 = 4, par2 = 2 ){
  
  if(siteid[1]==0){
    siteid <- as.numeric(Listes.True.quantities$index_SitesInteret)
  }
  
  if(type.plot=="total"){
    par(mfrow=c(par1,par2), mar = c(5, 5.3, 2, 2))
    
    Estimated.probability <- List.Kernel.Estimation$Aggregated.probability.estimation
    Estimated.lambda <- List.Kernel.Estimation$Aggregated.lambda.estimation
    True.probability <- Listes.True.quantities$ProbabiliteTotal.2toTnext_SitesInteret[-nrow(Listes.True.quantities$ProbabiliteTotal.2toTnext_SitesInteret), ]
    True.lambda <- Listes.True.quantities$lambdaTotal.2toTnext_SitesInteret[-nrow(Listes.True.quantities$ProbabiliteTotal.2toTnext_SitesInteret), ]
    nbtimesteps <- dim(True.lambda)[1]
    ylim.max.p <- max(True.probability,Estimated.probability)
    
    if(length(siteid)==1){
      plot(c(1:nbtimesteps),True.probability[,siteid],type="l", col="black",lwd=2,
           main=bquote("Record prob at site "~ .(siteid)),ylim=c(0,ylim.max.p),
           ylab=TeX(r'($P_t^{s_o}$)'),xlab="t from 2 to T")
      points(c(1:nbtimesteps),Estimated.probability[,siteid],col='red')
      if(plot.lambda == TRUE){
        ylim.max.l <- max(True.lambda, Estimated.lambda)
        plot(c(1:nbtimesteps),True.lambda[,siteid],type="l", col="black",lwd=2,
             main=bquote("Record prob at site "~ .(siteid)),ylim=c(0,ylim.max.l),
             ylab=TeX(r'($\lambda_t^{s_o}$)'),xlab="t from 2 to T")
        points(c(1:nbtimesteps),Estimated.lambda[,siteid],col='red')
      }
    }
    
    if(length(siteid)>1){
      for (i in 1:length(siteid)){
        siteid.i <- siteid[i]
        plot(c(1:nbtimesteps),True.probability[,siteid.i],type="l", col="black",lwd=2,
             main=bquote("Record prob at site "~ .(siteid.i)),ylim=c(0,ylim.max.p),
             ylab=TeX(r'($P_t^{s_o}$)'),xlab="t from 2 to T")
        points(c(1:nbtimesteps),Estimated.probability[,siteid.i],col='red')
        if(plot.lambda == TRUE){
          ylim.max.l <- max(True.lambda, Estimated.lambda)
          plot(c(1:nbtimesteps),True.lambda[,siteid.i],type="l", col="black",lwd=2,
               main=bquote("Associtaed lambda at site "~ .(siteid.i)), ylim=c(0,ylim.max.l),
               ylab=TeX(r'($\lambda_t^{s_o}$)'),xlab="t from 2 to T")
          points(c(1:nbtimesteps),Estimated.lambda[,siteid.i],col='red')
        }
      }
    }
  }
  
  
  if(type.plot=="marginal"){
    par(mfrow=c(par1,par2), mar = c(5, 5.3, 4, 2))
    if(site.other[1]==0){
      site.other <- as.numeric(Listes.True.quantities$index_SitesInteret)
    }
    Estimated.probability <- List.Kernel.Estimation$Marginal.probability.estimation
    Estimated.lambda <- List.Kernel.Estimation$Marginal.lambda.estimation
    True.probability <- Listes.True.quantities$ProbabiliteMarginal.2toTnext_SitesInteret[-(dim(Listes.True.quantities$ProbabiliteMarginal.2toTnext_SitesInteret)[1]),,]
    True.lambda <- Listes.True.quantities$lambdaMarginal.2toTnext_SitesInteret[-(dim(Listes.True.quantities$lambdaMarginal.2toTnext_SitesInteret)[1]),,]
    nbtimesteps <- dim(True.lambda)[1]
    ylim.max.p <- max(True.probability,Estimated.probability)
    # a[-(dim(a)[1]),,]
    if(length(siteid)==1){
      if(length(site.other)==1){
        plot(c(1:nbtimesteps),True.probability[,site.other,siteid],type="l", col="blue",lwd=2,
             main=bquote(atop("Marginal record prob at site "~ .(siteid)," w.r.t site " ~ .(site.other))),
             ylim=c(0,ylim.max.p), ylab=TeX(r'($P_t^{s,s_o}$)'),xlab="t from 2 to T")
        points(c(1:nbtimesteps),Estimated.probability[,site.other,siteid],col='red')
        if(plot.lambda == TRUE){
          ylim.max.l <- max(True.lambda, Estimated.lambda)
          plot(c(1:nbtimesteps),True.lambda[,site.other,siteid],type="l", col="blue",lwd=2,ylim=c(0,ylim.max.l),
               main=bquote(atop("Associated marginal "~ lambda ~ "at site "~ .(siteid)," w.r.t site " ~ .(site.other))),
               ylab=TeX(r'($\lambda_t^{s,s_o}$)'),xlab="t from 2 to T")
          points(c(1:nbtimesteps),Estimated.lambda[,site.other,siteid],col='red')
        }
      }
      if(length(site.other)>1){
        for (j in 1:length(site.other)){
          other.j <- site.other[j]
          plot(c(1:nbtimesteps),True.probability[,other.j,siteid],type="l", col="blue",lwd=2,
               main=bquote(atop("Marginal record prob at site "~ .(siteid)," w.r.t site " ~ .(other.j))),
               ylim=c(0,ylim.max.p), ylab=TeX(r'($P_t^{s,s_o}$)'),xlab="t from 2 to T")
          points(c(1:nbtimesteps),Estimated.probability[,other.j,siteid],col='red')
          if(plot.lambda == TRUE){
            ylim.max.l <- max(True.lambda, Estimated.lambda)
            plot(c(1:nbtimesteps),True.lambda[,other.j,siteid],type="l", col="blue",lwd=2,ylim=c(0,ylim.max.l),
                 main=bquote(atop("Associated marginal "~ lambda ~ "at site "~ .(siteid)," w.r.t site " ~ .(other.j))),
                 ylab=TeX(r'($\lambda_t^{s,s_o}$)'),xlab="t from 2 to T")
            points(c(1:nbtimesteps),Estimated.lambda[,other.j,siteid],col='red')
          }
        }
      }
    }
    
    if(length(siteid)>1){stop("please chose a site of interest")}
  }
}


Forecast.Holt.Winters.func <- function(input.traj){
  
  if(is.matrix(input.traj) == TRUE){ # on travaille avec les probas/lambdas aggregés.
    nbtimesteps <- dim(input.traj)[1]
    nbsites.interest <- dim(input.traj)[2]
    HW_2toT.fit.all <- matrix(NA,nrow=(nbtimesteps-2),ncol=nbsites.interest)
    HW_Tnext.fit.all <- rep(NA,nbsites.interest)
    for (i in 1:nbsites.interest){
      HW_2toT <- HoltWinters(ts(input.traj[,i], frequency = 1), gamma = FALSE)
      HW_2toT.fit  <- HW_2toT$fitted[,1]# green lines 
      HW_Tnext.fit <- predict(HW_2toT, n.ahead = 1, prediction.interval = TRUE)
      HW_2toT.fit.all[,i] <- HW_2toT.fit
      HW_Tnext.fit.all[i] <- HW_Tnext.fit[,1]
    }
  }
  
  if(is.matrix(input.traj) == FALSE){ # on travaille avec les probas/lambdas marginals
    
    nbtimesteps <- dim(input.traj)[1]
    nbsites.others <- dim(input.traj)[2]
    nbsites.interest <- dim(input.traj)[3]
    HW_2toT.fit.all <- array(NA,c((nbtimesteps-2),nbsites.others,nbsites.interest))
    HW_Tnext.fit.all <- matrix(NA,nrow=nbsites.others,ncol=nbsites.interest)
    
    for (i in 1:nbsites.interest){
      for (j in 1:nbsites.others){
        HW_2toT <- HoltWinters(ts(input.traj[,j,i], frequency = 1), gamma = FALSE)
        HW_2toT.fit  <- HW_2toT$fitted[,1]
        HW_Tnext.fit <- predict(HW_2toT, n.ahead = 1, prediction.interval = TRUE)
        HW_2toT.fit.all[,j,i] <- HW_2toT.fit
        HW_Tnext.fit.all[j,i] <- HW_Tnext.fit[,1]
      }
    }
  }
  return(list("HW.fit" = HW_2toT.fit.all,"HW.forecast" = HW_Tnext.fit.all))
}


##################################################################################
######################## All other functions #####################################
##################################################################################

K <- function(x) { # kernel 
  3/4 * (1 - x^2) * (abs(x) <= 1)
}

Gmaxtrue_untilt<- function(theta, sigma_t, ksi){ # true Gmax (prod of true G) "glissant" function
  function(x,t) {
    prod(pgev(x, theta, sigma_t[1:t], ksi))
    #pgev(x, theta, sigma_t, ksi)
  }
}

V <- function(sigmaMat,alph){ # logisitic function 
  sumv <- 0
  for (i in 1:(dim(sigmaMat)[2])){
    v <- 1/(sigmaMat[,i]^alph)
    sumv <- sumv + v }
  Vtotal <- sumv^(1/alph)
  return(Vtotal)
}

dEpan <- function(x){
  ## Function of Epanechnikov density distribution
  k <- (3/4)* (1-x^2)
  k[-1>x] <- 0
  k[x>1] <- 0
  return (k)
}

pgev <- function(x, theta, sigma, ksi) {# GEV CDF reparametrized with (xi,theta,sigma)
  exp( -pmax(0, (ksi/sigma * (x - theta)))^(-1/ksi) )
}

G_estimator <- function(t, xt, h) { # kernel nonparametric estimator of the CDF G
  function(x, t0) {
    Kvect <- K((t - t0) / h)
    sum(Kvect * (xt <= x)) / sum(Kvect)
  }
}

Gmax_glissant_estimator <- function(t, Ghat) { # nonparametric estimator of the CDF Gmax until time t-1, where t is our time of interest
  function(x) {  
    # sum()
    tmoins1 <- c(1:(t-1))
    prod(sapply(tmoins1, Ghat, x = x))
    # sapply(t, Ghat, x = x)
  }
}

RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

##################################################################################

# Funtions for TRUE prbability
GenerateStationarySmaple <- function(Ttoday, alpha, S, xi ,mu0, sigma0){
  dat <- evd::rmvevd(Ttoday, dep= alpha, model = c("log"), d = S, mar=c(mu0, sigma0,xi))
  return(dat)
}

AddTrendToData <- function(dat.mat, sigma.mat, support){
  datgev.mat <- matrix(NA, nrow=dim(dat.mat)[1], ncol=dim(dat.mat)[2])
  mu.mat <- support + sigma.mat/xi
  for ( i in 1:dim(dat.mat)[2]){
    datgev.mat[,i] <- dat.mat[,i]*sigma.mat[,i] + mu.mat[,i]
  }
  return(datgev.mat)
}

TruePnext <- function(sigma.mat, sigma.next, xi, alpha){
  sigma_s.mat <- t(t(apply(sigma.mat^(1/xi), 2,sum)))^xi
  lambda_next_s.mat <- (sigma_s.mat/sigma.next)^(1/xi)
  lambda_next_s_inv_sansV <- lambda_next_s.mat^(-1)
  lambda_next_s_inv_sansV <- t(as.numeric(lambda_next_s_inv_sansV))
  lambda_next <- V(lambda_next_s_inv_sansV,alpha)
  TruePrecord <- 1/(1+lambda_next)
  return(TruePrecord)
}

TrueP1toT <- function(sigma.mat, sigma.next, xi, alpha){
  TruePrecord <- matrix(NA, ncol=(dim(sigma.mat)[2]+1), nrow=dim(sigma.mat)[1])
  TrueLambdas <- matrix(NA, ncol=(dim(sigma.mat)[2]+1), nrow=dim(sigma.mat)[1])
  for (i in 1: dim(sigma.mat)[2]){
    sigma_maxst <- cumsum((sigma.mat[,i]^(1/xi)))^xi
    lambda_st <- (sigma_maxst/ c(sigma.mat[-1,1],sigma_next))^(1/xi)
    P_record <- 1/(1+lambda_st)
    TruePrecord[,i] <- P_record
    TrueLambdas[,i] <- lambda_st
  }
  lambdaMatMatinv <- TrueLambdas^(-1)
  lambda_Vt <- V(lambdaMatMatinv[,1:dim(sigma.mat)[2]],alpha)
  Truerecord_all <- 1/(1+lambda_Vt)
  TrueLambdas[,(dim(sigma.mat)[2]+1)] <- lambda_Vt
  TruePrecord[,(dim(sigma.mat)[2]+1)] <- Truerecord_all
  return(list(TruePrecord,TrueLambdas))
}

# Functions for estimated probability
HatGmaxs <- function(dat.mat, theta, sigma.mat, xi , h1, t){
  theta <- theta[1]
  Gmaxhat.mat <- matrix(NA, ncol=dim(sigma.mat)[2], nrow=dim(sigma.mat)[1])
  Gmax.mat <- matrix(NA, ncol=dim(sigma.mat)[2], nrow=dim(sigma.mat)[1])
  for (i in 1: dim(sigma.mat)[2]){
    TrueGmax <- Gmaxtrue_untilt(theta=theta, sigma_t=sigma.mat[,i], ksi=xi)
    Ghats <- G_estimator(t, dat.mat[,i], h1[i])
    for (j in 2: dim(sigma.mat)[1]){
      Gmax.mat[j,i] <- TrueGmax(dat.mat[j,1],(j-1))
      
      Gmaxhats_tj <- Gmax_glissant_estimator(t=j,Ghats)
      Gmaxhat_s <- Gmaxhats_tj(dat.mat[j,1]) 
      Gmaxhat.mat[j,i] <- Gmaxhat_s
    }
    
  }
  Gmaxhat.mat <-as.matrix(Gmaxhat.mat[-1,])
  Gmax.mat <- as.matrix(Gmax.mat[-1,])
  return(list(Gmaxhat.mat,Gmax.mat))
}

HatEGmaxs <- function(h2, t, Gmaxhat.mat, Gmax.mat){
  EGmaxhat.mat <- matrix(NA, ncol=dim(Gmaxhat.mat)[2],nrow=dim(Gmaxhat.mat)[1])
  EGmax.mat <- matrix(NA, ncol=dim(Gmax.mat)[2],nrow=dim(Gmax.mat)[1])
  for (i in 1:dim(Gmaxhat.mat)[2]){
    Kij <- outer(t[-1],t[-1],function(zz,z) dEpan((zz - z) / h2[i])) 
    W <- Kij / rowSums(Kij)
    EGmax.mat[,i] <- W %*% as.numeric(Gmax.mat[,i])
    EGmaxhat.mat[,i] <- W %*% as.numeric(Gmaxhat.mat[,i])
  }
  return(list(EGmaxhat.mat,EGmax.mat))
}

HatPnext <- function(EGmaxhat.mat,EGmax.mat){
  HWpred.mat <- rep(0, dim(EGmax.mat)[2])
  HWhatpred.mat <-rep(0, dim(EGmaxhat.mat)[2])
  Lambdapred.mat <- rep(0, dim(EGmax.mat)[2])
  HatLambdapred.mat <-rep(0, dim(EGmaxhat.mat)[2])
  for (i in 1:dim(EGmax.mat)[2]){
    HWGmax_s <- HoltWinters(ts(EGmax.mat[,i], frequency = 1), gamma = FALSE)
    HWGmaxpred_s <- predict(HWGmax_s , n.ahead = 1, prediction.interval = TRUE)
    HWpred.mat[i] <- HWGmaxpred_s[,1]
    
    HWGmaxhat_s <- HoltWinters(ts(EGmaxhat.mat[,i], frequency = 1), gamma = FALSE)
    HWGmaxhatpred_s <- predict(HWGmaxhat_s, n.ahead = 1, prediction.interval = TRUE)
    HWhatpred.mat[i] <- HWGmaxhatpred_s[,1]
    
    
    Lambdapred.mat[i] <- 1/(HWGmaxpred_s[,1]/(1-HWGmaxpred_s[,1]))
    HatLambdapred.mat[i] <- 1/(HWGmaxhatpred_s[,1]/(1-HWGmaxhatpred_s[,1]))
  }
  return(list(HWhatpred.mat,HWpred.mat,HatLambdapred.mat,Lambdapred.mat ))
}

LambdaSfunc <- function(Lambdapred.mat,alpha){
  Lambdapred.mat <- as.matrix(Lambdapred.mat, ncol=length(Lambdapred.mat) )
  lambdaMatMatinv_from_EGmax_hat <- Lambdapred.mat^(-1)
  lambda_next_from_EGmax_hat <- V(t(lambdaMatMatinv_from_EGmax_hat),alpha)
  Pnext_hat <- 1/(1+lambda_next_from_EGmax_hat)
  return(Pnext_hat)
}


############################################################################
############################################################################
# En cours 
############################################################################
############################################################################


# Ratio w.r.t theoretical stationary 
############################################################################

Compare_stationnary_fonction <- function(Estimated.probability2toTnext.mat, sigma.mat, xi, site.interest = 0, dependence.type="log", dep.param){
  nbtimesteps <- dim(Estimated.probability2toTnext.mat)[1]
  timesteps <- seq(1,nbtimesteps,by=1)
  nbsites <- dim(Estimated.probability2toTnext.mat)[2]
  
  sigma_stationary.vec <- as.numeric(sigma.mat[1,])
  
  if(length(site.interest) > nbsites){
    stop("the number of sites of study exceeds the number of columns of data ")}
  
  if(site.interest[1] == 0){site.interest <- seq(1,dim(Estimated.probability)[2],by=1)} # all the sites where we did the estimation are the interest sites
  
  p.stationary.mat <- matrix(NA, ncol=length(site.interest), nrow = nbtimesteps)
  ratio_site.mat <- matrix(NA, ncol=length(site.interest), nrow = nbtimesteps)
  
  t_i_minus1 <- (seq(from=1, to=(nbtimesteps), by=1))
  
  for (i in site.interest){
    lambda_stationary.mat <- matrix(NA, ncol=nbsites, nrow = nbtimesteps)
    sigma_intermediate <- (sigma_stationary.vec/sigma_stationary.vec[i])^(1/xi) 
    lambda.stationary.marginal.siteid <- outer(t_i_minus1,sigma_intermediate,FUN = "*")
    
    if (dependence.type == "log"){
      lambda.stationary.marginal.siteid.invSansV <- lambda.stationary.marginal.siteid^(-1)
      lam.aggregated.stationary.siteid <- V(lambda.stationary.marginal.siteid.invSansV,dep.param)
      p.stationary.siteid <- 1/(1+lam.aggregated.stationary.siteid)
    }
    
    p.stationary.mat[,i] <- as.numeric(p.stationary.siteid)
    ratio_site.i <- as.numeric(Estimated.probability2toTnext.mat[,i]/p.stationary.siteid)
    ratio_site.mat[,i] <- ratio_site.i # each colum is a site of interest 
  }
  
  return(list("ratio" = ratio_site.mat,
              "p.stationary" = p.stationary.mat))
}








