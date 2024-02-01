# Test

############################################################################
############################################################################
# Temperature
############################################################################
############################################################################

# Parameters
############################################################################
############################################################################
xi <- -0.2
support <- 50
alpha <- 0.3
S <- 8
lengthfrom1toTnext <- 86
h1 <- 35
h2 <- 45

# Generare Data
############################################################################
############################################################################
Simulated_obs_and_Tnext <- DataStructure.func(Tnext=lengthfrom1toTnext,S=S,xi=xi,dependence.param=alpha)
SimulatedData <- Simulated_obs_and_Tnext$SimulatedObservations
obsTnext <- Simulated_obs_and_Tnext$SimulatedNext

# Plot : Data
############################################################################

# Plot : the trajectories independly
PlotData.separate(SimulatedData) 
PlotData.separate(SimulatedData, add.max=TRUE) 
PlotData.separate(SimulatedData,siteid = 1, add.max=TRUE, par1 = 1, par2 = 1) 

# Plot : compare trajectories
PlotData.compare(SimulatedData, site.interest = 1,  max.plot = TRUE)
PlotData.compare(SimulatedData, site.interest = S,  max.plot = TRUE)

# Estimation - step by step
############################################################################
############################################################################

# Step 1 : Estimators GmaxHat(Y)
GmaxHatY <- GmaxY.estimation.function.step1(SimulatedData, bandwidth = h1)$Array.GmaxHatY
dim(GmaxHatY) # dim1 : times steps associated to estimations, t = 2 .... Ttoday
              # dim2 : GHatmax,s(Y_so) for s in 1,...,S and fixed so
              # dim3 : GHatmax,s(Y_so) for fixed s and so in 1,...,S 

# Step 2 : Estimator EGmaxHat(Y)
E.GmaxY.test <- EGmaxY.estimation.function.step2(GmaxHatY, bandwidth = h2)$Array.EGmaxHatY
dim(E.GmaxY.test) # Same dimentions as GmaxHatY 

# Step 3 : Estimator EGmaxHat(Y)
Estimated.aggregated.probability.and.lambda <- Marginal.aggregation.function.step3( Marginal.Array = E.GmaxY.test, dep.param = alpha) # we consider a logistic known dependance structure
Aggregated.lambdaHat<- Estimated.aggregated.probability.and.lambda$Kernel.aggregated.probability
Aggregated.lambdaHat <- Estimated.aggregated.probability.and.lambda$Kernel.aggregated.lambda

# Estimation : All steps in 1 function
############################################################################

Kernel.estimaton <- Estimation.RecordProbability(SimulatedData, bandwidth.step1 = h1,bandwidth.step2 = h2, dep.param = alpha)
Estimated.probability <- Kernel.estimaton$Aggregated.probability.estimation
Estimated.lambda <- Kernel.estimaton$Aggregated.lambda.estimation
Estimated.marginal.probability <- Kernel.estimaton$Marginal.probability.estimation
Estimated.marginal.lambda <- Kernel.estimaton$Marginal.lambda.estimation

# Theoretical results
############################################################################
############################################################################

# Theoretical parameters : ( Here scale is the one used by defalt
#                            in the generatin function DataStructure.func() )

mu.mat <- as.matrix(replicate(S,seq(20, 35, length.out = lengthfrom1toTnext )),ncol=S)
sigma.mat <- -xi*(support-mu.mat) 

True.quantities <- TrueP2toTnext(sigma.mat = sigma.mat, xi = xi, alpha=alpha) # list of 5 elements 

True.quantities$index_SitesInteret

True.marginal.probability <- True.quantities$ProbabiliteMarginal.2toTnext_SitesInteret
dim(True.marginal.probability) # dim1 : times steps , t = 2 .... Tnext
                               # dim2 : marginal P_{s,so} for s in 1,...,S and fixed so
                               # dim3 : marginal P_{s,so} for fixed s and so in 1,...,S 

True.marginal.lambda <- True.quantities$lambdaMarginal.2toTnext_SitesInteret
dim(True.marginal.lambda) # dimentions : same as before

True.probability <- True.quantities$ProbabiliteTotal.2toTnext_SitesInteret
dim(True.probability) # dim1 : nb times steps , t = 2 .... Tnext
                      # dim2 == 1 
                      # dim3 : sites of study , here $index_SitesInteret

True.lambda <- True.quantities$lambdaTotal.2toTnext_SitesInteret
dim(True.lambda) # same as before

# Plot : Theoretical (input is the list obtained directly from function TrueP2toTnext() )
############################################################################

Plot.ProbLambda(True.quantities)
Plot.ProbLambda(True.quantities, siteid=0, plot.type = "total", par1=4,par2=4)
Plot.ProbLambda(True.quantities, siteid=1, site.other = 1,plot.type = "marginal")
Plot.ProbLambda(True.quantities, siteid=1, site.other = 0,plot.type = "marginal", par1=4,par2=4)

# Plot : Compare estimation and theoretical values
############################################################################

Compare.estimation.plot.func(Kernel.estimaton, True.quantities,siteid = c(1,2), par1 = 1, par2 = 2)
Compare.estimation.plot.func(Kernel.estimaton, True.quantities,siteid = 1, plot.lambda = TRUE ,par1 = 1, par2 = 2)
Compare.estimation.plot.func(Kernel.estimaton, True.quantities)
Compare.estimation.plot.func(Kernel.estimaton, True.quantities,plot.lambda = TRUE )
Compare.estimation.plot.func(Kernel.estimaton, True.quantities,type.plot="marginal",siteid = 1,par1 = 3, par2 = 3)
# IL FAUT QUE J'AMELIORE CETTE FONCTION POUR QUE L"ECHELLE NE GENE PAS , 
# PEUT ETRE FAIRE UN ZOOM ???

# Forecast : Holt-Winters 
############################################################################
############################################################################

# Forecast of aggregated probabilities
Forecasted.probabilities.4toTnext <- Forecast.Holt.Winters.func(Estimated.probability)
HW.probabilities.observed.fit <- Forecasted.probabilities.4toTnext$HW.fit
HW.probabilities.Tnext.fit <- Forecasted.probabilities.4toTnext$HW.forecast

# Forecast of marginal probability 
Forecasted.marginal.probabilities.4toTnext <- Forecast.Holt.Winters.func(Estimated.marginal.probability)
HW.marginal.probabilities.observed.fit <- Forecasted.marginal.probabilities.4toTnext$HW.fit
HW.marginal.probabilities.Tnext.fit <- Forecasted.marginal.probabilities.4toTnext$HW.forecast

# Plot : forecast evolution
############################################################################
Compare.estimation.and.forecast.plot.func(List.Forecast.object = Forecasted.probabilities.4toTnext ,
                                          List.Kernel.Estimation = Kernel.estimaton,
                                          Listes.True.quantities = True.quantities)

Compare.estimation.and.forecast.plot.func(List.Forecast.object = Forecasted.probabilities.4toTnext ,
                                          List.Kernel.Estimation = Kernel.estimaton,
                                          Listes.True.quantities = True.quantities, siteid=c(1,2), par1 = 1,par2 = 2)

Compare.estimation.and.forecast.plot.func(List.Forecast.object = Forecasted.marginal.probabilities.4toTnext,
                                          List.Kernel.Estimation = Kernel.estimaton,
                                          Listes.True.quantities = True.quantities,
                                          type.plot="marginal", siteid=c(1,2), par1 = 1,par2 = 2)


Compare.estimation.and.forecast.plot.func(List.Forecast.object = Forecasted.marginal.probabilities.4toTnext,
                                          List.Kernel.Estimation = Kernel.estimaton,
                                          Listes.True.quantities = True.quantities,
                                          type.plot="marginal", siteid=c(1,2), site.others = c(1,2) ,par1 = 2,par2 = 2)

# Plot : Forecast vs True at Tnext 
############################################################################

# Aggregated
cat("Forecasted vs True probabilities : ","\n",
    True.quantities$index_SitesInteret,"\n", "True : ","\n",
    True.quantities$ProbabiliteTotal.2toTnext_SitesInteret[lengthfrom1toTnext-1,],"\n",
    "Forecast : ","\n",HW.probabilities.Tnext.fit)

# Marginal
indx.interest <- c(1,2,3) # True.quantities$index_SitesInteret
indx.others <- c(1,2)
for (i in indx.interest){
  cat("Forecasted vs True probabilities : ","\n",
      "indx interest : ", i,"\n", 
      "indx others : ", indx.others,"\n", 
      "True (col = site.interest, row :site.other): ","\n",
      as.matrix(True.quantities$ProbabiliteMarginal.2toTnext_SitesInteret[lengthfrom1toTnext-1,indx.others,i]),"\n",
      "Forecast : ","\n",as.matrix(HW.marginal.probabilities.Tnext.fit[indx.others,i]),"\n")
}

# Ratio w.r.t stationary 
############################################################################
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

############################################################################
############################################################################

Estimated.probability <- Kernel.estimaton$Aggregated.probability.estimation
Estimated.P2toTnext <- rbind(Estimated.probability,HW.probabilities.Tnext.fit)

EvolRation.estimated <- Compare_stationnary_fonction(Estimated.probability2toTnext.mat = Estimated.P2toTnext, xi = xi,
                                                     sigma.mat = sigma.mat, site.interest = 0, dependence.type="log", dep.param = alpha)
# P.theo.stationary <- EvolRation.estimated$p.stationary
Estimated.ratio <- EvolRation.estimated$ratio
True.probability # True.quantities <- TrueP2toTnext(sigma.mat = sigma.mat, xi = xi, alpha=alpha) # list of 5 elements 
                 # # True.probability <- True.quantities$ProbabiliteTotal.2toTnext_SitesInteret
Theoretical.ratio <- True.probability/EvolRation.estimated$p.stationary

par(mfrow=c(2,4), mar = c(5, 5.3, 2, 2))
ylim.min <- min(Theoretical.ratio,Estimated.ratio)
ylim.max <- max(Theoretical.ratio,Estimated.ratio)
for( i in 1:S){
  plot(c(1:85),Theoretical.ratio[,i],type="l",lwd=2,ylim=c(ylim.min,ylim.max), xlab="t from to T+1", ylab="Ratio Nonstationarity vs Stationarity",
       main=bquote("Ratio at site "~ .(i)))
  abline(h=1,col="gray",lwd = 2); abline(h=2,col="gray",lwd = 1,lty=2); abline(h=3,col="gray",lwd = 1,lty=2); abline(h=4,col="gray",lwd = 1,lty=2)
  abline(h=5,col="gray",lwd = 1,lty=2); abline(h=6,col="gray",lwd = 1,lty=2); abline(h=7,col="gray",lwd = 1,lty=2); abline(h=8,col="gray",lwd = 1,lty=2)
  points(c(1:85),Estimated.ratio[,i], col="red")
  legend('topleft',legend=c("Theoretical ratio", 'Estimated ratio'),col=c("black","red"), lty=c(1,0), pch=c(NA,1),bg="transparent")
  
}




############################################################################
############################################################################

# Compare multiple simulations 
############################################################################
############################################################################

N <- 200
S <- 8 
Nmean.Estimated.probability.from2toT <- array(NA,c(N,lengthfrom1toTnext-2,S))
N_forecast.aggregated.probability <- matrix(NA,ncol=N,nrow=S)
for (n in 1:N){
  Simulated_obs_and_Tnext <- DataStructure.func(Tnext=lengthfrom1toTnext,S=S,xi=xi,dependence.param=alpha)
  SimulatedData <- Simulated_obs_and_Tnext$SimulatedObservations

  Kernel.estimaton.n <- Estimation.RecordProbability(SimulatedData, bandwidth.step1 = h1,bandwidth.step2 = h2, dep.param = alpha)
  Estimated.probability.n <- Kernel.estimaton.n$Aggregated.probability.estimation
  Nmean.Estimated.probability.from2toT[n,,] <- Estimated.probability.n 
  
  HW.forecast.n <- Forecast.Holt.Winters.func(Estimated.probability.n)$HW.forecast
  N_forecast.aggregated.probability[,n] <- HW.forecast.n
}

# Boxplots: N forecasts 
############################################################################
# forcastast for site 1
N_forecast.aggregated.probability[1,]
# forcastast for site S
N_forecast.aggregated.probability[S,]

# mean forecats 
Nmean.forecast <- rowMeans(N_forecast.aggregated.probability)
Nmean.forecast 

# True
True.compare.boxplot <- True.quantities$ProbabiliteTotal.2toTnext_SitesInteret[85,]

# boxplot forecast ALL 
boxplot(t(N_forecast.aggregated.probability), xlab="sites", ylab="Probability",
        main="Fotecasted spatio-temporal record probability",col=rgb(1.0,0,0,alpha=0.3))
abline(h=True.compare.boxplot,col="cyan3",lwd = 2)
text(x=1,y=True.compare.boxplot*1.1,round(True.compare.boxplot,3),col="cyan4",cex=1.2)

# boxplot forecast site 1 
site.plot <- 1
boxplot(N_forecast.aggregated.probability[site.plot,],
        main=bquote(" Forecasted spatio-temporal record probability at site "~ .(site.plot)),
        xlab=bquote("site "~.(site.plot)),ylab="Probability",col=rgb(1.0,0,0,alpha=0.3))
abline(h=True.compare.boxplot,col="cyan3",lwd = 2)
text(x=0.8,y=True.compare.boxplot*1.1,round(True.compare.boxplot,3),col="cyan4",cex=1.2)


# Compare to probability in a stationnary situation
############################################################################
############################################################################

