### Adaptive MCMC to fit "shifting-kappa" model simultaneously to resource supply and periodic 
### starvation experiments using Biomphalaria glabrata and Schistosoma mansoni

library("adaptMCMC") # needed for MCMC
library("deSolve") # needed to simulate the models

# The model definition is written in C, so we need to use Rtools and load the model as a dll
rtools <- "C:\\Rtools\\bin"
gcc <- "C:\\Rtools\\gcc-4.6.3\\bin"
path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
new_path <- c(rtools, gcc, path)
new_path <- new_path[!duplicated(tolower(new_path))]
Sys.setenv(PATH = paste(new_path, collapse = ";"))

# compile my model from C definition
setwd("C:/RData")
dyn.unload("IndividualModel_YVE.dll") # unload dll
system("R CMD SHLIB IndividualModel_YVE.c")
dyn.load("IndividualModel_YVE.dll") # Load dll


#### Fixed information ####

# 1 - Data
setwd("C:/RData")

data = read.csv("ResourceSupplyExp.csv")
data2 = read.csv("PeriodicStarveExp.csv")

data = list(t = data$Date, L = data$Length, Negg = data$C_Eggs, Nworms = data$C_Worms, Alive=data$Last_Alive,
            L2 = data2$Length, Negg2 = data2$C_Eggs, Nworms2 = data2$C_Worms, Alive2=data2$Last_Alive)


# 2 - Initial conditions
setinits.Food<-function(F0 = 16.5, L0=0.7, e0=0.9, D0 = 0, RH0 = 0, P0 = 0, RP0 = 0, DAM0=0, HAZ0=0){
  inits<-c(F = F0, L = L0, e = e0, D = D0, RH = RH0, P = P0, RP = RP0, DAM=DAM0, HAZ=HAZ0)
  return(inits)
}

setinits.Starve<-function(F0 = 10.74, L0=0.7, e0=0.9, D0 = 0, RH0 = 0, P0 = 0, RP0 = 0, DAM0=0, HAZ0=0){
  inits<-c(F = F0, L = L0, e = e0, D = D0, RH = RH0, P = P0, RP = RP0, DAM=DAM0, HAZ=HAZ0)
  return(inits)
}


# 3 - Functions for feeding events (these match the timing of feeding events in the actual ecxperiments)

# For the starvation experiment
periodic.starvation.events = function(initial.food, starve.period, infection.date, weeks.duration){
  infection.value = 2.85e-5
  all.fed.dates = infection.date + c(1,4) # Snails were infected on a Thursday, All fed on Friday and Monday, then treatments
  potential.feedings = infection.date + 5 + sort(c((1:weeks.duration)*7 - 6,((1:weeks.duration)*7 - 4),((1:weeks.duration)*7 - 2)))
  #print(0:(weeks.duration/starve.period - 1)*6 + 1)
  potential.feedings.vals = rep(initial.food, times = length(potential.feedings))
  # work out starvation
  if( starve.period == 2){
    potential.feedings.vals[(potential.feedings - (infection.date + 5))  %% 14 <= 7] = 0
  }
  if (starve.period == 3){
    potential.feedings.vals[(potential.feedings - (infection.date + 5)) %% 21 <= 14] = 0
  }
  if (starve.period == 4){
    potential.feedings.vals[(potential.feedings - (infection.date + 5)) %% 28 <= 21] = 0
  }
  event.dates = c(infection.date, all.fed.dates, potential.feedings, infection.date)
  event.values = c(infection.value, initial.food, initial.food, potential.feedings.vals, 0)
  methods = rep("replace", times=length(event.dates))
  vars = c("P", rep("F", times = length(event.dates)-2), "HAZ")
  data.frame(var=vars, time= event.dates, value = event.values, method= methods)
}

# For the resource supply experiment
feeding.events <- function(dates, var="F", value, method="replace", Infected=0){
  # Assemble all data
  events = length(dates) # number of events
  vars = rep(var, times = events)
  values = rep(value, times = events)
  methods = rep(method, times = events)
  
  #build data.frame
  result = data.frame(var = vars, time = dates, value = values, method = methods)
  if(Infected > 0){
    result = rbind(result, data.frame(var="P", time=Infected, value=2.85e-5, method="replace"))
  }
  result = rbind(result, data.frame(var="HAZ", time=63, value=0, method="replace"))
  result
}

# Duration of experiments
dur.R = 245 # Resource supply
dur.P = 140 # Periodic starvation

in.R = setinits.Food()
in.P = setinits.Starve()
Feed.R = feeding.events(dates = sort(c((1:35)*7,((1:35)*7 - 3))), var="F", in.R[1], method="replace", Infected=28)


# 4 - Functions to solve DEBs
solve.DEB.supply<-function(params, inits, duration, feeding.events){
  feed.sup <- feeding.events
  feed.sup.U <- subset(feed.sup, var != "P")
  parms = as.numeric(params[1:24])
  params = c(iM=parms[1], k=parms[2], M=parms[3], EM=parms[4], Fh=parms[5], muD=parms[6],
             DR=parms[7], fe=parms[8], yRP=parms[9], ph=parms[10], yPE=parms[11], iPM=parms[12],
             eh=parms[13], mP=parms[14], alpha=parms[15], yEF=parms[16], LM=parms[17],kR=parms[18], 
             delta0=parms[19], hdelta=parms[20], hb=parms[21], theta=parms[22], mR=parms[23], yVE=parms[24], startage=28)
  
  Sup.6 <- data.frame(lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                            initfunc = "initmod",  nout=1, outnames="Survival", maxsteps=5e5,
                            params[1:25],  rtol=1e-6, atol=1e-6,   
                            events = list(data = feed.sup)))
  
  Sup.6U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                  initfunc = "initmod",  nout=1, outnames="Survival", maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.sup.U))
  
  feed.sup[1:70,3] <- 11
  feed.sup.U[1:70,3] <- 11
  Sup.5 <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                 initfunc = "initmod",  nout=1, outnames="Survival", maxsteps=5e5,
                 params[1:25],  rtol=1e-6, atol=1e-6,   
                 events = list(data = feed.sup))
  if(attributes(Sup.5)$istate[1] != 2)(return(Sup.6))
  Sup.5U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                  initfunc = "initmod",  nout=1, outnames="Survival", maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.sup.U))
  if(attributes(Sup.5U)$istate[1] != 2)(return(Sup.6))
  
  feed.sup[1:70,3] <- 5.5
  feed.sup.U[1:70,3] <- 5.5
  Sup.4 <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                 initfunc = "initmod",  nout=1, outnames="Survival", maxsteps=5e5,
                 params[1:25],  rtol=1e-6, atol=1e-6,   
                 events = list(data = feed.sup))  
  if(attributes(Sup.4)$istate[1] != 2)(return(Sup.6))
  Sup.4U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                  initfunc = "initmod",  nout=1, outnames="Survival", maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.sup.U))
  if(attributes(Sup.4U)$istate[1] != 2)(return(Sup.6))
  
  feed.sup[1:70,3] <- 2.75
  feed.sup.U[1:70,3] <- 2.75
  Sup.3 <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                 initfunc = "initmod",  nout=1, outnames="Survival", maxsteps=5e5,
                 params[1:25],  rtol=1e-6, atol=1e-6,   
                 events = list(data = feed.sup))
  if(attributes(Sup.3)$istate[1] != 2)(return(Sup.6))
  Sup.3U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                  initfunc = "initmod",  nout=1, outnames="Survival", maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.sup.U))
  if(attributes(Sup.3U)$istate[1] != 2)(return(Sup.6))
  
  feed.sup[1:70,3] <- 1.375
  feed.sup.U[1:70,3] <- 1.375
  Sup.2 <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                  initfunc = "initmod",  nout=1, outnames="Survival", maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.sup))
  if(attributes(Sup.2)$istate[1] != 2)(return(Sup.6))
  Sup.2U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                   initfunc = "initmod",  nout=1, outnames="Survival", maxsteps=5e5,
                   params[1:25],  rtol=1e-6, atol=1e-6,   
                   events = list(data = feed.sup.U))
  if(attributes(Sup.2U)$istate[1] != 2)(return(Sup.6))
  
  feed.sup[1:70,3] <- 0.6875
  feed.sup.U[1:70,3] <- 0.6875
  Sup.1 <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                  initfunc = "initmod",  nout=1, outnames="Survival", maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.sup))
  if(attributes(Sup.1)$istate[1] != 2)(return(Sup.6))
  Sup.1U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                   initfunc = "initmod",  nout=1, outnames="Survival", maxsteps=5e5,
                   params[1:25],  rtol=1e-6, atol=1e-6,   
                   events = list(data = feed.sup.U))
  if(attributes(Sup.1U)$istate[1] != 2)(return(Sup.6))
  
  result <- rbind(#Infecteds (n=66)
    Sup.1, Sup.1, Sup.1, Sup.1, Sup.1, Sup.1, Sup.1, Sup.1, Sup.1, Sup.1, 
    Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, Sup.2, 
    Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, Sup.3, 
    Sup.4, Sup.4, Sup.4, Sup.4, Sup.4, Sup.4, Sup.4, Sup.4, Sup.4, Sup.4, 
    Sup.5, Sup.5, Sup.5, Sup.5, Sup.5, Sup.5, Sup.5, Sup.5, Sup.5, Sup.5, 
    Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, Sup.6, 
    #Uninfecteds (n=30)
    Sup.1U, Sup.1U, Sup.1U, Sup.1U, Sup.1U, 
    Sup.2U, Sup.2U, Sup.2U, Sup.2U, Sup.2U, 
    Sup.3U, Sup.3U, Sup.3U, Sup.3U, Sup.3U, 
    Sup.4U, Sup.4U, Sup.4U, Sup.4U, Sup.4U, 
    Sup.5U, Sup.5U, Sup.5U, Sup.5U, Sup.5U,
    Sup.6U, Sup.6U, Sup.6U, Sup.6U, Sup.6U)
  
  result
  
}

# This function is customized to my model and data
solve.DEB.starve<-function(params, inits, duration){
  # Collect the params the way C likes them
  parms = as.numeric(params[1:25])
  params = c(iM=parms[1], k=parms[2], M=parms[3], EM=parms[4], Fh=parms[5], muD=parms[6],
             DR=parms[7], fe=parms[8], yRP=parms[9], ph=parms[10], yPE=parms[11], iPM=parms[12],
             eh=parms[13], mP=parms[14], alpha=parms[15], yEF=parms[25], LM=parms[17],kR=parms[18], 
             delta0=parms[19], hdelta=parms[20], hb=parms[21], theta=parms[22], mR=parms[23], yVE=parms[24], startage=14)
  
  # Simulate dynamics for 1-0
  food_1_0 = periodic.starvation.events(2.69, 0, 14, 18)
  food_1_0.U <- subset(food_1_0, var != "P")
  out_1_0 <- data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                            initfunc = "initmod",  nout=1, outnames="Survival", method="lsoda",
                            params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=500000,
                            events = list(data = food_1_0)))
  out_1_0U <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                  initfunc = "initmod",  nout=1, outnames="Survival", method="lsoda",
                  params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=500000,
                  events = list(data = food_1_0.U))
  if(attributes(out_1_0U)$istate[1] != 2)(return(out_1_0))
  
  # Simulate dynamics for 2-2
  food_2_2 = periodic.starvation.events(5.37, 2, 14, 18)
  food_2_2.U <- subset(food_2_2, var != "P")
  out_2_2 <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                 initfunc = "initmod",  nout=1, outnames="Survival", method="lsoda",
                 params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=500000,
                 events = list(data = food_2_2))
  if(attributes(out_2_2)$istate[1] != 2)(return(out_1_0))
  out_2_2U <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                  initfunc = "initmod",  nout=1, outnames="Survival", method="lsoda",
                  params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=500000,
                  events = list(data = food_2_2.U))
  if(attributes(out_2_2U)$istate[1] != 2)(return(out_1_0))
  
  # Simulate dynamics for 3-3
  food_3_3 = periodic.starvation.events(8.06, 3, 14, 18)
  food_3_3.U <- subset(food_3_3, var != "P")
  out_3_3 <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                 initfunc = "initmod",  nout=1, outnames="Survival", method="lsoda",
                 params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=500000,
                 events = list(data = food_3_3))
  if(attributes(out_3_3)$istate[1] != 2)(return(out_1_0))
  out_3_3U <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                  initfunc = "initmod",  nout=1, outnames="Survival", method="lsoda",
                  params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=500000,
                  events = list(data = food_3_3.U))
  if(attributes(out_3_3U)$istate[1] != 2)(return(out_1_0))
  # Simulate dynamics for 4-4
  food_4_4 = periodic.starvation.events(10.74, 4, 14, 18)
  food_4_4.U <- subset(food_4_4, var != "P")
  out_4_4 <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                 initfunc = "initmod",  nout=1, outnames="Survival", method="lsoda",
                 params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=500000,
                 events = list(data = food_4_4))
  if(attributes(out_4_4)$istate[1] != 2)(return(out_1_0))
  out_4_4U <- ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_YVE", 
                  initfunc = "initmod",  nout=1, outnames="Survival", method="lsoda",
                  params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=500000,
                  events = list(data = food_4_4.U))
  if(attributes(out_4_4U)$istate[1] != 2)(return(out_1_0))
  
  
  result <- rbind(
    out_1_0, out_1_0, out_1_0, out_1_0, out_1_0, # 5
    out_2_2, out_2_2, out_2_2, # 3
    out_3_3, out_3_3, out_3_3, out_3_3, out_3_3, out_3_3, out_3_3, out_3_3, # 8    
    out_4_4, out_4_4, out_4_4, out_4_4,     
    
    out_1_0U, out_1_0U, out_1_0U, out_1_0U, out_1_0U,
    out_2_2U, out_2_2U, out_2_2U, out_2_2U, out_2_2U,
    out_3_3U, out_3_3U, out_3_3U, out_3_3U, out_3_3U,
    out_4_4U, out_4_4U, out_4_4U, out_4_4U, out_4_4U)
  
  result
  
}

extract.data<-function(data, w.t=7, start.age=0){
  ww<- which(data$time%%w.t==0 & data$time>=start.age)
  data[ww,]
}

## This function takes the simulator, params, inits, etc, and runs a
## forward simulation, and then extracts a subset of the points from
## the simulator (determined by w.t), without noise
make.states<-function(params, inits.F, inits.S, duration.R, duration.P, feeding.events.F, w.t=7){
  result.F = solve.DEB.supply(params, inits.F, duration.R, feeding.events.F)
  result.F = extract.data(result.F, start.age=28)
  Survival.F = result.F$Survival
  
  # This fix for survival data is currently specific to the sample size (96) and duration (32) of the experiment
  Survival.F[-((1:96)*32)] =   Survival.F[-((1:96)*32)] - Survival.F[-(1+(0:95)*32)]
  
  result.S = solve.DEB.starve(params, inits.S, duration.P)
  result.S = extract.data(result.S, start.age=14)
  Survival.S = result.S$Survival
  # This fix for survival data is currently specific to the sample size (40) and duration (19) of the experiment
  Survival.S[-((1:40)*19)] =   Survival.S[-((1:40)*19)] - Survival.S[-(1+(0:39)*19)]
  
  return(list(time=result.F$time, L=result.F$L, RH=result.F$RH, RP=result.F$RP,
              L2=result.S$L, W2=result.S$RP, E2=result.S$RH, SurvF=Survival.F, SurvS=Survival.S))
}

# 5 - Prior likelihood
prior.likelihood = function(x){
  prior.lik = with(as.list(x),
                   sum(dbeta(c(yPE, yEF, yEF2, yRP, yVE, mP, eh, k, fe), 1, 1, log=T)) + 
                     sum(dunif(c(sd.LI1, sd.LU1, sd.EI1, sd.EU1, sd.W1, sd.LI2, sd.LU2, sd.EI2, sd.EU2, sd.W2), min=0, max=10, log=T)) +
                     sum(dunif(c(ph, alpha, iPM, EM, DR, Fh, muD, kR, delta0, hdelta, theta, mR, hb), min=0, max=1000000, log=T)) +
                     dnorm(iM, mean=0.0183, sd=0.0016, log=T) + dnorm(M, mean=0.004, sd=0.00047, log=T) + dnorm(LM, mean=35, sd=2, log=T)
  )
  return(prior.lik)
}


# 6 - data likelihood
full.likelihood<-function(x){
  # simulate data
  sim.data = make.states(x, in.R, in.P, dur.R, dur.P, Feed.R, w.t=7)
  
  # data likelihood
  e.c<-1
  
  ## observation model
  gammaH<-0.015 # C content of eggs
  gammaP<-4e-5 # C content of cercs
  
  ## convert predictions into correct count units
  l.temp<-sim.data$L
  n.temp<-sim.data$RH/gammaH
  w.temp<-sim.data$RP/gammaP
  
  
  l2.temp<-sim.data$L2
  n2.temp<-sim.data$E2/gammaH
  w2.temp<-sim.data$W2/gammaP
  
  SF.temp<-sim.data$SurvF
  SS.temp<-sim.data$SurvS
  
  sd.LI1<-as.numeric(x["sd.LI1"])
  sd.LU1<-as.numeric(x["sd.LU1"])
  sd.EI1<-as.numeric(x["sd.EI1"])
  sd.EU1<-as.numeric(x["sd.EU1"])
  sd.W1<-as.numeric(x["sd.W1"])
  
  sd.LI2<-as.numeric(x["sd.LI2"])
  sd.LU2<-as.numeric(x["sd.LU2"])
  sd.EI2<-as.numeric(x["sd.EI2"])
  sd.EU2<-as.numeric(x["sd.EU2"])
  sd.W2<-as.numeric(x["sd.W2"])
  
  # Avoids simulations that fell short
  NObs = length(data$L)
  NObs2 = length(data$L2)
  if(length(n.temp) != NObs){print("Simulation too short"); return(-1e6)}
  if(length(n2.temp) != NObs2){print("Simulation too short"); return(-1e6)}    
  
  if(anyNA(l.temp)){print("NaNs in l.temp"); return(-1e6)}
  if(anyNA(l2.temp)){print("NaNs in l2.temp"); return(-1e6)}
  if(anyNA(SF.temp)){print("NaNs in SF.temp");return(-1e6)}
  if(anyNA(SS.temp)){print("NaNs in SS.temp"); return(-1e6)}
  
  sds = c(rep(sd.LI1, times = 2112), rep(sd.LU1, times=960))
  sd.Eggs = c(rep(sd.EI1, times = 2112), rep(sd.EU1, times=960))
  
  ## likelihoods from food gradient
  llik.L<- sum(dnorm(log(data$L), mean=log(l.temp), sd=sds, log=TRUE), na.rm=T)
  llik.Negg<- sum(dnorm(log(data$Negg+e.c), mean=log(n.temp+e.c), sd=sd.Eggs, log=TRUE), na.rm=T)
  llik.Nworms<- sum(dnorm(log(data$Nworms+e.c), mean=log(w.temp+e.c), sd=sd.W1, log=TRUE), na.rm=T)
  SF =SF.temp[which(data$Alive == 1)]
  llik.Survival <- sum(log(SF))
  
  sds2 = c(rep(sd.LI2, times = 380), rep(sd.LU2, times=380))
  sd.Eggs2 = c(rep(sd.EI2, times = 380), rep(sd.EU2, times=380))
  
  ## likelihoods from density gradient
  llik.L2<- sum(dnorm(log(data$L2), mean=log(l2.temp), sd=sds2, log=TRUE), na.rm=T)
  llik.Negg2<- sum(dnorm(log(data$Negg2+e.c), mean=log(n2.temp+e.c), sd=sd.Eggs2, log=TRUE), na.rm=T)
  llik.Nworms2<- sum(dnorm(log(data$Nworms2+e.c), mean=log(w2.temp+e.c), sd=sd.W2, log=TRUE), na.rm=T)
  
  
  SS =SS.temp[which(data$Alive2 == 1)]
  llik.Survival2 <- sum(log(SS)) 
  
  llik<-(llik.L + llik.Negg + llik.Nworms + llik.Survival + llik.L2 + llik.Negg2 + llik.Nworms2 + llik.Survival2)
  
  if(is.na(llik)|!is.finite(llik)){
    print("Infinite NLL")
    return(-1e6)}
  
  lprior = prior.likelihood(x)
  
  return(llik + lprior)
}


### Tuning ###

setwd("C:/RData")
samps = readRDS("Full_Fitting_ShiftingK.Rda")
pars = samps$samples[which.max(samps$log.p),]
pars = as.numeric(pars[1:35])
names(pars) = c("iM", "k", "M", "EM", "Fh", "muD", "DR", "fe", "yRP",
                "ph", "yPE", "iPM", "eh", "mP", "alpha", "yEF", "LM",
                "kR", "delta0", "hdelta", "hb", "theta", "mR", "yVE", "yEF2",
                "sd.LI1", "sd.LU1", "sd.EI1", "sd.EU1", "sd.W1", 
                "sd.LI2", "sd.LU2", "sd.EI2", "sd.EU2", "sd.W2")

variances = samps$cov.jump


### running the mcmc ###
test = MCMC(full.likelihood, init=pars, scale=as.matrix(variances), adapt=50, acc.rate = 0.3, n=200)

### converting to coda
testc = convert.to.coda(test)
testc = cbind(testc, "lpost" = test$log.p)
round(1 - rejectionRate(mcmc(testc)), 3)
plot(1:length(testc[,"lpost"]), testc[,"lpost"])
saveRDS(test, file="Full_Fitting_ShiftingK.Rda")
