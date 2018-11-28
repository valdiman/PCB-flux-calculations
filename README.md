# Individual PCB flux calculations IHSC 2017
R code to calculate the air-water flux of individual PCB congeners from Indiana Harbor and Ship Canal from air and water passive samples from 2017. It includes a Monte Carlo simulation (rep = 1000).

# set working directory on our computer
setwd("Z:.../Flux")

# Start functions

final.result = function(MW.PCB, H0.mean, H0.error, 
         C.PCB.water.mean, C.PCB.water.error, C.PCB.air.mean, C.PCB.air.error, nOrtho.Cl)
{

## fixed parameters

R = 8.3144
T = 298.15
w = 3

## Monte Carlo simulation

F.PCB.aw <- NULL
for (replication in 1:1000)
{

# random parameters

a <- rnorm(1, 0.085, 0.007)
b <- rnorm(1, 1, 0.5)
c <- rnorm(1, 32.7, 1.6)
H0 <- rnorm(1, H0.mean, H0.error) # should be normalized distribution

## Specific condition for each deployment
P <- rnorm(1, P.mean, P.error)					
u <- abs(rnorm(1, u.mean, u.error)) #m/s
C.PCB.water <- rnorm(1, C.PCB.water.mean, C.PCB.water.error) #ng/L
C.PCB.air <- rnorm(1, C.PCB.air.mean, C.PCB.air.error) #ng/m3
T.water <- rnorm(1, T.water.mean, T.water.error) #C 
T.air <- rnorm(1, T.air.mean, T.air.error) #C
Q <- abs(rnorm(1, Q.mean, Q.error)) #m/s3
h <- abs(rnorm(1, h.mean, h.error)) #m

## Equations

DeltaUaw <- (a*MW.PCB-b*nOrtho.Cl+c)*1000
K <- 10^(H0)*101325/(R*T)
K.air.water <- K*exp(-DeltaUaw/R*(1/(T.water+273.15)-1/T))
K.final <- K.air.water*(T.water+273.15)/(T.air+273.15) # no units
	
D.water.air <- 10^(-3)*1013.25*((273.15+T.air)^1.75*((1/28.97)+(1/18.0152))^(0.5))/P/(20.1^(1/3)+9.5^(1/3))^2
D.PCB.air <- D.water.air*(MW.PCB/18.0152)^(-0.5)
V.water.air <- 0.2*u +0.3 #cm/s eq. 20-15
V.PCB.air <- V.water.air*(D.PCB.air/D.water.air)^(2/3) #cm/s
	
visc.water <- 10^(-4.5318-220.57/(149.39-(273.15+T.water)))
dens.water <- (999.83952+16.945176*T.water-7.9870401*10^-3*T.water^2-46.170461*10^-6*3+105.56302*10^-9*T.water^4-280.54253*10^-12*T.water^5)/(1+16.87985*10^-3*T.water)
v.water <- visc.water/dens.water*10000
diff.co2 <- 0.05019*exp(-19.51*1000/(273.15+T.water)/R)
D.PCB.water <- diff.co2*(MW.PCB/44.0094)^(-0.5)
Sc.PCB.water <- v.water/D.PCB.water
Sc.co2.water <- v.water/diff.co2
Flow.veloc <- Q*0.02/(w*h)
k600 <- 1.72*(Flow.veloc*100/h)^(0.5) #cm/h
V.PCB.water = k600/3600*(Sc.PCB.water/Sc.co2.water)^(-0.5) #cm/s

mtc.PCB <- ((1/V.PCB.water+1/(V.PCB.air*K.final)))^(-1) #cm/s

F.PCB.aw <- c(F.PCB.aw, mtc.PCB*(C.PCB.water-C.PCB.air/K.final/1000)*3600*24/100) #ug/m2/d

}

mmm <- mean(F.PCB.aw)	#ng/m2/day
sss <- sd(F.PCB.aw)	#ng/m2/day
q2.5 <- quantile(F.PCB.aw, 0.025)
q97.5 <- quantile(F.PCB.aw, 0.975)

c(mmm, sss, q2.5, q97.5)
}

# Individual chemical properties, meteorological and environmental conditions

## Chemical properties (e.g., deployment 1)

pars <- read.csv("Variables/Variables congenerD1.csv")
Congener <- pars$Congener
MW.PCB <- pars$MW.PCB
H0.mean <- pars$H0
H0.error <- pars$X
C.PCB.water.mean <- pars$C.PCB.water #ng/L
C.PCB.water.error <- pars$X.1
C.PCB.air.mean <- pars$C.PCB.air #ng/m3
C.PCB.air.error <- pars$X.2
nOrtho.Cl <- pars$nOrtho.Cl

 ## Meteorological and environmental conditions

parsM <- read.csv("Meteo_param.csv") #change number before braquets, for different deployment!
T.air.mean <- parsM$D1[1]
T.air.error <- parsM$X.1[1]
T.water.mean <- parsM$D1[2]
T.water.error <- parsM$X.1[2]
u.mean <- parsM$D1[3]
u.error <- parsM$X.1[3]
P.mean <- parsM$D1[4]
P.error <- parsM$X.1[4]
Q.mean <- parsM$D1[5]
Q.error <- parsM$X.1[5]
h.mean <- parsM$D1[6]
h.error <- parsM$X.1[6]

# Results

Num.Congener <- length(Congener)

result <- NULL
for (i in 1:Num.Congener)
{
	result <- rbind(result, final.result(MW.PCB[i], H0.mean[i], H0.error[i], 
         C.PCB.water.mean[i], C.PCB.water.error[i], C.PCB.air.mean[i], C.PCB.air.error[i], nOrtho.Cl[i]))
}

final.result = data.frame(Congener, result)
names(final.result) = c("Congener", "Mean", "Std", "2.5%CL", "97.5%CL")

## Example for deployment 1

write.csv(final.result, row.names=F, file="Results/Net/fluxD1.csv")

 # References:
  # Paper
  
  Martinez, Awad, Herkert, Hornbuckle (2019) Determination of PCB fluxes from Indiana Harbor and Ship Canal using dual-deployed air and water passive samplers Environmental Pollution 244, 469-476, https://doi.org/10.1016/j.envpol.2018.10.048
  
  # Available data:
  
  Martinez, Andres; Awad, Andrew M; Herkert, Nicholas J; Hornbuckle, Keri C (2018): PCB congener data of gas-phase, freely-dissolved water, air-water fugacity ratios and air-water fluxes in Indiana Harbor and Ship Canal, IN, USA. PANGAEA, https://doi.org/10.1594/PANGAEA.894908
  

