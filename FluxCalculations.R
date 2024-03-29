# Introduction ------------------------------------------------------------
# This code describes the individual PCB air-water flux calculation for one only
# deployment time.
# Change should be done in section Read Pangaea datasets
# for the others deployment times

# Read Pangaea datasets ---------------------------------------------------

# https://cran.microsoft.com/snapshot/2022-01-01/web/packages/pangaear/pangaear.pdf
# Install packages
install.packages('pangaear')

# Libraries
library(pangaear)

# set a different cache path from the default
pg_cache$cache_path_set(full_path = "/Users/andres/OneDrive - University of Iowa/work/ISRP/Project4/Old/Codes/PCBFluxesIHSC/PCBFluxesIHSC")

# Download the datasets from Pangaea
data.water <- pg_data(doi = '10.1594/PANGAEA.894906')
data.air <- pg_data(doi = '10.1594/PANGAEA.894905')
data.meteo <- pg_data(doi = '10.1594/PANGAEA.894919')

# Obtain just concentrations from Pangaea dataset
pars.water <- data.water[[1]]$data # pg/L
pars.air <- data.air[[1]]$data # ng/m3

# Extract first deployment time data
# i.e., 2016-11-23 to 2017-01-24
pars.water.2 <- pars.water[1,]
pars.air.2 <- pars.air[1,]

# Install package to work with data
install.packages('tidyverse')

# Library
library(tidyverse)

# Remove individual standard deviation concentrations
pars.water.3 <- select(pars.water.2, -contains("std"))
pars.air.3 <- select(pars.air.2, -contains("std"))

# Remove metadata
pars.water.4 <- subset(pars.water.3,
                       select = -c(`Method comm (Values = 0 indicates non-dete...)`:`Date/Time (end)`))
pars.air.4 <- subset(pars.air.3,
                     select = -c(`Method comm (Values = 0 indicates non-dete...)`:`Date/Time (end)`))

# Obtain just meteorological parameters from Pangaea dataset
pars.meteo <- data.meteo[[1]]$data # ng/m3

# Extract first deployment time data
# i.e., 2016-11-23 to 2017-01-24
pars.meteo.2 <- pars.meteo[1,]

# Create P-C properties matrix --------------------------------------------

# Create matrix to storage P-C data
pars <- data.frame(matrix(NA, nrow = 171, ncol = 5))

# Add column names
colnames(pars) <- c('Congener', 'MW.PCB', 'nOrtho.Cl', 'H0.mean', "H0.error")

# Add PCB names
pars[,1] <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
              '12+13', '15', '16', '17', '18+30', '19', '20+28', '21+33',
              '22', '23', '24', '25', '26+29', '27', '31', '32', '34',
              '35', '36', '37', '38', '39', '40+71', '41', '42', '43',
              '44+47+65', '45', '46', '48', '49+69', '50+53', '51', '52',
              '54', '55', '56', '57', '58', '59+62+75', '60', '61+70+74+76',
              '63', '64', '66', '67', '68', '72', '73', '77', '78',
              '79', '80', '81', '82', '83', '84', '85+116', '86+97+109+119',
              '88', '89', '90+101+113', '91', '92', '93+100', '94', '95',
              '96', '87+125', '98', '99', '102', '103', '104', '105',
              '106', '108', '107+124', '110', '111', '112', '114', '115',
              '117', '118', '120', '121', '122', '123', '126', '127',
              '129+138+163', '130', '131', '132', '133', '134', '135+151',
              '136', '137', '139+140', '141', '143', '142', '144', '145',
              '146', '147+149', '148', '150', '152', '153+168', '154', '155',
              '156+157', '158', '159', '160', '161', '162', '164', '165',
              '167', '169', '170', '171+173', '172', '174', '175', '176', '177',
              '178', '179', '180+193', '181', '182', '183', '184', '185',
              '186', '187', '188', '189', '190', '191', '192', '194', '195',
              '196', '197', '198+199', '200', '201', '202', '203', '205',
              '206', '207', '208', '209')

# Add MW of individual PCB congeners
pars[1:3,2] <- c(188.644)
pars[4:13,2] <- c(223.088)
pars[14:33,2] <- c(257.532)
pars[34:65,2] <- c(291.976)
pars[66:102,2] <- c(326.42)
pars[103:135,2] <- c(360.864)
pars[136:157,2] <- c(395.308)
pars[158:167,2] <- c(429.752)
pars[168:170,2] <- c(465.740544)
pars[171,2] <- c(498.64)

# Add ortho Cl of individual PCB congeners
pars[,3] <- c(1,	0,	0,	2,	1,	1,	1,	1,	1,	2,	0,	0,	0,	2,
              2,	2,	3,	1,	1,	1,	1,	2,	1,	1,	2,	1,	2,	1,
              0,	0,	0,	0,	0,	2,	2,	2,	2,	2,	3,	3,	2,	2,
              3,	3,	2,	4,	1,	1,	1,	1,	2,	1,	1,	1,	2,	1,
              1,	1,	1,	2,	0,	0,	0,	0,	0,	2,	2,	3,	2,	2,
              3,	3,	2,	3,	2,	3,	3,	3,	4,	2,	3,	2,	3,	3,
              4,	1,	2,	1,	1,	2,	1,	2,	1,	2,	2,	1,	1,	2,
              1,	1,	0,	0,	2,	2,	3,	3,	2,	3,	3,	4,	2,	3,
              2,	3,	3,	3,	4,	2,	3,	3,	4,	4,	2,	3,	4,	1,
              2,	1,	2,	2,	1,	2,	2,	1,	0,	2,	3,	2,	3,	3,
              4,	3,	3,	4,	2,	3,	3,	3,	4,	3,	4,	3,	4,	1,
              2,	2,	2,	1,	3,	3,	4,	3,	4,	4,	4,	3,	2,	3,
              4,	4,	4)

# Add Ho of individual PCB congeners
pars[,4] <- c(-3.526,	-3.544,	-3.562,	-3.483,	-3.622,	-3.486,	-3.424,
              -3.518,	-3.490,	-3.373,	-3.537,	-3.595,	-3.649,	-3.600,
              -3.428,	-3.495,	-3.355,	-3.544,	-3.62,	-3.719,	-3.497,
              -3.500,	-3.500,	-3.526,	-3.393,	-3.562,	-3.407,	-3.375,
              -3.3745, -3.473, -3.818,	-3.634,	-3.524,	-3.503,	-3.612,
              -3.592, -3.475,	-3.638,	-3.450,	-3.470,	-3.519,	-3.452,
              -3.366, -3.292,	-3.496,	-3.242,	-3.739,	-3.820,	-3.568,
              -3.602, -3.517,	-3.816,	-3.694,	-3.615,	-3.565,	-3.693,
              -3.631, -3.424,	-3.441,	-3.284,	-3.989,	-3.787,	-3.705,
              -3.426, -3.844,	-3.835,	-3.674,	-3.600,	-3.716,	-3.745,
              -3.415, -3.526,	-3.610,	-3.461,	-3.585,	-3.468,	-3.407,
              -3.523, -3.387,	-3.745,	-3.407,	-3.603,	-3.431,	-3.298,
              -3.130, -4.003,	-3.783,	-3.798,	-3.768,	-3.707,	-3.574,
              -3.574, -3.845,	-3.610,	-3.618,	-3.901,	-3.610,	-3.253,
              -3.901, -3.759,	-4.087,	-3.807,	-3.886,	-3.817,	-3.616,
              -3.693, -3.691,	-3.639,	-3.548,	-3.492,	-3.731,	-3.483,
              -3.760, -3.502,	-3.531,	-3.529,	-3.328,	-3.727,	-3.625,
              -3.367, -3.296,	-3.369,	-3.783,	-3.418,	-3.075,	-4.053,
              -3.782, -3.808,	-3.670,	-3.545,	-3.881,	-3.754,	-3.560,
              -3.959, -4.186,	-4.059,	-3.763,	-3.924,	-3.772,	-3.651,
              -3.527, -3.787,	-3.671,	-3.560,	-3.969,	-3.638,	-3.590,
              -3.696, -3.339,	-3.669,	-3.434,	-3.693,	-3.353,	-3.177,
              -3.950, -3.876,	-3.718,	-4.174,	-3.926,	-3.884,	-3.596,
              -3.644, -3.619,	-3.884,	-3.651,	-3.853,	-4.059,	-4.059,
              -3.772, -3.777,	-3.948)

# Add Ho error
pars[,5] <- c(0.662)

# Adjust names
Congener <- pars$Congener
MW.PCB <- pars$MW.PCB
H0.mean <- pars$H0.mean
H0.error <- pars$H0.error
nOrtho.Cl <- pars$nOrtho.Cl

# Flux calculations -------------------------------------------------------

final.result = function(MW.PCB, H0.mean, H0.error, 
                        C.PCB.water.mean, C.PCB.water.error,
                        C.PCB.air.mean, C.PCB.air.error, nOrtho.Cl) {
  # fixed parameters
  
  R = 8.3144
  T = 298.15
  w = 3
  
  F.PCB.aw <- NULL
  for (replication in 1:1000) {
    
    # random parameters
    
    a <- rnorm(1, 0.085, 0.007)
    b <- rnorm(1, 1, 0.5)
    c <- rnorm(1, 32.7, 1.6)
    H0 <- rnorm(1, H0.mean, H0.error)
    P <- rnorm(1, P.mean, P.error)					
    u <- abs(rnorm(1, u.mean, u.error)) #m/s
    C.PCB.water <- rnorm(1, C.PCB.water.mean, C.PCB.water.error) #ng/L
    C.PCB.air <- rnorm(1, C.PCB.air.mean, C.PCB.air.error) #ng/m3
    T.water <- rnorm(1, T.water.mean, T.water.error) #C 
    T.air <- rnorm(1, T.air.mean, T.air.error) #C
    Q <- abs(rnorm(1, Q.mean, Q.error))
    h <- abs(rnorm(1, h.mean, h.error)) #m
    
    # computed values
    
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
    
    # Two methods to determine k600. K600 = f(flow or wind).
    # For the particular conditions of the sampling location,
    # a better estimate is provided by the flow and not the wind speed
    # see k600 script
    
    # i) flow
    
    Sc.co2.water <- v.water/diff.co2
    Flow.veloc <- Q*0.02/(w*h)
    k600 <- 1.72*(Flow.veloc*100/h)^(0.5) #cm/h
    V.PCB.water = k600/3600*(Sc.PCB.water/Sc.co2.water)^(-0.5) #cm/s
    mtc.PCB <- ((1/V.PCB.water+1/(V.PCB.air*K.final)))^(-1) #cm/s
    F.PCB.aw <- c(F.PCB.aw, mtc.PCB*(C.PCB.water-C.PCB.air/K.final/1000)*3600*24*10) #pg/m2/d
    
  }
  
  mmm <- mean(F.PCB.aw)	#pg/m2/day
  sss <- sd(F.PCB.aw)	#pg/m2/day
  q2.5 <- quantile(F.PCB.aw, 0.025) #pg/m2/day
  q97.5 <- quantile(F.PCB.aw, 0.975) #pg/m2/day
  
  c(mmm, sss, q2.5, q97.5)
}

C.PCB.water.mean <- as.numeric(pars.water.4/1000) # 1000 to get ng/L from pg/L
C.PCB.water.error <- as.numeric(pars.water.4*0.2/1000) # 20% error
C.PCB.air.mean <- as.numeric(pars.air.4) # ng/m3
C.PCB.air.error <- as.numeric(pars.air.4*0.2) # 20% error
T.air.mean <- pars.meteo.2$`TTT [°C] (average)`
T.air.error <- pars.meteo.2$`TTT std dev [±]`
T.water.mean <- pars.meteo.2$`Temp [°C] (average)`
T.water.error <- pars.meteo.2$`Temp std dev [±]`
u.mean <- pars.meteo.2$`ff [m/s] (average)`
u.error <- pars.meteo.2$`ff std [±]`
P.mean <- pars.meteo.2$`PPPP [hPa]`
P.error <- pars.meteo.2$`PPPP std [±]`
Q.mean <- pars.meteo.2$`Q [m**3/s] (average)`
Q.error <- pars.meteo.2$`Q std dev [±]`
h.mean <- pars.meteo.2$`Depth water [m] (average)`
h.error <- pars.meteo.2$`Depth water std dev [±]`

Num.Congener <- length(Congener)

result <- NULL
for (i in 1:Num.Congener) {
  result <- rbind(result,
                  final.result(MW.PCB[i], H0.mean[i], H0.error[i],
                               C.PCB.water.mean[i], C.PCB.water.error[i],
                               C.PCB.air.mean[i], C.PCB.air.error[i], nOrtho.Cl[i]))
}

final.result = data.frame(Congener, result)
names(final.result) = c("Congener", "Mean (pg/m2/d)", "Std (pg/m2/d)",
                        "2.5%CL (pg/m2/d)", "97.5%CL (pg/m2/d)")
