
# Set working directory
#setwd("~/Box/Jeff Work/Jeff work/")

# load libraries that are needed
pacman::p_load('gss', 'copula', 'matrixcalc', 'orthopolynom', 
               'foreign', 'logspline', 'readstata13', 'plm', 
               'doBy', 'dplyr', 'usmap', 'viridis', 'ggplot2', 'purrr')


# additional source files that Ximing programmed
source('~/Box/Jeff Work/Jeff work/ximing/codeXiming/legendre.r')
source('~/Box/Jeff Work/Jeff work/ximing/codeXiming/svd_inv.r')
source('~/Box/Jeff Work/Jeff work/ximing/codeXiming/copula_tri4.r')
source('~/Box/Jeff Work/Jeff work/ximing/codeXiming/get_copula2.r')
source('~/Box/Jeff Work/Jeff work/ximing/codeXiming/pk3.r')
source('~/Box/Jeff Work/Jeff work/ximing/codeXiming/lp0.r')
source('~/Box/Jeff Work/Jeff work/ximing/codeXiming/lp1.r')
source('~/Box/Jeff Work/Jeff work/ximing/codeXiming/lp2.r')
source('~/Box/Jeff Work/Jeff work/ximing/codeXiming/lp3.r')

#-------------------------------------------------------------------------------------------------------------------------------;
# load data set
#-------------------------------------------------------------------------------------------------------------------------------;
      
dat=read.dta13('~/Box/Jeff Work/Jeff work/ximing/dataSTATA/corn/timeCountyFE/trimester2.dta')

#futdat7099 <- readRDS("~/Box/Jeff Work/futureClimate/residuals85.Rdata") %>% filter(year %in% 2070:2099)
#fut45dat3565 <- readRDS("~/Box/Jeff Work/futureClimate/residuals45.Rdata") %>% filter(year %in% 2035:2065)
#fut45dat7099 <- readRDS("~/Box/Jeff Work/futureClimate/residuals45.Rdata") %>% filter(year %in% 2070:2099)
dat8005 <- read.dta13('~/Box/Jeff Work/Jeff work/ximing/dataSTATA/corn/timeCountyFE/trimester2.dta') %>% filter(year %in% 1980:2005)

res.dat = cbind(dat$logYieldRes,dat$dday29CRes,dat$precRes, dat$year)
res.dat8005 <- cbind(dat8005$dday29CRes,dat8005$precRes, dat8005$year)
    
        
        #-------------------------------------------------------------------------------------------------------------------------------;
        # estimate model;
        #-------------------------------------------------------------------------------------------------------------------------------;
        res.fit = copula_tri4(res.dat,gam=.0001,r=6,m=40,tt=2)
        cbind(res.fit$moment,res.fit$m_hat)
        
        # estimate the marginal using log-spline
        res.fitMarginal = logspline(res.dat[,1], lbound=min(res.dat[,1])-0.1*(max(res.dat[,1])-min(res.dat[,1])), ubound=max(res.dat[,1])+0.1*(max(res.dat[,1])-min(res.dat[,1])))
        
        # calculate the conditional mean on grids of temperature and prec
        # get quantiles ranging from 2% to 98% in 2% steps, i.e., 49 values
        res.q.temp = quantile(res.dat[,2],prob=seq(.02,.98,.02))
        res.q.prec = quantile(res.dat[,3],prob=seq(.02,.98,.02))
        res.q.grid = expand.grid(res.q.temp,res.q.prec)
        res.q.gridID = expand.grid(seq(.02,.98,.02),seq(.02,.98,.02))
        
        # initalize matrix to save results with five columns
        res.out = matrix(0,nrow=dim(res.q.grid)[1],ncol=9)
        res.nObs = dim(res.dat)[1]
        
        # Gaussian quadrature
        res.y0 = gauss.quad(1000,c(min(res.dat[,1])-0.1*(max(res.dat[,1])-min(res.dat[,1])), max(res.dat[,1])+0.1*(max(res.dat[,1])-min(res.dat[,1])) ))
        res.dy0 = dlogspline(res.y0$pt,res.fitMarginal)
        res.py0 = res.y0$pt
        for (i in 1:length(res.y0$pt)) {
          res.py0[i]=sum(res.y0$pt[i]>=res.dat[,1])/(res.nObs+1)
        }
        
        res.py0b = plogspline(res.y0$pt,res.fitMarginal)	
        res.q.yield5 = quantile(res.dat[,1],.05)
        res.q.M5 = mean(res.dat[,1])-0.05
        res.q.M10 = mean(res.dat[,1])-0.1
        res.q.M20 = mean(res.dat[,1])-0.2
        res.q.M33 = mean(res.dat[,1])-1/3
        res.q.se = mean(res.dat[,1])-sd(res.dat[,1])
        #-------------------------------------------------------------------------------------------------------------------------------;
        # get results;
        # recall about conditional density f(y|x) = f(x,y) / fx(x)
        #        where f(x,y) is joint density
        #              fx(x) is marginal density
        #-------------------------------------------------------------------------------------------------------------------------------;
        for (i in 1:dim(res.out)[1]) {
          temp = cbind(res.py0, res.q.gridID[i,1], res.q.gridID[i,2])
          temp2 = get_copula2(temp,res.fit,pp=2)
          
          # integral of conditional density;
          res.out[i,1] = sum(temp2 * res.dy0 * res.y0$wt)
          
          # conditional mean (not normalized by integral of conditional density res.out[,1];
          res.out[i,2] = sum(res.y0$pt * temp2 * res.dy0 * res.y0$wt)
          
          # conditional second moment (not normalized yet); 
          res.out[i,3] = sum((res.y0$pt)^2 * temp2 * res.dy0 * res.y0$wt)
          
          # conditional third moment (not normalized yet);
          res.out[i,4] = sum((res.y0$pt)^3 * temp2 * res.dy0 * res.y0$wt)
          
          # for models with county fixed effects, log yields give directly yield shortfall
            res.out[i,5] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= -0.05))
            res.out[i,6] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= -0.1))
            res.out[i,7] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= -0.2))
            res.out[i,8] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= -1/3))
            res.out[i,9] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= res.q.se))
        } 

        
        res.z1  = res.out[,2] / res.out[,1]
        res.z2  = res.out[,3] / res.out[,1]
        res.zse = sqrt(res.z2 - res.z1^2)
        res.z3  = res.out[,4] / res.out[,1]
        res.z4  = res.out[,5] / res.out[,1]
        res.z5  = res.out[,6] / res.out[,1]
        res.z6  = res.out[,7] / res.out[,1]
        res.z7  = res.out[,8] / res.out[,1]
        res.z8 = res.out[,9] / res.out[,1]
        
#Subset for the earlier time period
modelmeans <- function(modelnum){
  
  futdat3565 <- readRDS("~/Box/Jeff Work/futureClimate/residuals85.Rdata") %>% filter(year %in% 2035:2065, modelNr == modelnum)
  res.datfut85_3565 = cbind(futdat3565$resid_temp,futdat3565$resid_prec, futdat3565$year)
  #res.datfut85_7099 = cbind(futdat7099$residdday29,futdat7099$residprec, futdat7099$year)
        
  #res.datfut45_3565 = cbind(fut45dat3565$residdday29,fut45dat3565$residprec, fut45dat3565$year)
  #res.datfut45_7099 = cbind(fut45dat7099$residdday29,fut45dat7099$residprec, fut45dat7099$year)
        
        
  index <- matrix(NA, nrow = nrow(res.datfut85_3565), ncol = 1)
        # fill the index by searching for closest temp/precip match for each observation in our data
        for (i in 1:nrow(res.datfut85_3565)){
          index[i] <- row.names(res.q.grid[which(abs(res.q.grid[,2] - res.datfut85_3565[i,2]) == min(abs(res.q.grid[,2] - res.datfut85_3565[i,2])) & abs(res.q.grid[,1] - res.datfut85_3565[i,1])==min(abs(res.q.grid[,1] - res.datfut85_3565[i,1]))),])
        }
  
        Pi <- matrix(NA, nrow = nrow(res.datfut85_3565), ncol=1)
        pred_shortfall <- as.data.frame(res.out)[,6]
        
        for (i in 1:nrow(index)) {
          j <- as.numeric(index[i])
          Pi[i] <- pred_shortfall[j]
        }
        
        Pi <- as_data_frame(Pi)
        Pi$modelNr <- modelnum
        return(Pi)
        #saveRDS(Pi, paste0("futPredShortfall", modelnum, ".Rdata"))
        #rm(Pi)
        #rm(res.datfut85_3565)
        #rm(pred_shortfall)
}  

Pi <- map_dfr(1:21, modelmeans)

#We also want the conditional mean for each year/county combo
#Right now, we have the probability that a yield is below the average given a particular combination of 
#temperature and precipitation
#I think I know how to pull the conditional yield for a given scenario
#Just call the information that corresponds to the closest temperature/precip combo
#But how do I translate that into a new yield shortfall scenario?
#Do I re-run the model using the predicted yield data from the historical time period to obtain the 
#yield shortfall?
# I can't come up with a better solution, so I might try it that way. 

futdat3565 <- readRDS("~/Box/Jeff Work/futureClimate/residuals85.Rdata") %>% filter(year %in% 2035:2065)

#Add in fips and year from the future data (I need to do something different with the fips and year)
Pi$fips <- futdat3565$fips
Pi$year <- futdat3565$year


Pi2 <- Pi %>% group_by(fips, year) %>% summarise(modelMean = mean(V1))

        
 #Here's the end where I need to take the mean across all of these
        
        #Let's condense this information
        Pi3 <- Pi2 %>% group_by(fips) %>% summarise(meanProb = mean(modelMean)) 
        
        map_with_data(Pi3, values = "meanProb")
        df <- us_map(regions = "counties")
        
Plot <- plot_usmap(regions = "counties", 
                   include = c(.east_north_central, .east_south_central,
                              .mid_atlantic, .south_atlantic,
                              .west_north_central, .west_south_central),
                   data = Pi3, values = "meanProb",
                   lines = "grey30", labels = FALSE,
                   label_color = "black") + 
          ggtitle("Average Probability of a 10% Yield Shortfall, RCP 8.5 (2035-2065)") + 
          guides(fill=guide_legend(title="Probability")) + 
          scale_fill_gradientn(colours = viridis(8), na.value = "grey80", limits = c(0,0.8))
Plot  


#Run this again for 2070-2099 ----
modelmeans <- function(modelnum){
  
  futdat7099 <- readRDS("~/Box/Jeff Work/futureClimate/residuals85.Rdata") %>% filter(year %in% 2070:2099, modelNr == modelnum)
  res.datfut85_7099 = cbind(futdat7099$resid_temp,futdat7099$resid_prec, futdat7099$year)
  
  index <- matrix(NA, nrow = nrow(res.datfut85_7099), ncol = 1)
  # fill the index by searching for closest temp/precip match for each observation in our data
  for (i in 1:nrow(res.datfut85_7099)){
    index[i] <- row.names(res.q.grid[which(abs(res.q.grid[,2] - res.datfut85_7099[i,2]) == min(abs(res.q.grid[,2] - res.datfut85_7099[i,2])) & abs(res.q.grid[,1] - res.datfut85_7099[i,1])==min(abs(res.q.grid[,1] - res.datfut85_7099[i,1]))),])
  }
  Pi <- matrix(NA, nrow = nrow(res.datfut85_7099), ncol=1)
  pred_shortfall <- as.data.frame(res.out)[,6]
  
  for (i in 1:nrow(index)) {
    j <- as.numeric(index[i])
    Pi[i] <- pred_shortfall[j]
  }
  
  Pi <- as_data_frame(Pi)
  Pi$modelNr <- modelnum
  return(Pi)
  #saveRDS(Pi, paste0("futPredShortfall", modelnum, ".Rdata"))
  #rm(Pi)
  #rm(res.datfut85_3565)
  #rm(pred_shortfall)
}  

Pi <- map_dfr(1:21, modelmeans)

futdat7099 <- readRDS("~/Box/Jeff Work/futureClimate/residuals85.Rdata") %>% filter(year %in% 2070:2099)

#Add in fips and year from the future data (I need to do something different with the fips and year)
Pi$fips <- futdat7099$fips
Pi$year <- futdat7099$year


Pi2 <- Pi %>% group_by(fips, year) %>% summarise(modelMean = mean(V1))


#Here's the end where I need to take the mean across all of these

#Let's condense this information
Pi3 <- Pi2 %>% group_by(fips) %>% summarise(meanProb = mean(modelMean)) 

map_with_data(Pi3, values = "meanProb")
df <- us_map(regions = "counties")

Plot <- plot_usmap(regions = "counties", 
                   include = c(.east_north_central, .east_south_central,
                               .mid_atlantic, .south_atlantic,
                               .west_north_central, .west_south_central),
                   data = Pi3, values = "meanProb",
                   lines = "grey30", labels = FALSE,
                   label_color = "black") + 
  ggtitle("Average Probability of a 10% Yield Shortfall, RCP 8.5 (2070-2099)") + 
  guides(fill=guide_legend(title="Probability")) + 
  scale_fill_gradientn(colours = viridis(8), na.value = "grey80", limits = c(0,0.8))
Plot  


        
#Let's try this again with RCP 4.5


modelmeans <- function(modelnum){
  
  futdat3565 <- readRDS("~/Box/Jeff Work/futureClimate/residuals45.Rdata") %>% filter(year %in% 2035:2065, modelNr == modelnum)
  res.datfut45_3565 = cbind(futdat3565$resid_temp,futdat3565$resid_prec, futdat3565$year)
  
  
  index <- matrix(NA, nrow = nrow(res.datfut45_3565), ncol = 1)
  # fill the index by searching for closest temp/precip match for each observation in our data
  for (i in 1:nrow(res.datfut45_3565)){
    index[i] <- row.names(res.q.grid[which(abs(res.q.grid[,2] - res.datfut45_3565[i,2]) == min(abs(res.q.grid[,2] - res.datfut45_3565[i,2])) & abs(res.q.grid[,1] - res.datfut45_3565[i,1])==min(abs(res.q.grid[,1] - res.datfut45_3565[i,1]))),])
  }
  Pi <- matrix(NA, nrow = nrow(res.datfut45_3565), ncol=1)
  pred_shortfall <- as.data.frame(res.out)[,6]
  
  for (i in 1:nrow(index)) {
    j <- as.numeric(index[i])
    Pi[i] <- pred_shortfall[j]
  }
  
  Pi <- as_data_frame(Pi)
  Pi$modelNr <- modelnum
  return(Pi)
  #saveRDS(Pi, paste0("futPredShortfall", modelnum, ".Rdata"))
  #rm(Pi)
  #rm(res.datfut85_3565)
  #rm(pred_shortfall)
}  

Pi <- map_dfr(1:21, modelmeans)

futdat3565 <- readRDS("~/Box/Jeff Work/futureClimate/residuals45.Rdata") %>% filter(year %in% 2035:2065)

#Add in fips and year from the future data (I need to do something different with the fips and year)
Pi$fips <- futdat3565$fips
Pi$year <- futdat3565$year


Pi2 <- Pi %>% group_by(fips, year) %>% summarise(modelMean = mean(V1))


#Here's the end where I need to take the mean across all of these

#Let's condense this information
Pi3 <- Pi2 %>% group_by(fips) %>% summarise(meanProb = mean(modelMean)) 

map_with_data(Pi3, values = "meanProb")
df <- us_map(regions = "counties")

Plot <- plot_usmap(regions = "counties", 
                   include = c(.east_north_central, .east_south_central,
                               .mid_atlantic, .south_atlantic,
                               .west_north_central, .west_south_central),
                   data = Pi3, values = "meanProb",
                   lines = "grey30", labels = FALSE,
                   label_color = "black") + 
  ggtitle("Average Probability of a 10% Yield Shortfall, RCP 4.5 (2035-2065)") + 
  guides(fill=guide_legend(title="Probability")) + 
  scale_fill_gradientn(colours = viridis(8), na.value = "grey80", limits = c(0,0.8))
Plot  


#Run this again for 2070-2099 ----
modelmeans <- function(modelnum){
  
  futdat7099 <- readRDS("~/Box/Jeff Work/futureClimate/residuals45.Rdata") %>% filter(year %in% 2070:2099, modelNr == modelnum)
  res.datfut45_7099 = cbind(futdat7099$resid_temp,futdat7099$resid_prec, futdat7099$year)
  
  index <- matrix(NA, nrow = nrow(res.datfut45_7099), ncol = 1)
  # fill the index by searching for closest temp/precip match for each observation in our data
  for (i in 1:nrow(res.datfut45_7099)){
    index[i] <- row.names(res.q.grid[which(abs(res.q.grid[,2] - res.datfut45_7099[i,2]) == min(abs(res.q.grid[,2] - res.datfut45_7099[i,2])) & abs(res.q.grid[,1] - res.datfut45_7099[i,1])==min(abs(res.q.grid[,1] - res.datfut45_7099[i,1]))),])
  }
  Pi <- matrix(NA, nrow = nrow(res.datfut45_7099), ncol=1)
  pred_shortfall <- as.data.frame(res.out)[,6]
  
  for (i in 1:nrow(index)) {
    j <- as.numeric(index[i])
    Pi[i] <- pred_shortfall[j]
  }
  
  Pi <- as_data_frame(Pi)
  Pi$modelNr <- modelnum
  return(Pi)
  #saveRDS(Pi, paste0("futPredShortfall", modelnum, ".Rdata"))
  #rm(Pi)
  #rm(res.datfut85_3565)
  #rm(pred_shortfall)
}  

Pi <- map_dfr(1:21, modelmeans)

futdat7099 <- readRDS("~/Box/Jeff Work/futureClimate/residuals45.Rdata") %>% filter(year %in% 2070:2099)

#Add in fips and year from the future data (I need to do something different with the fips and year)
Pi$fips <- futdat7099$fips
Pi$year <- futdat7099$year


Pi2 <- Pi %>% group_by(fips, year) %>% summarise(modelMean = mean(V1))


#Here's the end where I need to take the mean across all of these

#Let's condense this information
Pi3 <- Pi2 %>% group_by(fips) %>% summarise(meanProb = mean(modelMean)) 

map_with_data(Pi3, values = "meanProb")
df <- us_map(regions = "counties")

Plot <- plot_usmap(regions = "counties", 
                   include = c(.east_north_central, .east_south_central,
                               .mid_atlantic, .south_atlantic,
                               .west_north_central, .west_south_central),
                   data = Pi3, values = "meanProb",
                   lines = "grey30", labels = FALSE,
                   label_color = "black") + 
  ggtitle("Average Probability of a 10% Yield Shortfall, RCP 4.5 (2070-2099)") + 
  guides(fill=guide_legend(title="Probability")) + 
  scale_fill_gradientn(colours = viridis(8), na.value = "grey80", limits = c(0,0.8))
Plot  




#Run this again for 1980-2005----

index <- matrix(NA, nrow = nrow(res.dat8005), ncol = 1)
# fill the index by searching for closest temp/precip match for each observation in our data
for (i in 1:nrow(res.dat8005)){
  index[i] <- row.names(res.q.grid[which(abs(res.q.grid[,2] - res.dat8005[i,2]) == min(abs(res.q.grid[,2] - res.dat8005[i,2])) & abs(res.q.grid[,1] - res.dat8005[i,1])==min(abs(res.q.grid[,1] - res.dat8005[i,1]))),])
}
Pi <- matrix(NA, nrow = nrow(res.dat8005), ncol=1)
pred_shortfall <- as.data.frame(res.out)[,6]

for (i in 1:nrow(index)) {
  j <- as.numeric(index[i])
  Pi[i] <- pred_shortfall[j]
}

#Add in fips and year from the future data
Pi <- as_data_frame(Pi)
Pi$fips <- dat8005$fips
Pi$year <- dat8005$year


#Let's condense this information
Pi2 <- Pi %>% group_by(fips) %>% summarise(meanProb = mean(V1))


map_with_data(Pi2, values = "meanProb")
df <- us_map(regions = "counties")

Plot <- plot_usmap(regions = "counties", 
                   include = c(.east_north_central, .east_south_central,
                               .mid_atlantic, .south_atlantic,
                               .west_north_central, .west_south_central),
                   data = Pi2, values = "meanProb",
                   lines = "grey30", labels = FALSE,
                   label_color = "black") + 
  ggtitle("Average Probability of a 10% Yield Shortfall, (1980-2005)") + 
  guides(fill=guide_legend(title="Probability")) + 
  scale_fill_gradientn(colours = viridis(8), na.value = "grey80", guide = "legend",
                       breaks = c(0.2, 0.3, 0.4, 0.5, 0.6), limits = c(0.2,0.6))
Plot


Plot <- ggplot(Pi2) + 
  geom_histogram(data = Pi2, aes(meanProb), bins = 100, fill = "seagreen", color = "grey40") + 
  ggtitle("Ave. Probability of a 10% Yield Shortfall, (1980-2005)") + 
  ylab("Number of Counties") + 
  xlab("Mean Probability") +
  xlim(0, 0.8) +
  ylim(0, 150) +
  theme_bw()
Plot

