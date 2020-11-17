

#Subset for the earlier time period
modelmeans <- function(modelnum){
  
  futdat3565 <- readRDS("~/Box/Jeff Work/futureClimate/residuals85.Rdata") %>% filter(year %in% 2035:2065, modelNr == modelnum)
  res.datfut85_3565 = cbind(futdat3565$resid_temp,futdat3565$resid_prec, futdat3565$year)

  
  
  index <- matrix(NA, nrow = nrow(res.datfut85_3565), ncol = 1)
  # fill the index by searching for closest temp/precip match for each observation in our data
  for (i in 1:nrow(res.datfut85_3565)){
    index[i] <- row.names(res.q.grid[which(abs(res.q.grid[,2] - res.datfut85_3565[i,2]) == min(abs(res.q.grid[,2] - res.datfut85_3565[i,2])) & abs(res.q.grid[,1] - res.datfut85_3565[i,1])==min(abs(res.q.grid[,1] - res.datfut85_3565[i,1]))),])
  }
  
  fut <- matrix(NA, nrow = nrow(res.datfut85_3565), ncol=1)
  pred_shortfall <- as.data.frame(res.out)[,6]
  
  for (i in 1:nrow(index)) {
    j <- as.numeric(index[i])
    fut[i] <- pred_shortfall[j]
  }
  
  fut <- as_data_frame(fut)
  fut$modelNr <- modelnum
  return(fut)
}  

fut <- map_dfr(1:21, modelmeans)

futdat3565 <- readRDS("~/Box/Jeff Work/futureClimate/residuals85.Rdata") %>% filter(year %in% 2035:2065)

#Add in fips and year from the future data (I need to do something different with the fips and year)
fut$fips <- futdat3565$fips
fut$year <- futdat3565$year


fut2 <- fut %>% group_by(fips, year) %>% summarise(modelMean = mean(V1))


#Here's the end where I need to take the mean across all of these

#Let's condense this information
fut3 <- fut2 %>% group_by(fips) %>% summarise(meanProbFut = mean(modelMean)) 

#Now I need this information for 1980-2005
#Run this again for 1980-2005----
dat8005 <- read.dta13('~/Box/Jeff Work/Jeff work/ximing/dataSTATA/corn/timeCountyFE/trimester2.dta') %>% filter(year %in% 1980:2005)
res.dat8005 <- cbind(dat8005$dday29CRes,dat8005$precRes, dat8005$year)

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
Pi2 <- Pi %>% group_by(fips) %>% summarise(meanProbHist = mean(V1))

diff <- Pi2 %>% left_join(fut3, by = "fips")
diff <- diff %>% mutate(differ = meanProbFut - meanProbHist)

map_with_data(diff, values = "differ")
df <- us_map(regions = "counties")

Plot <- plot_usmap(regions = "counties", 
                   include = c(.east_north_central, .east_south_central,
                               .mid_atlantic, .south_atlantic,
                               .west_north_central, .west_south_central),
                   data = diff, values = "differ",
                   lines = "grey30", labels = FALSE,
                   label_color = "black") + 
  labs(title ="Difference in Probability of a 10% Yield Shortfall", 
       subtitle = "2035-2065 to 1980-2005 (RCP 8.5)") + 
  guides(fill=guide_legend(title="Probability")) + 
  scale_fill_gradientn(colours = viridis(8), na.value = "grey80", breaks = c(-0.2, 0, 0.2, 0.4, 0.6), limits = c(-0.2,0.6))
Plot  
