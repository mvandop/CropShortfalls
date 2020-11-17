
# Set working directory
#setwd("~/Box/Jeff Work/Jeff work/")

# load libraries that are needed
library('gss')
library('copula')
library('matrixcalc')
library('orthopolynom')
library('foreign')
library('logspline')
library('readstata13')
library('plm')
library('doBy')
library('tidyverse')
library('dplyr')
library('RColorBrewer')
library('scales')
library('viridis')

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

# loop over crop
for (c in 1:length(cropList)) {
  if (cropList[c] == 1) {
    # corn yields
    cropName = 'Corn'
    ddayUB = 29;
  } else if (cropList[c] == 2) {
    # soybeans yields
    cropName = 'Soybeans'
    ddayUB = 30;
  } else {
    print(paste0('Crop', cropList[c], ' is not defined!!!'))
    error
  }
  
  # loop over trimesters
  for (t in 1:length(trimesterList)) {
    
    # loop over models
    for (m in 1:length(modelList)) {
      #-------------------------------------------------------------------------------------------------------------------------------;
      # pick apprpriate model;
      #-------------------------------------------------------------------------------------------------------------------------------;
      if (modelList[m] == 0) {
        modelName = 'time'
      } else if (modelList[m] == 1) {
        modelName = 'soil'
      } else if (modelList[m] == 2) {
        modelName = 'weatherSoil'
      } else if (modelList[m] == 3) {
        modelName = 'time'
      } else if (modelList[m] == 4) {
        modelName = 'weatherSoilTime'
      } else if (modelList[m] == 11) {
        modelName = 'timeCountyFE'
      } else if (modelList[m] == 12) {
        modelName = 'weatherTimeCountyFE'
      } else {
        print(paste0('Model ', modelList[m], ' is not defined!!!'))
        error
      }
      
      # loop over subsets
      for (s in 1:length(subsetList)) {
        #-------------------------------------------------------------------------------------------------------------------------------;
        # pick apprpriate subset;
        #-------------------------------------------------------------------------------------------------------------------------------;
        if (subsetList[s] == 0) {
          # no subset
          subsetName = ''
        } else if (subsetList[s] == 1) {
          # hotter third
          subsetName = '_hotCounties'
        } else if (subsetList[s] == 2) {
          # colder third
          subsetName = '_coldCounties'
        } else if (subsetList[s] == 11) {
          # hotter third
          subsetName = '_year1950_1972'
        } else if (subsetList[s] == 12) {
          # colder third
          subsetName = '_year1973_1994'
        } else if (subsetList[s] == 13) {
          # colder third
          subsetName = '_year1995_2016'
        } else {
          print(paste0('Subset', subsetList[s], ' is not defined!!!'))
          error
        }
        
        #-------------------------------------------------------------------------------------------------------------------------------;
        # load data set
        #-------------------------------------------------------------------------------------------------------------------------------;
        dat=read.dta13(paste0('~/Box/Jeff Work/Jeff work/ximing/dataSTATA/',tolower(cropName),'/',modelName,'/trimester',trimesterList[t],subsetName,'.dta'))
        if (modelList[m] == 0) {
          fileName = 'noControl'
          if (cropList[c] == 1) {
            # corn yields
            res.dat = cbind(dat$logYield,dat$dday29C,dat$prec)
          } else if (cropList[c] == 2) {
            # soybeans yields
            res.dat = cbind(dat$logYield,dat$dday30C,dat$prec)
          }
        } else {
          fileName = modelName
          if (cropList[c] == 1) {
            # corn yields
            res.dat = cbind(dat$logYieldRes,dat$dday29CRes,dat$precRes)
          } else if (cropList[c] == 2) {
            # soybeans yields
            res.dat = cbind(dat$logYieldRes,dat$dday30CRes,dat$precRes)
          }
        }
        #get.varlabel(dat)
        res.dat2 = cbind(dat$logYield,dat$logYieldMean)
        
        #-------------------------------------------------------------------------------------------------------------------------------;
        # estimate model;
        #-------------------------------------------------------------------------------------------------------------------------------;
        res.fit = copula_tri4(res.dat,gam = .0001,r = 6,m = 40,tt = 2)
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
        res.y0 = gauss.quad(1000,c(min(res.dat[,1]) - 0.1*(max(res.dat[,1])-min(res.dat[,1])), max(res.dat[,1])+0.1*(max(res.dat[,1])-min(res.dat[,1])) ))

        res.dy0 = dlogspline(res.y0$pt,res.fitMarginal)
        res.py0 = res.y0$pt 
        for (i in 1:length(res.y0$pt)) {
          res.py0[i]=sum(res.y0$pt[i]>=res.dat[,1])/(res.nObs+1)
        }
        
        res.py0b = plogspline(res.y0$pt,res.fitMarginal)
        res.q.yield5 = quantile(res.dat[,1],.05)
        res.q.yield90 = quantile(res.dat[,1], .90)
        res.q.M5 = mean(res.dat[,1])-0.05
        res.q.M10 = mean(res.dat[,1])-0.1
        res.q.M15 = mean(res.dat[,1])-0.15
        res.q.M33 = mean(res.dat[,1])-1/3
        res.q.M90 = mean(res.dat[,1])-.9
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
          
          
          # numerical integration of the conditional density (not normalized)
          res.out[i,1] = sum(temp2 * res.dy0 * res.y0$wt)
          
          # conditional mean (not normalized)
          res.out[i,2] = sum(res.y0$pt * temp2 * res.dy0 * res.y0$wt)
          
          # conditional variance (not normalized)
          res.out[i,3] = sum((res.y0$pt-(mean(res.dat[,1])-0.1))^2 * temp2 * res.dy0 * res.y0$wt)
          
          # downside variance (un-normalized) for values of res.y0$pt < =z, where z can be an arbitrary cutoff point, for instance z=m1 
          res.out[i,4] = sum((res.y0$pt-(mean(res.dat[,1])-0.1))^2 * temp2 * res.dy0 * res.y0$wt)
          
          # for models with county fixed effects, log yields give directly yield shortfall
          if (modelList[m] > 10) {
            res.out[i,5] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= -0.05))
            res.out[i,6] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= -0.1))
            res.out[i,7] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= -0.15))
            res.out[i,8] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= -0.9))
            res.out[i,9] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= res.q.se))
          } else {
            res.out[i,5] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= res.q.M5))
            res.out[i,6] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= res.q.M10))
            res.out[i,7] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= res.q.M15))
            res.out[i,8] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= res.q.M90))
            res.out[i,9] = sum(temp2 * res.dy0 * res.y0$wt * (res.y0$pt <= res.q.se))
          } 
          
        }
        rm(i)
        
        res.z1  = res.out[,2] / res.out[,1]
        res.z2  = res.out[,3] / res.out[,1]
        res.zse = sqrt(res.z2 - res.z1^2)
        res.z3  = res.out[,4] / res.out[,1]
        res.z4  = res.out[,5] / res.out[,1]
        res.z5  = res.out[,6] / res.out[,1]
        res.z6  = res.out[,7] / res.out[,1]
        res.z7  = res.out[,8] / res.out[,1]
        res.z8 = res.out[,9] / res.out[,1]
        
        index <- matrix(NA, nrow = nrow(res.dat), ncol = 1)
        
        # fill the index by searching for closest temp/precip match for each observation in our data
        for (i in 1:nrow(res.dat)){
          index[i] <- row.names(res.q.grid[which(abs(res.q.grid[,2] - res.dat[i,3]) == min(abs(res.q.grid[,2] - res.dat[i,3])) & abs(res.q.grid[,1] - res.dat[i,2])==min(abs(res.q.grid[,1] - res.dat[i,2]))),])
        }
        Pi <- matrix(NA, nrow = nrow(res.dat), ncol=1)
        pred_prob <- as.data.frame(res.out)[,9]
        for (i in 1:nrow(index)) {
          j <- as.numeric(index[i])
          Pi[i] <- pred_prob[j]
        }
        
        # Di: indicator for whether the yield falls below 1 SD of historical mean
        Di <- as.matrix(ifelse(res.dat[,1] - mean(res.dat[,1]) <= -sd(res.dat[,1]),1,0))
        
        # Measure of Fit
        mof <- paste0("Measure of Fit ", sum((Di - Pi)^2))
        
        #-------------------------------------------------------------------------------------------------------------------------------;
        # save results;
        #-------------------------------------------------------------------------------------------------------------------------------;
        # Save Measure of Fit
        #capture.output(mof, file = paste0('~/Box/Jeff Work/Jeff work/ximing/dataR/log',cropName,'Yield_',fileName,'_trimester',trimesterList[t],subsetName,'mof_nonlinear90.txt'))
        
        # save data file;
        #save.image(file=paste0('~/Box/Jeff Work/Jeff work/ximing/dataR/log',cropName,'Yield_',fileName,'_trimester',trimesterList[t],subsetName,'90.RData')) 
        
        #flip color palette
        colpal <- viridis(30)
        colpal <- colpal[30:1]
     
        
        # mean;
        postscript(paste0('~/Box/Jeff Work/Jeff work/ximing/dataR/documents/figures/mean/log', cropName, 'Yield_',modelName,'_trimester',trimesterList[t],subsetName,'_contour.ps'))
        filled.contour(res.q.temp,res.q.prec,matrix(res.z1,49,49), nlevels = 20,xlab=paste0('Degree Days ',ddayUB,'C'),ylab = 'Precipitation (m)', main=paste0('Conditional Mean, ', cropName), color.palette = viridis, zlim = c(-0.5, 0.4))
        dev.off()
       

        # standard deviation
         postscript(paste0('~/Box/Jeff Work/Jeff work/ximing/dataR/documents/figures/sd/log', cropName, 'Yield_',modelName,'_trimester',trimesterList[t],subsetName,'_contour.ps'))
         filled.contour(res.q.temp,res.q.prec,matrix(res.zse,49,49),xlab=paste0('Degree Days ',ddayUB,'C'),ylab='Precipitation (m)',main=paste0('Conditional Standard Deviation, ', cropName), color.palette = viridis, zlim = c(0.08, 0.28))
         dev.off()
     
       # coefficient of variation
         postscript(paste0('~/Box/Jeff Work/Jeff work/ximing/dataR/documents/figures/sd/log', cropName, 'Yield_',modelName,'_trimester',trimesterList[t],subsetName,'_coeffVar_contour.ps'))
         filled.contour(res.q.temp,res.q.prec,matrix(res.zse/mean(res.dat2[,1]),49,49),xlab=paste0('Degree Days ',ddayUB,'C'),ylab='Precipitation (m)',main=paste0('Conditional Coefficient of Variation, ', cropName), color.palette = viridis, zlim = c(0.02, 0.08))
         dev.off()
        
       # skewness
         postscript(paste0('~/Box/Jeff Work/Jeff work/ximing/dataR/documents/figures/skewness/log', cropName, 'Yield_', modelName,'_trimester',trimesterList[t],subsetName,'_contourNoLim.ps'))
         filled.contour(res.q.temp,res.q.prec,matrix(-res.z3,49,49),xlab=paste0('Degree Days ',ddayUB,'C'),ylab='Precipitation (m)',main=paste0('Conditional Skewness, ', cropName), color.palette = viridis, zlim = c(-0.03, -0.09))
         dev.off()
 
       # conditional probability of 10% yield shortfall
         postscript(paste0('~/Box/Jeff Work/Jeff work/ximing/dataR/documents/figures/yieldShortFall/log',cropName,'YieldShortfall10_',fileName,'_trimester',trimesterList[t],subsetName,'_contour.ps'))
         filled.contour(res.q.temp,res.q.prec,matrix(res.z5,49,49),xlab=paste0('Degree Days ',ddayUB,'C'),ylab='Precipitation (m)',main=paste0('Conditional Probability of 10% Yield Shortfall, ', cropName), col = colpal, zlim = c(0,0.7))
         dev.off()


        # clean up data for current subset;
        rm(res.q.grid,res.q.gridID,res.q.temp,res.q.prec,res.q.yield5,res.q.M10,subsetName,fileName,res.dat,res.dat2);
        rm(res.out,res.dy0,res.fit,res.fitMarginal,res.nObs,res.py0,res.y0,res.z1,res.z2,res.z3,res.z4,res.z5,res.z6,res.z7,res.z8, res.zse);
      }
      rm(modelName)
    }
  }
  rm(cropName)
}
