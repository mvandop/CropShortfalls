library('gss')
library('copula')
library('matrixcalc')
library('orthopolynom')
source('legendre.r')
source('svd_inv.r')
#source('copula_tri3.r')
source('copula_tri4.r')
source('get_copula2.r')
source('pk3.r')
library('foreign')
library('logspline')

source('lp0.r')
source('lp1.r')
source('lp2.r')
source('lp3.r')


# library('foreign')
# dat=read.dta('weatherXiming.dta')

# install.packages('readstata13')
library('readstata13')

dat=read.dta13('weatherXiming.dta')


get.varlabel(dat)

attach(dat)
summary(dday10C)
summary(dday29C)
summary(dday30C)

summary(dday10C+dday29C+dday30C)

# use annual data
dat.a=dat[dat$trimester==0,]

# use only counties with longitude > -100
dat.a=dat.a[dat.a$longitude>-100,]



fips.freq=by(dat.a$fips,dat.a$fips,length)

# keep counties with at least 10 years of data
fips.freq2=rep(fips.freq,fips.freq)
dat.a2=dat.a[fips.freq2>=10,]

fit0=lm(cornYield~year+I((year-1971)*((year-1971)>0))+I((year-1992)*((year-1992)>0))+avgPlantingDate+avgHarvestDate,data=dat.a2)

# de-mean by county first; it seems R has a hard time
# dealing with de-mean (or county dummies) with large
# sample size

fips.freq=by(dat.a2$fips,dat.a2$fips,length)
dmean=by(dat.a2$cornYield,dat.a2$fips,mean)
dcornYield=dat.a2$cornYield-rep(dmean,fips.freq)

fit=lm(dcornYield~year+I((year-1971)*((year-1971)>0))+I((year-1992)*((year-1992)>0))+avgPlantingDate+avgHarvestDate,data=dat.a2)

y.res=fit$residual

dat.a3=cbind(y.res,dat.a2$tAvg,dat.a2$prec)
fit.cr=copula_tri4(dat.a3,gam=.0001,r=6,m=40,tt=2)
cbind(fit.cr$moment,fit.cr$m_hat)

# estimate the marginal using log-spline
library('logspline')
fit0=logspline(dat.a3[,1],lbound=min(dat.a3[,1])*1.1,ubound=max(dat.a3[,1])*1.1)

# calculate the conditional mean on grids of temperature and prec

q.temp=quantile(dat.a3[,2],prob=seq(.02,.98,.02))
q.prec=quantile(dat.a3[,3],prob=seq(.02,.98,.02))

q2=expand.grid(q.temp,q.prec)
id=expand.grid(seq(.02,.98,.02),seq(.02,.98,.02))


out=matrix(0,nrow=dim(q2)[1],ncol=5)
n=dim(dat.a3)[1]
y0=gauss.quad(1000,c(min(dat.a3[,1])*1.1,max(dat.a3[,1])*1.1))
dy0=dlogspline(y0$pt,fit0)
py0=y0$pt
for (i in 1:length(y0$pt))
  py0[i]=sum(y0$pt[i]>=dat.a3[,1])/(n+1)


py0b=plogspline(y0$pt,fit0)	
summary(dat.a3)
q05=quantile(dat.a3[,1],.05)
# around -36.1

for (i in 1:dim(out)[1])
{
  
  temp=cbind(py0,id[i,1],id[i,2])
  temp2=get_copula2(temp,fit.cr,pp=2)
  out[i,1]=sum(temp2*dy0*y0$wt)
  out[i,2]=sum(y0$pt*temp2*dy0*y0$wt)
  out[i,3]=sum((y0$pt)^2*temp2*dy0*y0$wt)	
  out[i,4]=sum(temp2*dy0*y0$wt*(y0$pt<=q05))
  out[i,5]=sum((y0$pt)^3*temp2*dy0*y0$wt)	
}
summary(out)



z1=out[,2]/out[,1]
z2=out[,3]/out[,1]
z.se=sqrt(z2-z1^2)
z3=out[,4]/out[,1]
z4=out[,5]/out[,1]

rm(dat)
rm(dat.a)
save.image(file='july_4_2015.RData') 

pdf('out.pdf')
# pdf('con_mean.pdf')
persp(q.temp,q.prec,matrix(z1,49,49),xlab='temperature',ylab='precipitation',zlab='',main='conditional mean')
# dev.off()
# pdf('con_mean_contour.pdf')


filled.contour(q.temp,q.prec,matrix(z1,49,49),xlab='temperature',ylab='precipitation',main='conditional mean')
# dev.off()


# pdf(file='con_se.pdf')
persp(q.temp,q.prec,matrix(z.se,49,49),xlab='temperature',ylab='precipitation',main='conditional standard deviation')
# dev.off()
# pdf(file='con_se_contour.pdf')
filled.contour(q.temp,q.prec,matrix(z.se,49,49),xlab='temperature',ylab='precipitation',main='conditional standard deviation')
# dev.off()


# coefficient of variation
# pdf(file='con_cov.pdf')
persp(q.temp,q.prec,matrix(z.se/(z1+mean(dat.a2$cornYield)),49,49),xlab='temperature',ylab='precipitation',main='conditional coefficient of variation')
# dev.off()
# pdf(file='con_cov_contour.pdf')
filled.contour(q.temp,q.prec,matrix(z.se/(z1+mean(dat.a2$cornYield)),49,49),xlab='temperature',ylab='precipitation',main='conditional cv')
# dev.off()

# pdf(file='con_prob.pdf')
persp(q.temp,q.prec,matrix(z3,49,49),xlab='temperature',ylab='precipitation',main='conditional extreme probability')
# dev.off()
# pdf(file='con_prob_contour.pdf')
filled.contour(q.temp,q.prec,matrix(z3,49,49),xlab='temperature',ylab='precipitation',main='conditional extreme probability')
# dev.off()

# conditional skewness
persp(q.temp,q.prec,matrix(-z4,49,49),xlab='temperature',ylab='precipitation',main='conditional (negative) skewness',zlab='')

filled.contour(q.temp,q.prec,matrix(-z4,49,49),xlab='temperature',ylab='precipitation',main='conditional (negative) skewness')
dev.off()

# conditional studentized skewness
persp(q.temp,q.prec,matrix(-z4/z.se^3,49,49),xlab='temperature',ylab='precipitation',main='conditional (negative) skewness',zlab='')

filled.contour(q.temp,q.prec,matrix(-z4/z.se^3,49,49),xlab='temperature',ylab='precipitation',main='conditional (negative) skewness')



# question: how to deal with heteroskedasticity?

# scatterplot
library(scatterplot3d)
qp=expand.grid(q.temp,q.prec)
dat3d=cbind(qp,z1)

temp=q.prec
temp[seq(5,49,5)]=0
qp.temp=expand.grid(q.temp,temp)
dat3d.temp1=dat3d[qp.temp[,2]==0,]


scatterplot3d(dat3d.temp1,type='h',color='blue',tick.marks = F,box=F,xlab = 'temperature',ylab='precipitation',zlab='conditional mean',angle=100,pch=20,scale.y=4)


qp=expand.grid(q.temp,q.prec)
dat3d=cbind(qp,z2)

temp=q.prec
temp[seq(5,49,5)]=0
qp.temp=expand.grid(q.temp,temp)
dat3d.temp2=dat3d[qp.temp[,2]==0,]


scatterplot3d(dat3d.temp2,type='h',color='blue',tick.marks = F,box=F,xlab = 'temperature',ylab='precipitation',zlab='conditional s.e.',angle=100,pch=20,scale.y=4)


# make the grid equally spaced rather than equal-quantile spaced
q.temp=quantile(dat.a3[,2],prob=c(.02,.98))
q.temp=seq(q.temp[1],q.temp[2],length=49)
q.prec=quantile(dat.a3[,3],prob=c(.02,.98))
q.prec=seq(q.prec[1],q.prec[2],length=49)


q2=expand.grid(q.temp,q.prec)

id1=apply(matrix(q.temp,ncol=1),1,function(x) mean(x>dat.a3[,2]))
id2=apply(matrix(q.prec,ncol=1),1,function(x) mean(x>dat.a3[,3]))
id=expand.grid(id1,id2)

out=matrix(0,nrow=dim(q2)[1],ncol=5)
n=dim(dat.a3)[1]
y0=gauss.quad(1000,c(-120,120))
dy0=dlogspline(y0$pt,fit0)
py0=y0$pt
for (i in 1:length(y0$pt))
  py0[i]=sum(y0$pt[i]>=dat.a3[,1])/(n+1)


py0b=plogspline(y0$pt,fit0)	
summary(dat.a3)
q05=quantile(dat.a3[,1],.05)
# around -36.6

for (i in 1:dim(out)[1])
{
  
  temp=cbind(py0,id[i,1],id[i,2])
  temp2=get_copula2(temp,fit.cr,pp=2)
  out[i,1]=sum(temp2*dy0*y0$wt)
  out[i,2]=sum(y0$pt*temp2*dy0*y0$wt)
  out[i,3]=sum((y0$pt)^2*temp2*dy0*y0$wt)	
  out[i,4]=sum(temp2*dy0*y0$wt*(y0$pt<=q05))
  out[i,5]=sum((y0$pt)^3*temp2*dy0*y0$wt)	
}
summary(out)



z1=out[,2]/out[,1]
z2=out[,3]/out[,1]
z.se=sqrt(z2-z1^2)
z3=out[,4]/out[,1]
z4=out[,5]/out[,1]

filled.contour(q.temp,q.prec,matrix(z1,49,49),xlab='temperature',ylab='precipitation',main='conditional mean')
# dev.off()


# pdf(file='con_se.pdf')
# persp(q.temp,q.prec,matrix(z.se,19,19),xlab='temperature',ylab='precipitation',zlab='conditional standard errors')
# dev.off()
# pdf(file='con_se_contour.pdf')
filled.contour(q.temp,q.prec,matrix(z.se,49,49),xlab='temperature',ylab='precipitation',main='conditional standard deviation')
# dev.off()


# coefficient of variation
# pdf(file='con_cov.pdf')
# persp(q.temp,q.prec,matrix(z.se/(z1+mean(dat.a2$cornYield)),19,19),xlab='temperature',ylab='precipitation',zlab='conditional coefficient of variation')
# dev.off()
# pdf(file='con_cov_contour.pdf')
filled.contour(q.temp,q.prec,matrix(z.se/(z1+mean(dat.a2$cornYield)),49,49),xlab='temperature',ylab='precipitation',main='conditional cv')
# dev.off()

# pdf(file='con_prob.pdf')
# persp(q.temp,q.prec,matrix(z3,19,19),xlab='temperature',ylab='precipitation',zlab='conditional probability')
# dev.off()
# pdf(file='con_prob_contour.pdf')
filled.contour(q.temp,q.prec,matrix(z3,49,49),xlab='temperature',ylab='precipitation',main='conditional extreme probability')
# dev.off()

# conditional skewness
filled.contour(q.temp,q.prec,matrix(-z4,49,49),xlab='temperature',ylab='precipitation',main='conditional skewness')


filled.contour(q.temp,q.prec,matrix(-z4/z.se^3,49,49),xlab='temperature',ylab='precipitation',main='conditional skewness')


qp=expand.grid(q.temp,q.prec)
dat3d=cbind(qp,z1)

temp=q.prec
temp[seq(5,49,5)]=0
qp.temp=expand.grid(q.temp,temp)
dat3d.temp1=dat3d[qp.temp[,2]==0,]


pdf('mean_scatterplot.pdf')
scatterplot3d(dat3d.temp1,type='h',color='blue',tick.marks = F,box=F,xlab = 'temperature',ylab='precipitation',zlab='conditional mean',angle=100,pch=20,scale.y=4)
dev.off()


qp=expand.grid(q.temp,q.prec)
dat3d=cbind(qp,z2)

temp=q.prec
temp[seq(5,49,5)]=0
qp.temp=expand.grid(q.temp,temp)
dat3d.temp2=dat3d[qp.temp[,2]==0,]


pdf('se_scatterplot.pdf')
scatterplot3d(dat3d.temp2,type='h',color='blue',tick.marks = F,box=F,xlab = 'temperature',ylab='precipitation',zlab='conditional s.e.',angle=100,pch=20,scale.y=4)
dev.off()



qp=expand.grid(q.temp,q.prec)
dat3d=cbind(qp,-z4)

temp=q.prec
temp[seq(5,49,5)]=0
qp.temp=expand.grid(q.temp,temp)
dat3d.temp2=dat3d[qp.temp[,2]==0,]

pdf('skew_scatterplot.pdf')
scatterplot3d(dat3d.temp2,type='h',color='blue',tick.marks = F,box=F,xlab = 'temperature',ylab='precipitation',main='conditional (negative) skewness',angle=100,pch=20,scale.y=4)
dev.off()



qp=expand.grid(q.temp,q.prec)
dat3d=cbind(qp,-z4/z.se^3)

temp=q.prec
temp[seq(5,49,5)]=0
qp.temp=expand.grid(q.temp,temp)
dat3d.temp2=dat3d[qp.temp[,2]==0,]


scatterplot3d(dat3d.temp2,type='h',color='blue',tick.marks = F,box=F,xlab = 'temperature',ylab='precipitation',zlab='',main='conditional (negative) skewness',angle=100,pch=20,scale.y=4)



qp=expand.grid(q.temp,q.prec)
dat3d=cbind(qp,z3)

temp=q.prec
temp[seq(5,49,5)]=0
qp.temp=expand.grid(q.temp,temp)
dat3d.temp2=dat3d[qp.temp[,2]==0,]


pdf('prob_scatterplot.pdf')
scatterplot3d(dat3d.temp2,type='h',color='blue',tick.marks = F,box=F,xlab = 'temperature',ylab='precipitation',zlab='',main='conditional extreme low yield',angle=100,pch=20,scale.y=4)
dev.off()




qp=expand.grid(q.temp,q.prec)
dat3d=cbind(qp,z.se/(z1+mean(dat.a2$cornYield)))

temp=q.prec
temp[seq(5,49,5)]=0
qp.temp=expand.grid(q.temp,temp)
dat3d.temp2=dat3d[qp.temp[,2]==0,]


scatterplot3d(dat3d.temp2,type='h',color='blue',tick.marks = F,box=F,xlab = 'temperature',ylab='precipitation',main='conditional coefficient of variation',angle=100,pch=20,scale.y=4)


# put in more weather variables in the regressions
fit2=lm(dcornYield~year+I((year-1971)*((year-1971)>0))+I((year-1992)*((year-1992)>0))+avgPlantingDate+avgHarvestDate+prec2+I(tMax-tMin),data=dat.a2)

y2.res=fit2$residual

dat.a4=cbind(y2.res,dat.a2$tAvg,dat.a2$prec)
fit.cr2=copula_tri4(dat.a4,gam=.0001,r=6,m=40,tt=2)

cbind(fit.cr2$moment,fit.cr2$m_hat)

# estimate the marginal using log-spline
library('logspline')
summary(dat.a4[,1])
fit0=logspline(dat.a4[,1],lbound=min(dat.a4[,1])*1.1,ubound=max(dat.a4[,1])*1.1)

# calculate the conditional mean on grids of temperature and prec

q.temp=quantile(dat.a4[,2],prob=seq(.02,.98,.02))
q.prec=quantile(dat.a4[,3],prob=seq(.02,.98,.02))

q2=expand.grid(q.temp,q.prec)
id=expand.grid(seq(.02,.98,.02),seq(.02,.98,.02))


out=matrix(0,nrow=dim(q2)[1],ncol=5)
n=dim(dat.a4)[1]
y0=gauss.quad(1000,c(min(dat.a4[,1])*1.1,max(dat.a4[,1])*1.1))
dy0=dlogspline(y0$pt,fit0)
py0=y0$pt
for (i in 1:length(y0$pt))
  py0[i]=sum(y0$pt[i]>=dat.a4[,1])/(n+1)


py0b=plogspline(y0$pt,fit0)	
summary(dat.a4)
q05=quantile(dat.a4[,1],.05)
# around -36.1

for (i in 1:dim(out)[1])
{
  
  temp=cbind(py0,id[i,1],id[i,2])
  temp2=get_copula2(temp,fit.cr2,pp=2)
  out[i,1]=sum(temp2*dy0*y0$wt)
  out[i,2]=sum(y0$pt*temp2*dy0*y0$wt)
  out[i,3]=sum((y0$pt)^2*temp2*dy0*y0$wt)	
  out[i,4]=sum(temp2*dy0*y0$wt*(y0$pt<=q05))
  out[i,5]=sum((y0$pt)^3*temp2*dy0*y0$wt)	
}
summary(out)



z1=out[,2]/out[,1]
z2=out[,3]/out[,1]
z.se=sqrt(z2-z1^2)
z3=out[,4]/out[,1]
z4=out[,5]/out[,1]



# pdf('con_mean.pdf')
#persp(q.temp,q.prec,matrix(z1,19,19),xlab='temperature',ylab='precipitation',zlab='conditional mean')
# dev.off()
# pdf('con_mean_contour.pdf')
pdf('out.pdf')
filled.contour(q.temp,q.prec,matrix(z1,49,49),xlab='temperature',ylab='precipitation',main='conditional mean')
# dev.off()


# pdf(file='con_se.pdf')
# persp(q.temp,q.prec,matrix(z.se,19,19),xlab='temperature',ylab='precipitation',zlab='conditional standard errors')
# dev.off()
# pdf(file='con_se_contour.pdf')
filled.contour(q.temp,q.prec,matrix(z.se,49,49),xlab='temperature',ylab='precipitation',main='conditional standard deviation')
# dev.off()


# coefficient of variation
# pdf(file='con_cov.pdf')
# persp(q.temp,q.prec,matrix(z.se/(z1+mean(dat.a2$cornYield)),19,19),xlab='temperature',ylab='precipitation',zlab='conditional coefficient of variation')
# dev.off()
# pdf(file='con_cov_contour.pdf')
filled.contour(q.temp,q.prec,matrix(z.se/(z1+mean(dat.a2$cornYield)),49,49),xlab='temperature',ylab='precipitation',main='conditional cv')
# dev.off()

# pdf(file='con_prob.pdf')
# persp(q.temp,q.prec,matrix(z3,19,19),xlab='temperature',ylab='precipitation',zlab='conditional probability')
# dev.off()
# pdf(file='con_prob_contour.pdf')
filled.contour(q.temp,q.prec,matrix(z3,49,49),xlab='temperature',ylab='precipitation',main='conditional extreme probability')
# dev.off()

# conditional skewness
filled.contour(q.temp,q.prec,matrix(-z4,49,49),xlab='temperature',ylab='precipitation',main='conditional negative skewness')
dev.off()


fips=unique(dat.a2$fips)
state=floor(fips/1000)
length(unique(state))

state_code=read.csv('state_fips_code.csv')
state_code=state_code[1:53,]
state_code

state.list=unique(state)
state_code_list=state_code[state_code[,3]%in%state.list,]
state_code_list

# one can use rought 36.5 latitude for north and south
summary(dat.a2$tAvg[dat.a2$latitude>36.5])
summary(dat.a2$tAvg[dat.a2$latitude<36.5])




# evalute the conditional density at 2014 average temperature and precipitation
summary(dat.a2$tAvg[dat.a2$year==2014])
summary(dat.a2$prec[dat.a2$year==2014])

summary(dat.a2$tAvg[dat.a2$year==2014&dat.a2$latitude>36.5])
summary(dat.a2$tAvg[dat.a2$year==2014&dat.a2$latitude<36.5])

summary(dat.a2$prec[dat.a2$year==2014&dat.a2$latitude>36.5])
summary(dat.a2$prec[dat.a2$year==2014&dat.a2$latitude<36.5])


library('logspline')
fit0=logspline(dat.a4[,1],lbound=min(dat.a4[,1])*1.1,ubound=max(dat.a4[,1])*1.1)

# calculate the conditional mean on grids of temperature and prec

q.temp=c(mean(dat.a2$tAvg[dat.a2$year==2014&dat.a2$latitude>36.5]),mean(dat.a2$tAvg[dat.a2$year==2014&dat.a2$latitude<36.5]))

q.temp=c(q.temp[1],q.temp[1]+c(1,2,3,5),q.temp[2],q.temp[2]+c(1,2,3,5))

q.prec=c(mean(dat.a2$prec[dat.a2$year==2014&dat.a2$latitude>36.5]),mean(dat.a2$prec[dat.a2$year==2014&dat.a2$latitude<36.5]))
q.prec=rep(q.prec,each=5)


id1=apply(matrix(q.temp,ncol=1),1,function(x) mean(x>dat.a4[,2]))
id2=apply(matrix(q.prec,ncol=1),1,function(x) mean(x>dat.a4[,3]))
id=cbind(id1,id2)


#q2=expand.grid(q.temp,q.prec)
#id=expand.grid(seq(.02,.98,.02),seq(.02,.98,.02))
q2=cbind(q.temp,q.prec)

out=matrix(0,nrow=dim(q2)[1],ncol=5)
n=dim(dat.a4)[1]
y0=gauss.quad(1000,c(min(dat.a4[,1])*1.1,max(dat.a4[,1])*1.1))
dy0=dlogspline(y0$pt,fit0)
py0=y0$pt
for (i in 1:length(y0$pt))
  py0[i]=sum(y0$pt[i]>=dat.a4[,1])/(n+1)


py0b=plogspline(y0$pt,fit0)	
summary(dat.a4)
q05=quantile(dat.a4[,1],.05)

# around -36.6

for (i in 1:dim(out)[1])
{
  
  temp=cbind(py0,id[i,1],id[i,2])
  temp2=get_copula2(temp,fit.cr2,pp=2)
  out[i,1]=sum(temp2*dy0*y0$wt)
  out[i,2]=sum(y0$pt*temp2*dy0*y0$wt)
  out[i,3]=sum((y0$pt)^2*temp2*dy0*y0$wt)	
  out[i,4]=sum(temp2*dy0*y0$wt*(y0$pt<=q05))
  out[i,5]=sum((y0$pt)^3*temp2*dy0*y0$wt)	
}
summary(out)



z1=out[,2]/out[,1]
z2=out[,3]/out[,1]
z.se=sqrt(z2-z1^2)
z3=out[,4]/out[,1]
z4=out[,5]/out[,1]


z1=matrix(z1,ncol=2)
z1[,1]=z1[,1]+mean(dat.a2$cornYield[dat.a2$year==2014&dat.a2$latitude>36.5])
z1[,2]=z1[,2]+mean(dat.a2$cornYield[dat.a2$year==2014&dat.a2$latitude<36.5])

colnames(z1)=c('north','south')
rownames(z1)=c('current','+1','+2','+3','+5')

barplot(z1,beside=T,main='Predicted mean under 0 and 1,2,3,5 degree increase in temperature')
z1.b=z1[-1,]/t(matrix(t(z1[1,]),ncol=4,nrow=2))-1
barplot(z1.b,beside=T,legend=rownames(z1.b),main='Predicted percentage change in mean yield')


z.se=matrix(z.se,ncol=2)
barplot(z.se,beside=T,main='Predicted s.e. under 0 and 1,2,3,5 degree increase in temperature')
zse.b=z.se[-1,]/t(matrix(t(z.se[1,]),ncol=4,nrow=2))-1
colnames(zse.b)=c('north','south')
barplot(zse.b,beside=T,legend=rownames(z1.b),main='Predicted percentage change in s.e.')


z3=matrix(z3,ncol=2)
barplot(z.se,beside=T,main='Predicted s.e. under 0 and 1,2,3,5 degree increase in temperature')
z3.b=z3[-1,]/t(matrix(t(z3[1,]),ncol=4,nrow=2))-1
colnames(z3.b)=c('north','south')
barplot(z3.b,beside=T,legend=rownames(z1.b),main='Predicted percentage change in extreme probability')

z.cv=z.se/z1
barplot(z.cv,beside=T)
zcv.b=z.cv[-1,]/t(matrix(t(z.cv[1,]),ncol=4,nrow=2))-1
colnames(zcv.b)=c('north','south')
barplot(zcv.b,beside=T,legend=rownames(z1.b),main='Predicted percentage change in c.v.')


z4=matrix(z4,ncol=2)
barplot(z.se,beside=T,main='Predicted s.e. under 0 and 1,2,3,5 degree increase in temperature')
z4.b=z4[-1,]/t(matrix(t(z4[1,]),ncol=4,nrow=2))-1
colnames(z4.b)=c('north','south')
barplot(z4.b,beside=T,legend=rownames(z1.b),main='Predicted percentage change in extreme probability')

pdf(file='out2.pdf')
barplot(z1.b,beside=T,legend=rownames(z1.b),main='Predicted percentage change in mean yield')

barplot(zse.b,beside=T,legend=rownames(z1.b),main='Predicted percentage change in s.e.')


barplot(zcv.b,beside=T,legend=rownames(z1.b),main='Predicted percentage change in c.v.')

barplot(z3.b,beside=T,legend=rownames(z1.b),main='Predicted percentage change in extreme probability')

barplot(z4.b,beside=T,legend=rownames(z1.b),main='Predicted percentage change in negative skewness')

dev.off()

