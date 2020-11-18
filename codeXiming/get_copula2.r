get_copula2=function(dat,fit,pp)
{
# enter data as rank already
dat=matrix(dat,ncol=3)

u1=dat[,1]
u2=dat[,2]
u3=dat[,3]

r=fit$r

uu=NULL
for (i1 in 0:r)
{
	for (i2 in 0:(r-i1))
	{
		for (i3 in 0:(r-i1-i2))
		{
		 uu=cbind(uu,legendre(u1,i1,pp)*legendre(u2,i2,pp)*legendre(u3,i3,pp))
		}
	}
}

lam=fit$lam
den=exp(uu%*%lam)
return(den=den)
}