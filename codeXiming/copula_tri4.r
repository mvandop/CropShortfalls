
copula_tri4=function(dat,r,gam,m,tol=1e-10,maxit=100,lam=0,tt=1,aa=0)
{
# tri-variate maxent density, with penalty on 2nd derivatives
# r: order of total moments
# m: level of smolyak quadrature
# legendre polynomials
# if tt=2, use diagonal penalty matrix
# aa: ridge for the penalty matrix to make it positive semidefinite

# m=30
# dat=matrix(runif(300),nrow=100,ncol=3)
#xrange=c(0,1)
#yrange=c(0,1)
# r=4
# gam=.1
# tol=1e-10
# maxit=100
# library('gss')
# library('copula')
# library('orthopolynom')
# pp=2
# source('legendre.r')
# source('svd_inv.r')

# n.cop=normalCopula(c(.1,.1,.1),dim=3,dispstr='un')
# fit=fitCopula(n.cop,dat)
# n.cop=normalCopula(fit@estimate,dim=3,dispstr='un')
# u0=dCopula(dat,n.cop)
n=dim(dat)[1]
pp=2 # legendre polynomials
# use peseudo data
u1=rank(dat[,1])/(n+1)
u2=rank(dat[,2])/(n+1)
u3=rank(dat[,3])/(n+1)


x0=smolyak.quad(3,m)
x1=x0$pt[,1]
x2=x0$pt[,2]
x3=x0$pt[,3]
w=x0$wt
# f0=dCopula(x0$pt,n.cop)

# f0=1
# u0=1

# xx1=NULL
# xx2=NULL
# xx3=NULL
# uu1=NULL
# uu2=NULL
# uu3=NULL



# for (i in 1:r)
# {
	# xx1=cbind(xx1,legendre(x1,i,pp))
	# xx2=cbind(xx2,legendre(x2,i,pp))
	# xx3=cbind(xx3,legendre(x3,i,pp))
	# uu1=cbind(uu1,legendre(u1,i,pp))
	# uu2=cbind(uu2,legendre(u2,i,pp))
	# uu3=cbind(uu3,legendre(u3,i,pp))
# }

xx=NULL
uu=NULL
id=NULL
for (i1 in 0:r)
{
	for (i2 in 0:(r-i1))
	{
		for (i3 in 0:(r-i1-i2))
		{
		 id=rbind(id,c(i1,i2,i3))
		 xx=cbind(xx,legendre(x1,i1,pp)*legendre(x2,i2,pp)*legendre(x3,i3,pp))
		 uu=cbind(uu,legendre(u1,i1,pp)*legendre(u2,i2,pp)*legendre(u3,i3,pp))
		}
	}
}

id=id[-1,]
xx=xx[,-1]
uu=uu[,-1]

mm=matrix(apply(uu,2,mean),ncol=1)

nr=length(mm)
lam=matrix(0,nrow=nr,ncol=1)

H=matrix(0,nr,nr)

change=1
it=0
H=matrix(0,nr,nr)
Change=NULL

# id2=apply(id,1,sum)
# id3=id[id2==3,]

# fac2=function(k,r)
# {
	# out=k
	# for (i in 1:length(k))
	# {
	# if (k[i]>=r)
		# out[i]=(factorial(k[i]+r))/(factorial(k[i]-r))
	# else
		# out[i]=0
	# }
	# return(out)	
# }		




# pk=matrix(0,nrow=dim(id)[1],ncol=1)
# for (i in 1:dim(id)[1])
# {
	# if (sum(id[i,])>2)
	# {
		# temp=0
		# for (j in 1:dim(id3)[1])
		# {
			# temp=temp+fac2(id[i,1],id3[j,1])*fac2(id[i,2],id3[j,2])*fac2(id[i,3],id3[j,3])
		# }
	# pk[i]=temp	
	# }

# }

# K=diag(nr)
# diag(K)=pk

K=pk3(r,1000)

K=K+diag(nr)*aa

if (tt==2)
{
	temp=diag(K)
	K=diag(nr)
	diag(K)=temp
}

while (change>tol&it<maxit)
{
	it=it+1
	d=exp(xx%*%lam)
	temp=sum(d*w)
	d=d/temp
	dw=d*w
	m_hat=matrix(0,nrow=r,ncol=1)
	for (i in 1:nr)
	{
		m_hat[i]=sum(xx[,i]*dw)
	}

	for (i in 1:nr)
	{
		temp=xx[,i]*dw
		for (j in (1:i))
		{
			t=sum(xx[,j]*temp)-m_hat[i]*m_hat[j]
			H[i,j]=t
			H[j,i]=t
		}
	}
	

	b=mm-m_hat-gam*K%*%lam
	H_inv=svd_inv(H+gam*K)
	incr=-H_inv%*%b
	lam=lam-incr
	change=-t(b)%*%incr/nr
	# change=sqrt(abs(change))
	Change=c(Change,change)
	
}

d=exp(xx%*%lam)
temp=sum(d*w)
d=d/temp
	
ent=sum(lam*mm)-log(temp)

if (it==maxit)
	ent=NULL

bic=-2*length(u1)*ent+nr*log(length(u1))
aic=-2*length(u1)*ent+2*nr

lam=matrix(c(-log(temp),lam),ncol=1)
den=exp(uu%*%lam[-1]+lam[1])

# H=f$hessian+gam*f$K
	
	ones=matrix(1,ncol=n)
	# temp=-matrix.trace((uu)%*%H_inv%*%t(uu))/(n-1)+((ones)%*%uu)%*%H_inv%*%t(ones%*%uu)/n/(n-1)
	# temp=eigen((uu)%*%H_inv%*%t(uu))
	# temp=-Re(sum(temp))/(n-1)+((ones)%*%uu)%*%H_inv%*%t(ones%*%uu)/n/(n-1)
	# temp=((ones)%*%uu)%*%H_inv%*%t(ones%*%uu)/n/(n-1)
	
	# uu_hat=apply(uu,2,mean)
	uu_diff=apply(uu,2,function(x) x-mean(x))/(n-1)
	temp=0
	for (i in 1:n)
		temp=temp+uu[i,]%*%H_inv%*%(uu_diff[i,])
	
	
	cv=sum(log(den))-temp

list(lam=lam,hessian=H,moment=mm,ent=ent,bic=bic,aic=aic,it=it,den=den,change=change,m_hat=m_hat,Change=Change,r=r,d_ev=d,K=K,gam=gam,xx=xx,uu=uu,H_inv=H_inv,cv=cv)
}
