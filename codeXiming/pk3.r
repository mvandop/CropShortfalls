pk3=function(r,n)
{
# calculate the penalty matrix
# r: order of legendre polynomials
# n: number of gauss quadrature points
xx=gauss.quad(n,c(0,1))
id=NULL
for (i1 in 0:r)
{
	for (i2 in 0:(r-i1))
	{
		for (i3 in 0:(r-i1-i2))
		{
		 id=rbind(id,c(i1,i2,i3))
		}
	}
}

id=id[-1,]

nr=dim(id)[1]

id2=apply(id,1,sum)
id3=id[id2==3,]


df2=function(x,m,r)
{
	if (r==0)
		return(lp0(x,m))
	if (r==1)
		return(lp1(x,m))
	if (r==2)
		return(lp2(x,m))
	if (r==3)
		return(lp3(x,m))
}

d1=d2=d3=list()
for (i in 1:nr)
{
	temp1=temp2=temp3=NULL
	
	a=id[i,]
b=NULL
for (i1 in 0:a[1])
{
	for (i2 in 0:a[2])
	{
		for (i3 in 0:a[3])
			b=rbind(b,c(i1,i2,i3))
	}
}

c=apply(b,1,sum)
d=b[c==3,]
d=matrix(d,ncol=3)
cc=sum(c==3)
if (cc==0)
	temp1=temp2=temp3=0
if (cc>0)
{
for (j in 1:cc)
{
		temp1=cbind(temp1,df2(xx$pt,id[i,1],d[j,1]))
		temp2=cbind(temp2,df2(xx$pt,id[i,2],d[j,2]))
		temp3=cbind(temp3,df2(xx$pt,id[i,3],d[j,3]))

}
}
d1[[i]]=temp1
d2[[i]]=temp2
d3[[i]]=temp3	
}


pk=matrix(0,nrow=nr,ncol=nr)
for (i1 in 1:nr)
{
	s11=d1[[i1]]
	s12=d2[[i1]]
	s13=d3[[i1]]
	for (i2 in 1:i1)
	{
		s21=d1[[i2]]
		s22=d2[[i2]]
		s23=d3[[i2]]
		
		temp=0
		
		n1=dim(s11)[2]
		n2=dim(s21)[2]
		
		if (!is.null(n1)&!is.null(n2)>0)
		{
		for (j1 in 1:n1)
		{
			for (j2 in 1:n2)
				temp=temp+sum(s11[,j1]*s21[,j2]*xx$wt)*sum(s12[,j1]*s22[,j2]*xx$wt)*sum(s13[,j1]*s23[,j2]*xx$wt)
		}
		}
		pk[i1,i2]=temp
		pk[i2,i1]=temp
	
	}
}


return(pk)
}