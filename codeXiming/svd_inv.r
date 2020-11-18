
svd_inv=function(x,svdtol=sqrt( .Machine$double.xmin ))
{
	svdh <- svd( x )
  if( min(svdh$d) < max(svdh$d)*svdtol )
    svdh$d <- svdh$d + max(svdh$d)*svdtol
  svdh$v %*% (t(svdh$u)/svdh$d)
}