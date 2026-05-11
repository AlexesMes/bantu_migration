#Function to extract quantiles from MCMC output 
extract <- function(x)
{
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^a\\[',colnames(tmp))]
  qta = apply(tmp2, 2, quantile, prob=c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
  return(qta)
}

#Function to plot posterior bars
post.bar <- function(x, i, h, barcolours=c("skyblue","dodgerblue","darkblue","darkgreen"))
{
  rect(xleft = x[2], xright = x[6], ybottom = i - h/5, ytop = i + h/5, border = NA, col = barcolours[1]) # 95% interval rectangle
  segments(x[2], i-h/3.5, x[2], i+h/3.5, lwd = 2, col = barcolours[1]) # horizontal ticks for 95%
  segments(x[6], i-h/3.5, x[6], i+h/3.5, lwd = 2, col = barcolours[1]) # horizontal ticks for 95%
  
  rect(xleft=x[3], xright=x[5], ybottom=i-h/3, ytop=i+h/3, border=NA, col=barcolours[2]) #50% interval
  segments(x[3], i-h/2.5, x[3], i+h/2.5, lwd = 2, col = barcolours[2])   # horizontal ticks for 50%
  segments(x[5], i-h/2.5, x[5], i+h/2.5, lwd = 2, col = barcolours[2])   # horizontal ticks for 50%
  
  points(x[4], i, pch = 16, col = barcolours[3], cex = 2) #posterior median 
}

#Function to plot posterior estimates of pairwise difference of arrival times
diffDens  <- function(x,y,prob=0.9,...)
{
  require(coda)
  z  <- x - y 
  nsample = length(z)
  left  <- c(HPDinterval(mcmc(z),prob=prob)[1],0)
  right  <- c(0,HPDinterval(mcmc(z),prob=prob)[2])
  plotRight=plotLeft=TRUE
  if (any(right<0)){left[2]=right[2];plotRight=FALSE}
  if (any(left>0)){right[1]=left[1];plotLeft=FALSE}
  dens = density(z)
  hpdi.left.x = dens$x[which(dens$x>=left[1]&dens$x<=left[2])]
  hpdi.left.y = dens$y[which(dens$x>=left[1]&dens$x<=left[2])]
  hpdi.right.x = dens$x[which(dens$x>=right[1]&dens$x<=right[2])]
  hpdi.right.y = dens$y[which(dens$x>=right[1]&dens$x<=right[2])]
  plot(dens$x,dens$y,type='n',xlab='Years',ylab='Probability Density',axes=FALSE,...)
  if(plotLeft){polygon(x=c(hpdi.left.x,rev(hpdi.left.x)),y=c(hpdi.left.y,rep(0,length(hpdi.left.y))),border=NA,col='lightpink')}
  if(plotRight){polygon(x=c(hpdi.right.x,rev(hpdi.right.x)),y=c(hpdi.right.y,rep(0,length(hpdi.right.y))),border=NA,col='lightblue')}
  lines(dens)
  abline(v=0,lty=2,lwd=1.5)
  axis(1,cex.axis=0.85,padj=-0.5)
  xlim = range(axTicks(1))
  axis(1,at=seq(xlim[1],xlim[2],100),labels=NA,tck=-.01)
}

orderPPlot  <- function(x, name.vec, addCoef=NULL)
{
  order.prob <- matrix(NA, nrow=ncol(x), ncol=ncol(x), dimnames=list(name.vec, name.vec))
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      if (i>j)
      {
        order.prob[i,j] = round(sum(x[,i] < x[,j])/nrow(x),2)
      }
    }
  }
  corrplot(order.prob,
           is.corr=FALSE,
           diag=FALSE,
           col=COL2('RdBu'),
           addCoef.col = addCoef,
           number.cex=0.8,
           method='color',
           type='lower',
           tl.col = "black",
           tl.cex = 0.7)
}