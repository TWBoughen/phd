source('../scripts/igp_functions.R')

auto.mfrow <- function(nplots, setup=TRUE) {
  
  if(setup) {
    if(nplots <= 3) par(mfrow=c(1, nplots))
    else if(nplots <= 4)  par(mfrow=c(2,2))
    else if(nplots <= 6)  par(mfrow=c(2,3))
    else if(nplots <= 9)  par(mfrow=c(3,3))
    else if(nplots <= 12) par(mfrow=c(3,4))
    else if(nplots <= 16) par(mfrow=c(4,4))
    else if(nplots <= 20) par(mfrow=c(4,5))
    else if(nplots <= 25) par(mfrow=c(5,5))
    else if(nplots <= 30) par(mfrow=c(5,6))
    else if(nplots <= 36) par(mfrow=c(6,6))
    else if(nplots <= 42) par(mfrow=c(6,7))
    else if(nplots <= 49) par(mfrow=c(7,7))
    else if(nplots <= 56) par(mfrow=c(7,8))
    else if(nplots <= 64) par(mfrow=c(8,8))
    else {
      stop("Too many plots")
    }
  }
  else {
    nblankplots <- par("mfrow")[1] * par("mfrow")[2] - nplots
    if(nblankplots > 0)
      for(i in 1:nblankplots)
        plot_blank()
  }
}
plot_dir = function(dir, lower=1, lines=T){
  dat_list_names = konect_to_df_dir(dir,lower=lower, return_names=T)
  n = length(dat_list_names$dat)
  mx_dat = 0
  min_prob = 1
  
  for(i in 1:n){
    d = dat_list_names$dat[[i]]
    mx_dat = max(mx_dat, max(d[,1]))
    min_prob = min(min_prob, d[nrow(d),2]/sum(d[,2]))
  }
  k = 1:mx_dat
  m=1
  pk_ba = 2*m*(m+1) / (k*(k+1)*(k+2))
  pk_ua = (1/2)^k
  auto.mfrow(n)
  x = par()$mfrow[1]
  y = par()$mfrow[2]
  i=1
  for(dat in dat_list_names$dat){
    x_ax = 'n'
    y_ax='n'
    marge = c(0,0,0,0)
    if(i %% x ==1){
      marge[2] = 0
      y_ax = NULL
    }
    if(i>n-x){
      x_ax=NULL
    }
    if(i>(y-1)*x){
      marge[1] = 2
    }
    par(mar=marge)
    
    Fk = cumsum(dat$Freq)/sum(dat$Freq)
    plot(dat$x,1-c(0,Fk[-length(Fk)]), log='xy',pch=20,
         xlim = c(1,mx_dat), ylim=c(min_prob, 1),cex = 0.5,
         xaxt = x_ax, yaxt=y_ax,las=1, ylab='')
    if(lines){
      lines(k,1-c(0,cumsum(pk_ba)[-length(pk_ba)]),lty=2, col='darkgreen')
      lines(k,1-c(0,cumsum(pk_ua)[-length(pk_ua)]),lty=3, col='red')
      legend('bottomleft',
             legend=c(paste0('n=',sum(dat[,2])),strsplit(dat_list_names$names[i], 'out.')[[1]][2]),
             pch=c(NA,20),
             col=c(NA,1))
      legend('topright', legend = c('BA','UA'), lty=c(2,3), col=c('darkgreen','red'))
    }else{
      legend('bottomleft',
             legend=c(paste0('n=',sum(dat[,2])),strsplit(dat_list_names$names[i], 'out.')[[1]][2]),
             pch=c(NA,20),
             col=c(NA,1))
    }
    
    if(i>0){
      
    }else{
      legend('bottomleft',
             legend=c(paste0('n=',sum(dat[,2])),strsplit(dat_list_names$names[i], 'out.')[[1]][2]),
             pch=c(NA,20),
             col=c(NA,1))
    }
    
    i=i+1
  }
  return(dat_list_names$names)
}

