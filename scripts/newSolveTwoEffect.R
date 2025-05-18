require(nleqslv)
dnorminv<-function(y,s=1) sqrt(-2*s*log(sqrt(2*pi)*y))
pPoisConv <- function(t,lambda,norm.sd,alphal=1,risk.allele=FALSE){
  ## recover()
  ## t = threshold position
  ## lambda = genetic mean
  ## norm.sd = standard deviation of normal component
  ## alphal = effect size of large effect allele
  ## risk.allele = T/F risk or protective allele
  my.range <- seq(qpois(0,lambda),qpois(1-1e-8,lambda))
  gen.dist <- dpois(my.range,lambda=lambda)
  prevs <- sapply(t,function(t)pnorm(t,alphal*(my.range+risk.allele-lambda),norm.sd,lower.tail=FALSE))
  gen.dist%*%prevs
}
dPoisConv <- function(t,lambda,norm.sd,alphal=1,risk.allele=FALSE){
  ## recover()
  ## t = threshold position
  ## lambda = genetic mean
  ## norm.sd = standard deviation of normal component
  ## alphal = effect size of large effect allele
  ## risk.allele = T/F risk or protective allele
  my.range <- seq(qpois(0,lambda),qpois(1-1e-8,lambda))
  gen.dist <- dpois(my.range,lambda=lambda)
  dens <- sapply(t,function(t)dnorm(t,alphal*(my.range+risk.allele-lambda),norm.sd))
  gen.dist%*%dens
}
dltstarDiff <- function(log.dl,log.tstar,gs,L,C,al,h2,Vas,ft){
  ## recover()
  dl <- exp(log.dl)
  tstar <- exp(log.tstar)
  lambda <- 2*L*(1-gs)*u / (dl*C)
  Val <- al^2 * lambda
  Va <- Vas + Val
  Ve <- (1-h2)/h2*Va
  Vt <- Va + Ve
  Vnorm <- Vas + Ve
  ft2 <- dPoisConv(tstar,lambda,sqrt(Vnorm),al)
  Fta <- pPoisConv(tstar,lambda,sqrt(Vnorm),al,TRUE)
  Ft <- pPoisConv(tstar,lambda,sqrt(Vnorm),al,FALSE)
  dl2 <- Fta - Ft
  output <- c(dl2 - dl,ft2-ft)
  ## print(output)
  output
}
dlDiff <- function(thr,L,gs,u,dl,C,al) {
  lambda <- 2*L*(1-gs)*u / (dl*C)
  Val <- al^2 * lambda
  Va <- Vas + Val
  Ve <- (1-h2)/h2*Va
  Vt <- Va + Ve
  Vnorm <- Vas + Ve
  ft - dPoisConv(thr,lambda,sqrt(Vnorm),al)
}
solveTwoEffect1D <- function(bt,
                             Ne,
                             as,
                             al,
                             L,
                             gs,
                             h2 = NULL,
                             Ve = NULL,
                             u,
                             C,
                             Bval=1,
                             init.dl=NULL,
                             init.tstar=NULL
) {
  ## recover()
  amean <- as*gs + al*(1-gs)
  bs <- (bt*amean - (1-gs)*al) / (gs*as)
  if(bs < 0 ) {
    return(c(dl=NA,tstar=NA,lambda=NA,prev=NA))
  }
  gammaS <- 0.5*log ((1+bs)/(1-bs))
  ft <- (4*Ne*C*as)^(-1) * log ((1+bs)/(1-bs))
  Vas <- 2*L*gs*u*as*bs / (C*ft)
  
  log.init.dl <- log(init.dl)
  log.init.tstar <- log(init.tstar)
  
  out <- nleqslv(
    x = c(log.init.dl,log.init.tstar) ,
    fn = function(X) dltstarDiff(X[1],X[2],gs,L,C,al,h2,Vas,ft) ,
    control = list(maxit=500)
  )
  if(out$termcd == 1){
    dl <- exp(out$x[1])
    tstar <- exp(out$x[2])
    lambda <- 2*L*(1-gs)*u / (dl*C)
    Val <- al^2 *lambda
    Va <- Vas + Val
    Ve <- (1-h2) / h2 * Va
    Vt <- Va + Ve
    Vnorm <- Vas + Ve
    this.prev <- pPoisConv(tstar,lambda,sqrt(Vnorm),al,FALSE)
    phit <- ft*sqrt(Vt)
    norm.prev <- 1-pnorm(dnorminv(phit))
    return(c(al=al,dl=dl,tstar=tstar,lambda=lambda,prev=this.prev,Vt=Vt,Vas=Vas,Val=Val,Ve=Ve,gammaS=gammaS,norm.prev=norm.prev))
  } else {
    return(c(al=NA,dl=NA,tstar=NA,lambda=NA,prev=NA,Vt=NA,Vas=NA,Val=NA,Ve=NA,gammaS=NA,norm.prev=NA))
  }
}
