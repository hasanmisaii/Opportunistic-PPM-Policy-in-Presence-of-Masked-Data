   # Weibull Distribution 
   alpha = c(0.5, 0.25, 0.7)
   beta = c(0.25, 0.5, 1.5)
   
   beta[1] * gamma(1+1/alpha[1])
   beta[2] * gamma(1+1/alpha[2])
   beta[3] * gamma(1+1/alpha[3])
   
   
   J = length(alpha)
   
   
   Rel  <- function(t) {
     
     fu <- function(u) {
         dweibull(u,alpha[1],beta[1])*(1-pweibull(u,alpha[2],beta[2]))*
         (1-pweibull(u,alpha[3],beta[3])) + 
         dweibull(u,alpha[2],beta[2])*(1-pweibull(u,alpha[1],beta[1]))*
         (1-pweibull(u,alpha[3],beta[3])) +
         dweibull(u,alpha[3],beta[3])*(1-pweibull(u,alpha[2],beta[2]))*
         (1-pweibull(u,alpha[1],beta[1]))
     }
     
    
       #R = integrate(fu, lower = t, upper = Inf)$value
       
       R=function(t){sapply(t, function(t) {
          integrate(fu, lower = t, upper = Inf)$value
       })}
       
       return(R(t))
    
   }
   
   Rel <- Vectorize(Rel)
    
   Rel(-Inf)
   Rel(1000)
   
   # 
   
  x= seq(0,1,0.01)
   
  Rel(x)
   
  plot(Rel, 0, 1, type = "l" )
  
  plot(Rel,0,1, lty = 1, col = 1,lwd=3, xlab="Time", ylab="Reliability function", cex.lab=1.55, cex.axis=1.55, cex.main=1.55, cex.sub=1.55)
  
  h <- function(t) {
    lambda = 0
    for (j in 1:J) {
     lambda = lambda + dweibull(t, alpha[j], beta[j])/(1-pweibull(t, alpha[j], beta[j]))
    }
    return(lambda)
  }
  
  plot(h, 0, 0.5, type = "l" )
  plot(h,0,0.5, lty = 1, col = 1,lwd=3, xlab="Time", ylab="Failure rate function", cex.lab=1.55, cex.axis=1.55, cex.main=1.55, cex.sub=1.55)
  
  
   ############################ Maintenance 
  RR = 1
  
  f_TK  <- function(t,l) {
     for (j in 1:J) {
       if(j != l) {
           RR=RR*(1-pweibull(t,alpha[j],beta[j]))
       }
        else next;
     }
   dweibull(t,alpha[l],beta[l])*RR
        
  }
  
  f_TK(2,1)
  
  PMl = c(1,1,1)
  S = 0
  
  p_j  <- function(k,tau,jj) {
     f_TK1 <- function(l) {sapply(l, function(l){
        integrate(f_TK, lower = (k-1)*tau, upper = k*tau, l=l)$value
     })
     }
     for (j in 1:J) {
        S = S + PMl[j]*f_TK1(j)
     }
     PMl[jj]*f_TK1(jj)/S
  }
  
  p_j(10,1,1)
  
 
  
  K = 100; tau = 0.9
  p_j <- Vectorize(p_j)
  
  p_j(1,(1:K)*tau,3)
  
  p_j(1:K,tau,3)
  
  
  c_ins = 1; c_r = c(2,4,5); c_pm = c(1.5,2.5,3.5)
  rho = 0.1
 
  
  
  PCM  <- function(tau,k,l) {
    RR = 1
    for (j in 1:J) {
      if(j != l) {
        RR=RR*(1-pweibull(k*tau,alpha[j],beta[j]))
      }
      else next;
    }
     (pweibull(k*tau,alpha[l],beta[l])-pweibull((k-1)*tau,alpha[l],beta[l])) * RR
  }
  
  PCM(0.02,1:100,1)
  
  
  
  PPM  <- function(tau,k,l,rho) {
    pp = 0
    for (i in 1:J) {
      if (i!=l) {
      for (j in 1:J) {
        if(j!=l & j!=i){
          pp = pp + (pweibull(k*tau,alpha[i],beta[i])-pweibull((k-1)*tau,alpha[i],beta[i]))*(1-pweibull(k*tau,alpha[j],beta[j]))
        }
        }
      } 
    }
     (1-pweibull(k*tau,alpha[l],beta[l]))*I(p_j(k,tau,l)>rho)*pp
  }
  
  PPM(0.02,1:100,3,rho)
  
  EC <- function(tau,rho) {
     ec = c(0,0,0)
     for (j in 1:J) {
        for (k in 1:K) {
           ec[j] = ec[j] + (k*c_ins + c_r[j])*PCM(tau,k,j) + (k*c_ins + c_pm[j])*PPM(tau,k,j,rho)
        } 
     }
     return(ec)
  }
  
  EC(0.1,0.2)
  
  rho = seq(0.1,0.9,0.1)
  
  EC <- Vectorize(EC)
  tau = 0.1
  EC(tau,rho)
  
  rho = 0.4
  
  tau = seq(0.000001,2,0.1)
  EC(tau,rho)
  
  
  ET <- function(tau,rho) {
     et = c(0,0,0)
     for (j in 1:J) {
        for (k in 1:K) {
           et[j] = et[j] + k*tau*(PCM(tau,k,j)+PPM(tau,k,j,rho))
        } 
     }
     return(et)
  }
  
  ET(0.2,rho)
  
  
  
  
  ET <- Vectorize(ET)
  rho = seq(0,1,0.1)
  tau = 0.9
  EC(tau,rho)
  
 
   CT <- function(tau,rho) {
     ct = 0
     for (j in 1:J) {
       ct = ct + EC(tau,rho)[j]/ET(tau,rho)[j]
     }
      return(ct)
   }
   
   CT(0.02,rho)
   
   K=20
   rho = seq(0,1,0.1)
   tau1 = 1.11
   
   
   CT <- Vectorize(CT)
   C = CT(tau1,rho)
   C
   
   plot(rho,C,type = "l")
   
   
   tau1 = tau1[is.na(C) == F & C != Inf]
   tau1
   
   C = C[is.na(C) == F & C != Inf]
   
   tau1[which.min(C)]
   
   plot(tau1,C,type = "l")
   
   
   tau = tau1[which.min(C)]
   
   ####### Optimization of tau 
   K=20
   tau1 = seq(0.01,2,0.1)
   rhoo = seq(0,1,0.1)
   T = matrix(0,length(tau1),3)
   
   for (i in 1:length(tau1)) {
     tau = tau1[i]
     crho= c(rep(0,length(rhoo)))

       crho = CT(tau,rhoo)
    
     ru = rhoo[which.min(crho)]
     cru = crho[which.min(crho)]
     
     T[i,]=c(tau, ru, cru)
   }
   
   T
   
   ##########################
   

   
   
   rhoo = seq(0,1,0.1)
   crho= c()
   
   tau = 1.11
   
   
   for (i in 1:length(rhoo)) {
    crho[i] = CT(tau,rhoo[i])
   }
   
   crho
   rhoo[which.min(crho)]
   crho[which.min(crho)]

   plot(rhoo,crho, type = "l")   
   
   
   plot(rhoo, crho, type = "l", lty = 1, xaxt = "n", yaxt="n",col = 1,lwd=3,xlab=expression(rho), ylab="Rate of cost", cex.lab=1.55, cex.axis=1.55, cex.main=1.55, cex.sub=1.55)
   
   
   axis(1,at=c(seq(0,1,0.1)),cex.axis=2)
   axis(2,at=c(0,round(c(crho[1],crho[2],crho[3],crho[4], crho[which.min(crho)],crho[6],crho[7],crho[8],crho[10]),3)),cex.axis=1.2)
   
   
   abline(v=rhoo[which.min(crho)],lwd=2,col="red")
   abline(h=crho[which.min(crho)],lwd=2,col="red")
   
   yy = seq(tau,K*tau,K*tau/10)
   plot(c(0,yy), crho, type = "l")
   
   library("plot3D")
   
   scatter3D(x = rhoo, y = c(0,yy), z = crho,
             phi = 0, bty = "g",  type = "h",
             ticktype = "detailed",
             pch = 19, cex = 0.5, xlab ="p", ylab = "Time", zlab = "", zlim = NULL,
             cex.lab=1.55, cex.axis=1.5, cex.main=1.55, cex.sub=1.5)
   
   
   


      
   
  