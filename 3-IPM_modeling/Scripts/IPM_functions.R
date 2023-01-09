################################################################################
###### File for IPM and other analysis functions ###############################
################################################################################
year <- function(String){
  a <- str_split(String, "-")
  years <- c()
  for(i in 1:length(a)){
    years[i] <- a[[i]][1]
  }
  return(years)
}

SmithsonVerkuilen2006=function(yi){
  n=length(yi)
  yo=(yi * (n-1) + 0.5) / n
  return(yo)
}

Inv.SmithsonVerkuilen2006=function(yo){
  n=length(yo)
  yi=(n*yo - 0.5)/(n-1)
  return(yi)
}


inv.logit <- function(x) {
  exp(x)/(1+exp(x))
}

################################################################################
###### The growth model ########################################################
################################################################################

growth_function <- function(y, x, g.int, g.slp, g.var) {
  mu <- (g.int + g.slp * x)
  sig <- sqrt(g.var)
  dnorm(y, mean=mu, sd=sig)
}

################################################################################
###### The survival model ######################################################
################################################################################


survival_function <- function(x, s.int, s.slp,interval_years) {
  #, s.int, s.slp)
  Line <- s.slp * x + s.int
  u <- inv.logit(Line)^(1/interval_years) # From binomial glm, 
                                          # link = logit function
  return(u)
}

survival_regional_function <- function(x, s.int, s.slp) {
  Line <- s.slp * x + s.int
  u <- inv.logit(Line)
  return(u)
}

################################################################################
###### The reproduction model ##################################################
################################################################################

reproduction_function <- function(y, x, rec, range, rec.size=rec.size, ind=F) {
  # 2*sqrt((10^-2.5)/pi)  # approx 6 cm  diameter coral
  if (ind == TRUE){
    out <- rep(0,length(x)) 
    out[1:range] <- rec/range} 
  
  if (ind == FALSE){
    out <-  (10^x) * rec
    out[x < rec.size | y >= rec.size] <- 0 
    #if x is below recruitment size then it doesn't count 
  }
  return(out)
}


reproduction_function_josh <- function(y, x, rec, rec.size=rec.size) {
  out <-  (10^x) * rec
  out[x < rec.size | y >= rec.size] <- 0
  return(out)
}

################################################################################
########## ... Running the IPM ... #############################################
################################################################################

bigmatrix <- function(ModelParameters, sensitivity = T) {
  # Take in the model parameters as a list and runs the IPM
  n          <- ModelParameters$n
  y          <- ModelParameters$y
  rec        <- ModelParameters$rec
  rec.size   <- ModelParameters$rec.size
  g.int      <- ModelParameters$g.int
  g.slp      <- ModelParameters$g.slp
  g.var      <- ModelParameters$g.var
  s.int      <- ModelParameters$s.int
  s.slp      <- ModelParameters$s.slp
  iy         <- ModelParameters$Interval_Years
  delta_size <- ModelParameters$ds
  
  
  Reproduction_kernel <- delta_size * outer(y, y, 
                                            reproduction_function_josh,rec=rec, 
                                            rec.size=rec.size)
 
  
  Growth_kernel <- delta_size * outer(y, y, growth_function, g.int=g.int, 
                                      g.slp=g.slp, g.var=g.var)
  
  
  Survival_kernel <- survival_function(y, s.int=s.int, s.slp=s.slp,
                                       interval_years = iy)
  
  
  P <- Growth_kernel
  i <- 1:n
  
  P[,i] <- Growth_kernel[,i]*Survival_kernel[i]
  
  K <- P + Reproduction_kernel
  
  lam         <- Re(eigen(K)$values[1])
  w           <- Re(eigen(K)$vectors[,1])
  stabledist  <- w/sum(w) # right eigenvector/stable distribution
  v.eigen     <- Re(eigen(t(K))$vectors[,1]) #left eigenvector
  reprodvalue <- v.eigen/v.eigen[1] #reproductive value
  
  if(sensitivity==T){
 #compute sensitivity and elasticity matrices using the eigenvectors and lambda
    #v.dot.w <- sum(stabledist * reprodvalue) * delta_size
    v.dot.w <- sum(stabledist * reprodvalue * delta_size) # JOSH CHANGE
    sens <- outer(reprodvalue,stabledist) / v.dot.w
    elas <- matrix(as.vector(sens) * as.vector(K) / lam, nrow = n)
    K.elas <- sens * (K) / lam * delta_size #overall kernel elasticity
    P.elas <- sens * (P) / lam * delta_size #elas for growth/survival
    R.elas <- sens * (Reproduction_kernel) / lam * delta_size #elas for reproduction
    # Don't think it's possible to do this for G and S, looking into it -- JOSH
    #G by matrix - seems like this a possible analogue - I'm going to try it a bit...
    G.elas <- sens * (Growth_kernel) / lam * delta_size # TOM CHANGE?
    
    return(list(K=K, Gk=Growth_kernel,Sk=Survival_kernel,Pk=P,lam=lam, 
                w=stabledist, v=reprodvalue,
                v.dot.w=v.dot.w, sens=sens, 
                K.elas=K.elas, P.elas=P.elas, R.elas=R.elas, #elas contributions
                eK=sum(K.elas), eP=sum(P.elas), eR=sum(R.elas) #eK = eP + eR = 1
                )) # JOSH CHANGE
  }
  
  if(sensitivity==F){
    return(list(K=K, Gk=Growth_kernel,Sk=Survival_kernel,Pk=P,lam=lam, 
                w=stabledist, v=reprodvalue))
  }
}

################Regional IPM###################

bigmatrix_regional <- function(ModelParameters, sensitivity = T) {
  # Take in the model parameters as a list and runs the IPM
  n          <- ModelParameters$n
  y          <- ModelParameters$y
  rec        <- ModelParameters$rec
  rec.size   <- ModelParameters$rec.size
  g.int      <- ModelParameters$g.int
  g.slp      <- ModelParameters$g.slp
  g.var      <- ModelParameters$g.var
  s.int      <- ModelParameters$s.int
  s.slp      <- ModelParameters$s.slp
  delta_size <- ModelParameters$ds
  
  
  Reproduction_kernel <- delta_size * outer(y, y, 
                                            reproduction_function_josh,rec=rec, 
                                            rec.size=rec.size)
  
  
  Growth_kernel <- delta_size * outer(y, y, growth_function, g.int=g.int, 
                                      g.slp=g.slp, g.var=g.var)
  
  
  Survival_kernel <- survival_regional_function(y, s.int=s.int, s.slp=s.slp)
  
  
  P <- Growth_kernel
  i <- 1:n
  
  P[,i] <- Growth_kernel[,i]*Survival_kernel[i]
  
  K <- P + Reproduction_kernel
  
  lam         <- Re(eigen(K)$values[1])
  w           <- Re(eigen(K)$vectors[,1])
  stabledist  <- w/sum(w) # right eigenvector/stable distribution
  v.eigen     <- Re(eigen(t(K))$vectors[,1]) #left eigenvector
  reprodvalue <- v.eigen/v.eigen[1] #reproductive value
  
  if(sensitivity==T){
    #compute sensitivity and elasticity matrices using the eigenvectors and lambda
    #v.dot.w <- sum(stabledist * reprodvalue) * delta_size
    v.dot.w <- sum(stabledist * reprodvalue * delta_size) #updated based on Josh's rec
    sens <- outer(reprodvalue,stabledist) / v.dot.w
    elas <- matrix(as.vector(sens) * as.vector(K) / lam, nrow = n)
    K.elas <- sens * (K) / lam * delta_size #overall kernel elasticity
    P.elas <- sens * (P) / lam * delta_size #elas for growth/survival
    R.elas <- sens * (Reproduction_kernel) / lam * delta_size #elas for reproduction
    
    return(list(K=K, Gk=Growth_kernel,Sk=Survival_kernel,Pk=P,lam=lam, 
                w=stabledist, v=reprodvalue,
                v.dot.w=v.dot.w, sens=sens,
                K.elas=K.elas, P.elas=P.elas, R.elas=R.elas, #elas contributions
                eK=sum(K.elas), eP=sum(P.elas), eR=sum(R.elas) #eK = eP + eR = 1
                ))
  }
  
  if(sensitivity==F){
    return(list(K=K, Gk=Growth_kernel,Sk=Survival_kernel,Pk=P,lam=lam, 
                w=stabledist, v=reprodvalue))
  }
}

