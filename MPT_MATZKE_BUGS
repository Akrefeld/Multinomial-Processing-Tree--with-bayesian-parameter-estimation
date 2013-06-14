Multinomial-Processing-Tree--with-bayesian-parameter-estimation
===============================================================
model{
#stimulus loop
for (s in 1:S){
  
  #mpt category probabilities for complete presented stimuli
  # theta[,1] is D, theta[,2] is d, and theta[,3] is N
  pC[s,1] <- theta[s,1]*theta[s,2] + theta[s,1]*(1-theta[s,2])*(1-c) + (1-theta[s,1])*a*(1-c) # category probability Complete presentation -> response "complete"
  pC[s,2] <- theta[s,1]*(1-theta[s,2])*c + (1-theta[s,1])*a*c # category probability Complete presentation -> response "halved"
  pC[s,2] <- (1-theta[s,1])*(1-a)# category probability Complete presentation -> response "new"
  
  # Given the category probabilities, the data for word pairs 
  # (i.e., category counts for pairs) follow a multinomial distribution       
  n1[s,1:3] ~ dmulti(pC[s,1:3],C[s])

  
  pH[s,1] <- theta[s,1]*theta[s,2] + theta[s,1]*(1-theta[s,2])*c + (1-theta[s,1])*a*c # category probability halved presentation -> response "halved"
  pH[s,2] <- theta[s,1]*(1-theta[i,2])*(1-c) + (1-theta[s,1])*a*(1-c) # category probability halved presentation -> response "complete"
  pH[s,3] <- (1-theta[s,1])*(1-a) # category probability halved presentation -> response "new"
  
  n2[s,1:3] ~ dmulti(pH[s,1:3],H[s])
  
  pN[s,1] <- theta[s,3] + (1-theta[s,3])*(1-a) # category probability NO presentation -> response "new"
  pN[s,2] <- (1-theta[s,3])*a*(1-c) # category probability NO presentation -> response "complete"
  pN[s,3] <- (1-theta[s,3])*a*c
  
  n3[s,1:3] ~ dmulti(pN[s,1:3], N[s])
  # Build parameters c, r, and u parameters on probit scale
  for (p in 1:P){
    #probit transformation
    theta[s,p] <- phi(theta.probit[s,p])
    #############
    #here must be the regression equation to predict the parameter
    #scale parameter to normal distribution
    theta[S, p] <- mu[p] + xi.stim[p]*delta.stim.raw[s,p]
  }
  delta.stim.raw[s,1:P] ~dmnomr(mu.delta.raw[1:P], T.prec.stim[1:P,1:P])
  # Hyperpriors
  for(p in 1:P){
    mu.delta.raw[p] <- 0
    mu[p] ~ dnorm(0,1)    	
    xi.stim[p] ~ dunif(0,100)				
}
  
  T.prec.stim[1:P,1:P] ~ dwish(W[1:P,1:P],df)
  df<-P+1
  T.stim[1:P,1:P] <- inverse(T.prec.stim[1:P,1:P])
  
  # Scale sigma's and compute parameter correlations
  for(p in 1:P){
    for(p.prime in 1:P){
      # Off-diagonal elements of S
      rho.stim[p,p.prime] <- T.stim[p,p.prime]/sqrt(T.stim[p,p]*T.stim[p.prime,p.prime])
    }
    # Diagonal elements of S
    sigma.stim[p]<-xi.stim[p]*sqrt(T.stim[p,p])
  } 
  a~ dbeta(2,2)
  c~ dbeta(2,2)
}
