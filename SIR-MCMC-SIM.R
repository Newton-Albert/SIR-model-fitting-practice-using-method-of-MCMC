
set.seed(156829)
library("sfsmisc")
library("deSolve")
##################################################################################
##################################################################################
##################################################################################
SIRfunc=function(t, x, vparameters){
  S = x[1]  # the value of S at time t
  I = x[2]  # the value of I at time t
  R = x[3]  # the value of R at time t
  if (I<0) I=0 # this is a cross check to ensure that we always have sensical values of I
  
  with(as.list(vparameters),{
    npop = S+I+R   # the population size is always S+I+R because there are no births or deaths in the model
    dS = -beta*S*I/npop            # the derivative of S wrt time
    dI = +beta*S*I/npop - gamma*I  # the derivative of I wrt time
    dR = +gamma*I                  # the derivative of R wrt time
    out = c(dS,dI,dR)
    list(out)
  })
}

##################################################################################
##################################################################################
# Let's set up the model parameters, and some initial conditions at time t=0
##################################################################################
gamma = 1/3         # recovery period in days^{-1}
R0    = 1.5         # R0 of R of the disease
beta = gamma*R0

N = 1000   # population size
I_0 = 10     # number intially infected people in the population
S = N-I_0
R = 0

vt = seq(0,1000,1)
vparameters=c(gamma=gamma,beta=beta)
inits=c(S=S,I=I,R=R)
sirmodel = as.data.frame(lsoda(inits, vt, SIRfunc, vparameters))

##################################################################################
##################################################################################
##################################################################################
time = 0
vstate = c(S,I,R)

##################################################################################
# work out the possible transitions, and put them in the JxK matrix, lambda
##################################################################################
K = length(vstate)  # number of compartments
J = 2               # number of possible state changes
lambda = matrix(0,nrow=J,ncol=length(vstate))
lambda[1,] = c(-1,1,0)
lambda[2,] = c(0,-1,1)

##################################################################################
##################################################################################
# The list object, zstate, will have all the S,I, and R values over
# time for each realisation, along with the time, and realisation number
# (such that each element of the list is a vector of length 5)
#
# Using a list to store this information, instead of appending to a data frame at
# each step, is computationally more expedient
##################################################################################
zstate = list()
i = 1
while(vstate[2]>0&vstate[1]>0){
  zstate[[i]] = c(vstate,time)
  
  S = vstate[1]
  I = vstate[2]
  R = vstate[3]
  
  vec_p = c(beta*S*I/N
            ,gamma*I)
  
  delta_t = 1/sum(vec_p)
  
  vec_l = rpois(length(vec_p),vec_p*delta_t)
  vstate = vstate + vec_l%*%lambda
  
  vstate[vstate<0] = 0
  i = i+1
  time = time + delta_t
  
  if (i%%100==1) cat(i,time,vstate,"\n")
}

##################################################################################
##################################################################################
# the sapply functions here extract the elements of the zstate list into vectors
##################################################################################
vS = sapply(zstate, "[[", 1)
vI = sapply(zstate, "[[", 2)
vR = sapply(zstate, "[[", 3)
vtime = sapply(zstate, "[[", 4)

##################################################################################
##################################################################################
# plot the results of the stochastic realisation, 
# with the deterministic model overlaid
##################################################################################
par(mfrow=c(2,2))
plot(vtime,vS/N,lwd=3,type="l",xlab="Time, in Days",ylab="S",main=paste("N=",N," and I_0=",vI[1],sep=""),ylim=c(0,1))
lines(sirmodel$time,sirmodel$S/N,col=2,lwd=2)
legend("topright",legend=c("MCMC","Deterministic"),col=c(1,2),lwd=4,bty="n")
plot(vtime,vI/N,lwd=3,type="l",xlab="Time, in Days",ylab="I",main=paste("N=",N," and I_0=",vI[1],sep=""),ylim=c(0,0.1))
lines(sirmodel$time,sirmodel$I/N,col=2,lwd=2)
plot(vtime,vR/N,lwd=3,type="l",xlab="Time, in Days",ylab="R",main=paste("N=",N," and I_0=",vI[1],sep=""),ylim=c(0,1))
lines(sirmodel$time,sirmodel$R/N,col=2,lwd=2)

