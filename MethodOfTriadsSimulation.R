# Simulate data sets and store
rm(list=ls())
library(MASS)
set.seed(125)

# Set up storage
n.sim = 500
n.substudy = 1000
datasets.substudy = array(0, c(n.substudy, 4, n.sim))
# Col1: biomarker (M). Col2 = FR (R). Col3 = FFQ (W). Col4 = TDI (X)

# Mean and sd parameters drawn from Table 3 Prentice et al 2011b
# Truth X~N(mu.x, sig.x^2)
mu.x = log(2141)
sig.x = log(14)
for(i in 1:n.sim)
{
  datasets.substudy[,4,i] = rnorm(n.substudy, mu.x, sig.x)
}


# Simulation settings:
# Considering lambda.wx = 0.8 and sig.m = alpha*sig.x for alpha = 0.5
# M, R and W are MVN_3; mu.m = 0 + lambda.mx*X; mu.r = 0 + lambda.rx*X; mu.w = mu.w + lambda.wx*X; 
# Sigma: covariance matrix of M, R and W errors. Key assumption is this is diagonal!
mu.m = 0
mu.r = 0
mu.w = log(1485)
lambda.mx = 1
lambda.rx = 0.8 #0.95, 0.8
lambda.wx = 0.5 #0.9, 0.5
lambda = matrix(c(lambda.mx, lambda.rx, lambda.wx), nrow=3)
alpha = 0.5; sig.m = alpha*sig.x           # Consider best case scenario 
sig.r = log(21)
sig.w = log(28)
rho.mr = 0
rho.mw = 0
rho.rw = 0 #0, 0.1, 0.8
Sigma = matrix(NA, 3,3)
Sigma[1,] = c(sig.m^2, sig.m*rho.mr*sig.r, sig.m*rho.mw*sig.w)
Sigma[2,] = c(sig.m*rho.mr*sig.r, sig.r^2, sig.r*rho.rw*sig.w)
Sigma[3,] = c(sig.m*rho.mw*sig.w, sig.r*rho.rw*sig.w, sig.w^2)

for(i in 1:n.sim)
{
  datasets.substudy[,1:3,i] = mvrnorm(n.substudy, c(mu.m, mu.r, mu.w), Sigma) + t(lambda%*%datasets.substudy[,4,i])
}

# True validity coefficient!
vc.true = sqrt(1/(1 + sig.w^2/(lambda.wx^2*sig.x^2)))


###### METHOD OF TRIADS APPROACH TO ASSESSING VALIDITY OF DIETARY INTAKE INSTRUMENT #############

vc.est.store = matrix(as.double(NA), nrow=n.sim, ncol=1)

for(i in 1:n.sim)
{
  # Validity coefficient = rho_MW*rho_WR/rho_MR
  # Col1: biomarker (M). Col2 = FR (R). Col3 = FFQ (W). Col4 = TDI (X)
  dat = datasets.substudy[,-4,i]
  rho_MW = cor(dat[,1], dat[,3])
  rho_WR = cor(dat[,2], dat[,3])
  rho_MR = cor(dat[,1], dat[,2])
  vc.wx = sqrt(rho_MW*rho_WR/rho_MR)
  vc.est.store[i,] = vc.wx
}


## Plot some results
library(ggplot2)
dat = data.frame(vc.est.store)
names(dat) = c("vc.est")

ggplot()+
  geom_histogram(aes(vc.est), data=dat, binwidth = 0.005)+
  labs(x=expression(VC[WX]), y="Count")+
  geom_vline(xintercept=vc.true, colour="#000099", linetype="dashed")



