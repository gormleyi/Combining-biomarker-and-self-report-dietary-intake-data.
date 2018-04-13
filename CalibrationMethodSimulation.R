# Simulate data sets and store
rm(list=ls())
library(MASS)
set.seed(125)

# Set up storage
n.sim = 500
n.sub = 1000
n.study = 10000
datasets.sub = array(0, c(n.sub, 3, n.sim))
datasets.study = array(0, c(n.study, 3, n.sim))
# Col1: biomarker (M). Col2 = FFQ (W). Col3 = TDI (X)

# Mean and sd parameters drawn from Table 3 Prentice et al. (2011)
# Truth X~N(mu.x, sig.x^2)
mu.x = log(2141)
sig.x = log(14)
for(i in 1:n.sim)
{
  datasets.sub[,3,i] = rnorm(n.sub, mu.x, sig.x)
  datasets.study[,3,i] = rnorm(n.study, mu.x, sig.x)
}

# Simulation settings:
# Considering lambda.wx = 0.1, 0.5, 0.8 and sig.m = alpha*sig.x for alpha = 0.5, 1, 2.
# (M, FFQ) are MVN_2; mu.m = 0 + lambda.mx*X; mu.w = mu.w + lambda.wx*X. 
# Sigma: rho.mw = 0, 0.1, 0.8.
alpha = 0.5; sig.m = alpha*sig.x
sig.w = log(28)
rho.mw = 0
Sigma = matrix(c(sig.m^2, sig.m*rho.mw*sig.w, sig.m*rho.mw*sig.w, sig.w^2), nrow=2)
mu.m = 0
mu.w = log(1485)
lambda.mx = 1
lambda.wx = 0.8
lambda = matrix(c(lambda.mx, lambda.wx), nrow=2, ncol=1)
for(i in 1:n.sim)
{
  datasets.sub[,1:2,i] = mvrnorm(n.sub, c(mu.m, mu.w), Sigma) + t(lambda%*%datasets.sub[,3,i])
  datasets.study[,1:2,i] = mvrnorm(n.study, c(mu.m, mu.w), Sigma) + t(lambda%*%datasets.study[,3,i])
}

# True parameters! X|W deriveable via Bayes Theorem
Psi.tilde = (sig.w^2*sig.x^2)/(sig.w^2 + lambda.wx^2*sig.x^2)
beta0.tilde = (sig.w^2*mu.x - sig.x^2*lambda.wx*mu.w)/(sig.w^2 + lambda.wx^2*sig.x^2)
beta1.tilde = (lambda.wx*sig.x^2)/(sig.w^2 + lambda.wx^2*sig.x^2)

# Calibration. Exact parameters from regressing M|W (after obtaining their marginal)
mu_m = mu.m + lambda.mx*mu.x
mu_w = mu.w + lambda.wx*mu.x
Omega = matrix(c(lambda.mx^2*sig.x^2 + sig.m^2, sig.m*rho.mw*sig.w + sig.x^2*lambda.mx*lambda.wx, sig.m*rho.mw*sig.w + sig.x^2*lambda.mx*lambda.wx, lambda.wx^2*sig.x^2 + sig.w^2),2,2)
denom = (sig.w^2 + lambda.wx^2*sig.x^2)
beta0 = (sig.w^2*mu_m -sig.x^2*lambda.wx*(lambda.mx*mu.w - lambda.wx*mu.m) - (sig.m*rho.mw*sig.w)*mu_w)/denom
beta1 = (lambda.mx*lambda.wx*sig.x^2 + sig.m*rho.mw*sig.w)/denom
# More concise calculation:
#beta0 = mu_m - Omega[1,2]*solve(Omega[2,2])*mu_w
#beta1 = Omega[1,2]*solve(Omega[2,2])

# When rho.mw=0, all are equal. Not the case when rho.mw \neq 0.
c(beta0.tilde, beta0)
c(beta1.tilde, beta1)


################################################################################
# Run calibration by regressing biomarker on FFQ from sub study
beta.est.store = matrix(as.double(NA), nrow=n.sim, ncol=2)
pred.store = matrix(as.double(NA), nrow=n.sim, ncol=n.study)
pred.diffs = matrix(as.double(NA), nrow=n.sim, ncol=n.study)

for(i in 1:n.sim)
{
  # Regress biomarker on FFQ from sub study
  res = lm(datasets.sub[,1,i]~datasets.sub[,2,i])
  beta.est.store[i,] = res$coef
  
  # Predict true intake from FFQ for study data
  pred.store[i,] = res$coef[1] + res$coef[2]*datasets.study[,2,i]
  pred.diffs[i,] = pred.store[i,] - datasets.study[,3,i]
}

# Plot some beta results
library(ggplot2)
dat = data.frame(beta.est.store[,2]/c(beta1.tilde))
names(dat) = c("BetaRatio")
ggplot()+
  geom_boxplot(aes(x=1, y=BetaRatio, fill="red", color="red"), data=dat)+
  labs(y = expression(hat(beta)[1]/tilde(beta)[1]), x="")+
  theme(legend.position="none", axis.title=element_text(size=12), axis.text=element_text(size=10))+
  geom_hline(yintercept=1, colour="#000099", linetype="dashed")


# Plot some prediction results
dat = data.frame(apply(abs(pred.diffs),1,mean)) # Mean absolute error
names(dat) = c("mae")

ggplot()+
  geom_histogram(aes(mae, fill = "blue"), binwidth=0.001, data=dat)+
  labs(x="Mean Absolute Prediction Error", y="Count")+
  theme(legend.position="none")


