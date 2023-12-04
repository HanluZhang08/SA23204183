## ---- include=FALSE-----------------------------------------------------------
library(ggplot2)
library(rpart)
library(rpart.plot)

## -----------------------------------------------------------------------------
data(mpg)
mydata = mpg[,c('displ','year','cyl','cty','hwy')]
summary(mydata)

## -----------------------------------------------------------------------------
cor(mydata, method = c("pearson"))

## -----------------------------------------------------------------------------
mydata$year <- as.factor(mydata$year)
t.test(displ~year, mydata, var.equal=T)

## -----------------------------------------------------------------------------
lm_model = lm(hwy~displ+cyl, mydata)
par(mfrow=c(2,2), mar = c(2, 2, 2, 2))
plot(lm_model)

## -----------------------------------------------------------------------------
tree_model = rpart(hwy~., data = mydata)
rpart.plot(tree_model)

## -----------------------------------------------------------------------------
set.seed(1)
n = 1000
u = runif(n)
generate_sample = sign(0.5-u)*log(1-2*abs(u-0.5))

## -----------------------------------------------------------------------------
# 绘制图像
hist(generate_sample, breaks=30, xlim=c(-10,10), prob=TRUE, 
     col=rgb(1,0,0,0.5), xlab="sample", ylab="Density",
     main = "Histogram of standard laplace")
y = seq(-10, 10, 0.1)
lines(y, 0.5*exp(-abs(y)))

## -----------------------------------------------------------------------------
set.seed(1)

generate_beta = function(a, b, n){
  k = 0 # counter for accept
  j = 0 # iterations
  beta_sample = numeric(n)
  while (k < n) {
    u = runif(1) # upper bound
    j = j + 1
    x = runif(1) # r.v. from g
    if (x^(a-1)*(1-x)^(b-1)>u){
      # accept x
      k = k + 1
      beta_sample[k] = x
    }
  }
  return(beta_sample)
}

## -----------------------------------------------------------------------------
# 生成
a = 3
b = 2
Beta_Sample = generate_beta(a,b,1000)

# 画图
hist(Beta_Sample, breaks=30, xlim=c(0,1), prob=TRUE, col=rgb(1,0,0,0.5), xlab="sample", 
     ylab="Density", main = "Histogram of Beta(3,2)")
y = seq(0, 1, 0.01)
lines( y, factorial(a+b-1)/(factorial(a-1)*factorial(b-1)) * y^(a-1) * (1-y)^(b-1) )

## -----------------------------------------------------------------------------
set.seed(1)

generate_fe = function(n){
  k = 0 # counter for accept
  j = 0 # iterations
  fe_sample = numeric(n)
  while (k < n) {
    u1 = runif(1, min = -1, max = 1)
    u2 = runif(1, -1, 1)
    u3 = runif(1, -1, 1)
    k = k + 1
    j = j + 1
    if ( abs(u3)>=abs(u2) & abs(u3)>=abs(u1) ){
      # accept u2
      fe_sample[k] = u2
    }else{
      fe_sample[k] = u3
    }
  }
  return(fe_sample)
}

## -----------------------------------------------------------------------------
fe_Sample = generate_fe(10000)

hist(fe_Sample, breaks=30, prob=TRUE, col=rgb(1,0,0,0.5), xlab="sample", 
     ylab="Density", main = "Histogram of rescaled Epanechnikov kernel")

## -----------------------------------------------------------------------------
my_sample = function(x, n = length(x), prob = rep(1/length(x), length(x))){
  cp = cumsum(prob)
  u = runif(n)
  random_sample = x[findInterval(u,cp)+1]
  return(random_sample)
}

## -----------------------------------------------------------------------------
set.seed(1)
x = c(0,1)
n = 10
s = sample(x, n, replace = TRUE)
ms = my_sample(x, n)
print(s)
print(ms)

## -----------------------------------------------------------------------------
set.seed(2)
x = letters
s = sample(x,replace = TRUE)
ms = my_sample(x)
print(s)
print(ms)

## -----------------------------------------------------------------------------
set.seed(3)
x = c(1:3)
n = 100
prob = c(0.2, 0.3, 0.5)
s = sample(x, n, prob, replace = TRUE)
ms = my_sample(x, n, prob)
count_sample <- as.vector(table(s))
count_my_sample <- as.vector(table(ms))
print(count_sample)
print(count_my_sample)

## ----eval = FALSE-------------------------------------------------------------
#  n = 1e6
#  K = 100
#  l = c(2, 1.6, 1)
#  d = 2
#  pihat = data.frame(matrix(NA, K, 3))
#  
#  for (i in 1:3) {
#    for (j in 1:K){
#      set.seed(j)
#      X = runif(n,0,d/2)
#      Y = runif(n,0,pi/2)
#      pihat[j,i] = 2*l[i]/d/mean(l[i]/2*sin(Y)>X)
#    }
#  }

## ----eval = FALSE-------------------------------------------------------------
#  # 计算方差
#  var = apply(pihat, 2, var)

## -----------------------------------------------------------------------------
# 函数参数：生成的随机数数量，是否要用对偶变量法
# 函数功能：如果要用对偶变量法，1/2的随机数要做变换；否则不用做变换，即简单MC法
MC_anti = function(n, antithetic = TRUE){
  u = runif(n/2, 0, 1)
  if(!antithetic) v = runif(n/2, 0, 1) else
    v = 1 - u 
  u = c(u, v)
  theta = mean(exp(u))
  theta
}

# 开始计算theta及其方差
K = 100
n = 10000
theta_simple = theta_anti = numeric(K)
for (i in 1:K) {
  set.seed(i)
  theta_simple[i] = MC_anti(n, antithetic = FALSE)
  theta_anti[i] = MC_anti(n, antithetic = TRUE)
}
var_simple = var(theta_simple)
var_anti = var(theta_anti)
print(var_simple)
print(var_anti)
print((var_simple-var_anti)/var_simple)

## ----echo = FALSE-------------------------------------------------------------
x <- seq(1, 6, .01)
g = x^2 / sqrt(2*pi) * exp(-x^2/2)
f1 = exp(-(x-1))
f2 = x^(-3/2) / 2

#figure (a)
plot(x, g, type = "l", ylab = "", ylim = c(0,2), lwd = 2, col=1, main='(A)')
lines(x, f1, lty = 2, lwd = 2, col=2)
lines(x, f2, lty = 3, lwd = 2, col=3)
legend("topright", legend = c(expression(g(x)==x^2*e^{-x^2/2}/sqrt(2*pi)),
                              expression(f[1](x)==e^{-(x-1)}),
                              expression(f[2](x)==x^{-3/2} / 2)),
       lty = 1:3, lwd = 2, inset = 0.02, col=1:3, cex = 0.8)

#figure (b)
plot(x, g/f1, type = "l", ylab = "",
     ylim = c(0,5), lwd = 2, lty = 2, col=2, main='(B)')
lines(x, g/f2, lty = 3, lwd = 2, col=3)
legend("topright", legend = c(expression(f[1](x)==e^{-(x-1)}),
                              expression(f[2](x)==x^{-3/2} / 2)),
       lty = 2:6, lwd = 2, inset = 0.02, col=2:6, cex = 0.8)

## ---- echo=FALSE--------------------------------------------------------------
m <- 1e6
set.seed(123)
est <- sd <- numeric(2)
u = runif(m)
g <- function(x) {
  x^2 / sqrt(2*pi) * exp(-x^2/2) * (x > 1)
}
f1 <- function(x) {exp(-(x-1))* (x > 1)}
f2 <- function(x) {x^(-3/2) / 2* (x > 1)}

#using f1
x = 1 - log(1-u)
fg <- g(x) / f1(x)
est[1] <- mean(fg)
sd[1] <- sd(fg)

#using f2
x = (1-u)^(-2)
fg <- g(x) / f2(x)
est[2] <- mean(fg)
sd[2] <- sd(fg)

## -----------------------------------------------------------------------------
m <- 1e6
set.seed(123)
g <- function(x) {
  x^2 / sqrt(2*pi) * exp(-x^2/2) * (x > 1)
}
f1 <- function(x) {exp(-(x-1))* (x > 1)}

u = runif(m)
x = 1 - log(1-u)
est <- mean(g(x)/f1(x))
cat('g(x)的MC估计为：', est)

## -----------------------------------------------------------------------------
m = 1000
theta_j = numeric(5)
theta_SI = numeric(10)
g = function(x){
  exp(-x - log(1+x^2))*(x<1)*(x>0)
}
for(i in 1:10){
  set.seed(i)
  theta_j[1] = mean(g(runif(m/5, 0, 0.2)))
  theta_j[2] = mean(g(runif(m/5, 0.2, 0.4)))
  theta_j[3] = mean(g(runif(m/5, 0.4, 0.6)))
  theta_j[4] = mean(g(runif(m/5, 0.6, 0.8)))
  theta_j[5] = mean(g(runif(m/5, 0.8, 1)))
  theta_SI[i] = mean(theta_j)
}
mu = mean(theta_SI)
sd = sd(theta_SI)

mu
sd

## -----------------------------------------------------------------------------
n = 20
m = 1000
alpha = 0.05

CI = function(n, alpha) {
  x = rchisq(n, 2)
  upCL = mean(x) + sqrt(var(x)/n) * qt(1-alpha/2, n-1)
  lowCL = mean(x) - sqrt(var(x)/n) * qt(1-alpha/2, n-1)
  CL = c(upCL, lowCL)
  return (CL)
}
set.seed(123)
myCI = replicate(m, expr = CI(n, alpha))

prob = (mean(myCI[1,]<2)+mean(myCI[2,]>2))
cat('随机样本的t区间覆盖概率为：', 1-prob)

## -----------------------------------------------------------------------------
m = 1000
alpha = 0.05
mu0 = 1
n = 20
p = matrix(0, 3, m)
p.hat = error = numeric(3)

for (j in 1:m) {
  set.seed(j)
  ttest1 = t.test(rchisq(n, 1), mu = mu0)
  ttest2 = t.test(runif(n, 0, 2), mu = mu0)
  ttest3 = t.test(rexp(n, 1), mu = mu0)
  p[1, j] = ttest1$p.value
  p[2, j] = ttest2$p.value
  p[3, j] = ttest3$p.value
}

for (i in 1:3) {
  p.hat[i] = mean(p[i,])
  error[i] = mean(p[i,]<alpha)
}
p.hat
error

## -----------------------------------------------------------------------------
M = 1000
m = 1000
alpha = 0.1

p_Bonf = p_BH = matrix(NA, M, m)
value_Bonf = value_BH = matrix(NA, M, 3)

for (i in 1:M) {
  set.seed(i)
  p_H0 = runif(0.95*m, 0, 1)
  p_H1 = rbeta(0.05*m, shape1 = 0.1, shape2 = 1)
  p = c(p_H0, p_H1)
  p_Bonf[i,] = p.adjust(p, method = 'bonferroni')
  p_BH[i,] = p.adjust(p, method = 'BH')
  # FWER
  value_Bonf[i, 1] = any(p_Bonf[i,1:950] < alpha)
  value_BH[i, 1] = any(p_BH[i,1:950] < alpha)
  # FDR
  value_Bonf[i, 2] = sum(p_Bonf[i, 1:950] < alpha) / 
    (sum(p_Bonf[i, 1:950] < alpha) + sum(p_Bonf[i, 951:1000] < alpha))
  value_BH[i, 2] = sum(p_BH[i, 1:950] < alpha) / 
    (sum(p_BH[i, 1:950] < alpha) + sum(p_BH[i, 951:1000] < alpha))
  # TPR
  value_Bonf[i, 3] = sum(p_Bonf[i, 951:1000] < alpha) / 
    (sum(p_Bonf[i, 951:1000] < alpha) + sum(p_Bonf[i, 951:1000] >= alpha))
  value_BH[i, 3] = sum(p_BH[i, 951:1000] < alpha) / 
    (sum(p_BH[i, 951:1000] < alpha) + sum(p_BH[i, 951:1000] >= alpha))
}

## -----------------------------------------------------------------------------
value = rbind(colMeans(value_Bonf), colMeans(value_BH))
colnames(value) = c('FWER', 'FDR', 'TPR')
rownames(value) = c('Bonferroni', 'B-H')
knitr::kable(round(value, 3), format = "latex",align='c')

## -----------------------------------------------------------------------------
lambda_bootstrap = function(n){
  # 理论值
  bias_the = 2*n/(n-1) - 2
  sd_the = 2*n/((n-1)*sqrt(n-2))
  
  m = 1000
  B = 1000
  bias = sd = numeric(m)
  for (i in 1:m) {
    x = rexp(n, 2)
    set.seed(i)
    lambda_star = numeric(B)
    # 采样值
    for (b in 1:B) {
      xstar = sample(x, replace = TRUE)
      lambda_star[b] = 1 / mean(xstar)
    }
    bias[i] = mean(lambda_star) - 1 / mean(x)
    sd[i] = sd(lambda_star)
  }
  value = c(mean(bias), bias_the, mean(sd), sd_the)
  return(value)
}

## -----------------------------------------------------------------------------
value = rbind(lambda_bootstrap(5),lambda_bootstrap(10),lambda_bootstrap(20))
colnames(value) = c('Bias-bootstrap', 'Bias-theory', 'Sd-bootstrap', 'Sd-theory')
rownames(value) = c('n=5', 'n=10', 'n=20')
knitr::kable(round(value, 3), format = "latex",align='c')

## ----echo = FALSE-------------------------------------------------------------
library(bootstrap)

## -----------------------------------------------------------------------------
correlation = function(data){
  cor(data[,1],data[,2])
}

## -----------------------------------------------------------------------------
Boot.t.CI = function(x, B = 500, R = 100, level = 0.95, statistic){
  # 初始化
  x = as.matrix(x)
  n = nrow(x)
  stat = numeric(B)
  se = numeric(B)
  alpha = 1 - level
  
  # estimate of standard error
  boot.se <- function(x, R, f) {
    # 初始化
    x = as.matrix(x)
    m = nrow(x)
    # 计算hat theta并返回sd
    theta = replicate(R, expr = {
      i = sample(1:m, size = m, replace = TRUE) 
      f(x[i, ])
      })
    return(sd(theta))
  }
  
  # B times bootstrap
  for (b in 1:B) {
    j = sample(1:n, size = n, replace = TRUE)
    y = x[j, ]
    stat[b] = statistic(y)
    se[b] = boot.se(y, R = R, f = statistic)
  }
  
  # CI
  stat0 = statistic(x)
  t = (stat - stat0) / se
  se0 = sd(stat)
  tstar = quantile(t, c(alpha/2, 1-alpha/2), type = 1)
  names(tstar) = rev(names(tstar))
  CI = rev(stat0 - tstar * se0)
}

## -----------------------------------------------------------------------------
set.seed(5)
data = cbind(law$LSAT,law$GPA)
corCI = Boot.t.CI(data, statistic = correlation, B=2000, R=200)
cat('correlation statistic的95%t置信区间是：', corCI)

## -----------------------------------------------------------------------------
# 导入数据
library(boot)
data = aircondit
lam.hat = mean(data$hours)
# 编写函数
lambda = function(dat, ind) {
  mean(dat[ind, 1]) }

## -----------------------------------------------------------------------------
set.seed(5)
boot_result = boot(data, statistic = lambda, R=2000)
boot_CI = boot.ci(boot_result, type=c("norm", "basic", "perc", "bca"))
print(boot_CI)

## -----------------------------------------------------------------------------
# global variables
lam0 = mean(data$hours)
conf = 0.05
alpha = c(1 - conf/2, conf/2)
zalpha = qnorm(alpha)
B = 2000
n = nrow(data)

lam.normal = lam.basic = lam.bca = numeric(B)
for (b in 1:B) {
  i = sample(1:n, size = n, replace = TRUE)
  dat = data$hours[i]
  lam.normal[b] = lam.basic[b] = lam.bca[b] = mean(dat)
}

## -----------------------------------------------------------------------------
# normal
se.normal = sd(lam.normal)
CI.normal = lam0 - zalpha*se.normal
cat('95% standard normal interval = (', CI.normal[1], ',', CI.normal[2], ')')

## -----------------------------------------------------------------------------
# basic
lam.quant = quantile(lam.basic, alpha)
CI.basic = 2*lam0 - lam.quant
cat('95% basic interval = (', CI.basic[1], ',', CI.basic[2], ')')

## -----------------------------------------------------------------------------
# percentile
cat('95% percentile interval = (', lam.quant[2], ',', lam.quant[1], ')')

## -----------------------------------------------------------------------------
# bca
BCa = function(dat, lam0, lam, stat, conf){
  dat = as.matrix(dat)
  n = nrow(dat)
  N = 1:n
  alpha = c(conf/2, 1 - conf/2)
  zalpha = qnorm(alpha)
  z0 = qnorm(sum(lam < lam0) / length(lam)) # length(lam) = B
  
  lam.jack = numeric(n)
  for (i in 1:n) {
    J = N[1:(n-1)] # 为了和dat[-i, ]的index对应
    jack.data = as.data.frame(dat[-i,1])
    lam.jack[i] <- stat(jack.data, J)
  }
  L = mean(lam.jack) - lam.jack
  a = sum(L^3)/(6 * sum(L^2)^1.5)
  hat.alpha = pnorm(z0 + (z0+zalpha)/(1-a*(z0+zalpha)))
  
  limits = quantile(lam, hat.alpha, type=6) # the formula below 7.9
  return("BCa"=limits)
}

BCa(data, lam0 = lam0, lam = lam.bca, stat = lambda, conf = conf)

## -----------------------------------------------------------------------------
library(bootstrap)
data = scor
n = nrow(data)
cov.fun = function(data, ind){cov(data[ind, ], method = "pearson")}
theta.fun = function(lambda){lambda[1]/sum(lambda)}

## -----------------------------------------------------------------------------
jack_th = numeric(n)
orig_sigma = cov.fun(data)
orig_th = theta.fun(eigen(orig_sigma)$values)
for (i in 1:n) {
    jack_sigma = cov.fun(data, -i)
    jack_th[i] = theta.fun(eigen(jack_sigma)$values)
}
bias = (n - 1) * (mean(jack_th) - orig_th)
se <- sqrt((n-1) * mean((jack_th - mean(jack_th))^2))
cat('jackknife estimate of bias is:', bias)
cat('jackknife estimate of standard error is:', se)

## -----------------------------------------------------------------------------
library(DAAG)
magnetic = ironslag$magnetic
chemical = ironslag$chemical
n = length(magnetic)
e1 = e2 = e3 = e4 = matrix(NA, choose(n, 2), 2)

for (i in 1:n-1) {
  for (j in (i+1):n) {
    k = (i-1)*(n-i/2) + (j-i)
    y <- magnetic[-c(i,j)]
    x <- chemical[-c(i,j)]
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[c(i,j)]
    e1[k, ] <- magnetic[c(i,j)] - yhat1
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[c(i,j)] +
      J2$coef[3] * chemical[c(i,j)]^2
    e2[k, ] <- magnetic[c(i,j)] - yhat2
    
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(i,j)]
    yhat3 <- exp(logyhat3)
    e3[k, ] <- magnetic[c(i,j)] - yhat3
    
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(i,j)])
    yhat4 <- exp(logyhat4)
    e4[k, ] <- magnetic[c(i,j)] - yhat4
  }
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## -----------------------------------------------------------------------------
L2 <- lm(magnetic ~ chemical + I(chemical^2))
L2

## -----------------------------------------------------------------------------
# 导入数据
attach(chickwts)
x = as.vector(weight[feed=='soybean'])
y = as.vector(weight[feed=='linseed'])
detach(chickwts)

## -----------------------------------------------------------------------------
R = 999
z = c(x, y)
n = length(x)
m = length(y)
# 计算每一个元素的经验累计密度函数值，前n个元素为x，后m个元素为y
ecdf_x = rank(x)/n - (rank(z)[1:14]-rank(x))/m
ecdf_y = (rank(z)[15:26]-rank(y))/n - rank(y)/m
K = 1:(m+n)
W_rep = numeric(R)
W0 = m*n*((sum(ecdf_x^2)+sum(ecdf_y^2)))/(m+n)^2
set.seed(5)
for (i in 1:R) {
  k = sample(K, size = 14, replace = FALSE)
  x1 = z[k]
  y1 = z[-k]
  z1 = c(x1, y1)
  ecdf_x1 = rank(x1)/n - (rank(z1)[1:14]-rank(x1))/m
  ecdf_y1 = (rank(z1)[15:26]-rank(y1))/n - rank(y1)/m
  W_rep[i] = m*n*((sum(ecdf_x1^2)+sum(ecdf_y1^2)))/(m+n)^2
}
p = mean(c(W0,W_rep)>=W0)
p

## ----echo=FALSE---------------------------------------------------------------
hist(W_rep, main = 'Permutation distribution of replicates',
     freq = FALSE, xlab = 'W2(p = 0.391)', breaks = 'scott')
points(W0, 0, cex = 1, pch = 16)

## -----------------------------------------------------------------------------
count5test = function(x, y){
  x = x - mean(x)
  y = y - mean(y)
  # count
  extrem_x = sum(x > max(y)) + sum(x < min(y))
  extrem_y = sum(y > max(x)) + sum(y < min(x))
  return(max(c(extrem_x, extrem_y)))
}

## -----------------------------------------------------------------------------
R = 999
n1 = 20
mu1 = mu2 = 0
sigma1 = sigma2 = 1
y_seq = seq(20,100,10)
p1 = numeric(length(y_seq))

for (j in 1:length(y_seq)) {
  set.seed(6)
  n2 = y_seq[j]
  K = 1:(n1+n2)
  x = rnorm(n1, mu1, sigma1)
  y = rnorm(n2, mu2, sigma2)
  C0 = count5test(x, y)
  
  C_rep = numeric(R)
  for (i in 1:R) {
    k = sample(K, size = n1, replace = FALSE)
    z = c(x, y)
    x1 = z[k]
    y1 = z[-k]
    C_rep[i] = count5test(x1, y1)
  }
  p1[j] = mean(c(C0, C_rep) >= C0)
}
p1

## -----------------------------------------------------------------------------
R = 999
n1 = 20
mu1 = mu2 = 0
sigma1 = 1
sigma2 = 1.5
y_seq = seq(20,100,10)
p2 = numeric(length(y_seq))

for (j in 1:length(y_seq)) {
  set.seed(6)
  n2 = y_seq[j]
  K = 1:(n1+n2)
  x = rnorm(n1, mu1, sigma1)
  y = rnorm(n2, mu2, sigma2)
  C0 = count5test(x, y)
  
  C_rep = numeric(R)
  for (i in 1:R) {
    k = sample(K, size = n1, replace = FALSE)
    z = c(x, y)
    x1 = z[k]
    y1 = z[-k]
    C_rep[i] = count5test(x1, y1)
  }
  p2[j] = mean(c(C0, C_rep) >= C0)
}
p2

## -----------------------------------------------------------------------------
value = rbind(p1, p2)
colnames(value) = y_seq
rownames(value) = c('方差相等', '方差不相等')
knitr::kable(round(value, 3), format = "latex", align='c')

## -----------------------------------------------------------------------------
alpha.solve = function(alpha, N, b1, b2, b3, f0){
  x1 = rpois(N, 1)
  x2 = rexp(N, 1)
  x3 = rbinom(N, 1, 0.5)
  g.alpha = function(alpha){
  temp = exp(-alpha-b1*x1-b2*x2-b3*x3)
  mean(1/(1+temp)) - f0
  }
  solution = uniroot(g.alpha, c(-20, 20))
  alpha = round(solution$root, 5)
  return(alpha)
}

set.seed(7)
f = c(0.1, 0.01, 0.001, 0.0001)
alpha.root = numeric(length(f))
for (i in 1:length(f)) {
  alpha.root[i] = alpha.solve(N=1e6, b1=0, b2=1, b3=-1, f0=f[i])
}
alpha.root

## ---- echo=FALSE--------------------------------------------------------------
# plot
f1 = seq(0.0001, 0.1, 0.001)
alpha.root1 = numeric(length(f1))
for (i in 1:length(f1)) {
  alpha.root1[i] = alpha.solve(N=1e6, b1=0, b2=1, b3=-1, f0=f1[i])
}
log.f0 = -log(f1)
plot(log.f0, alpha.root1)

## -----------------------------------------------------------------------------
# Metropolis函数，返回生成的样本x和接受率rate
laplace = function(x){exp(-abs(x))/2}
Metropolis = function(N, x0, sigma) {
  x = numeric(N)
  u = runif(N)
  x[1] = x0
  rate = 0
  for (i in 2:N) {
    y = rnorm(1, x[i-1], sigma)
    if(u[i] <= laplace(y) / laplace(x[i-1])){
      x[i] = y
      rate = rate + 1
    }
    else{
      x[i] = x[i-1]
    }
  }
  return(list(x = x, rate = rate/N))
}

## -----------------------------------------------------------------------------
set.seed(7)
N = 5000
sigma = c(0.05, 0.5, 2,  16)
x0 = 30
rw1 = Metropolis(N, x0, sigma[1])
rw2 = Metropolis(N, x0, sigma[2])
rw3 = Metropolis(N, x0, sigma[3])
rw4 = Metropolis(N, x0, sigma[4])

# 输出接受率
accepted = data.frame(sigma = sigma, 
                      accepted.rate = c(rw1$rate, rw2$rate, rw3$rate, rw4$rate))
knitr::kable(accepted,format='latex')

## ---- echo=FALSE--------------------------------------------------------------
# par(mfrow=c(2,2))
quantile = c(-3, 3)
rw = cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
for (j in 1:4) {
  plot(rw[,j], type="l",
       xlab = bquote(sigma == .(round(sigma[j],3))),
       ylab = "X", ylim=range(rw[,j]))
  abline(h = quantile, col = 'red')
}

## -----------------------------------------------------------------------------
# 定义参数
N = 5000
burn = 1000
S = matrix(0, N, 2)
rho = 0.9
mu1 = mu2 = 0
sigma1 = sigma2 = 1
k1 = sqrt(1-rho^2)*sigma1
k2 = sqrt(1-rho^2)*sigma2
# 初始化
S[1, ] = c(mu1, mu2)
for (i in 2:N) {
  y = S[i-1, 2]
  m1 = mu1 + rho * (y - mu2) * sigma1/sigma2
  S[i, 1] = rnorm(1, m1, k1)
  x = S[i, 1]
  m2 = mu2 + rho * (x - mu1) * sigma2/sigma1
  S[i, 2] = rnorm(1, m2, k2)
}
s = S[(burn + 1):N, ]

## ---- echo=FALSE--------------------------------------------------------------
plot(s[,1], type='l', col=1, lwd=2, xlab='Index', ylab='Random numbers')
lines(s[,2], col=2, lwd=2)
legend('topright', c(expression(X_t), expression(Y_t)), col=1:2, lwd=2)

## ---- echo=FALSE--------------------------------------------------------------
X = s[,1]
Y = s[,2]
model = lm(Y ~ X)
residuals <- resid(model)
# 残差图判断方差齐性
plot(residuals, main = "Residual Plot", xlab = "Observation", ylab = "Residuals")
# QQ图判断正态性
qqnorm(residuals)
qqline(residuals)

## -----------------------------------------------------------------------------
# Rayleigh function returns pdf
Rayleigh = function(x, sigma) {
  if (any(x < 0)) return (0)
  stopifnot(sigma > 0)
  return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}
# Rchain function generates chain for Rayleigh dist
Rchain = function(sigma, N, X1) {
  x = rep(0, N)
  x[1] = X1
  u = runif(N)
  for (i in 2:N) {
    xt = x[i-1]
    y = rchisq(1, df = xt)
    r1 = Rayleigh(y, sigma) * dchisq(xt, df = y)
    r2 = Rayleigh(xt, sigma) * dchisq(y, df = xt)
    if (u[i] <= r1/r2){
      x[i] = y
    }else{
      x[i] = xt
    }
  }
  return(x)
}

## -----------------------------------------------------------------------------
# 函数Gelman.Rubin以psi作为输入
Gelman.Rubin = function(psi) {
  psi = as.matrix(psi)
  n = ncol(psi)
  k = nrow(psi)
  B = n * var(rowMeans(psi))
  W = mean(apply(psi, 1, "var"))
  v.hat = W*(n-1)/n + B/n
  r.hat = v.hat / W
  return(r.hat)
}

## -----------------------------------------------------------------------------
# 参数
set.seed(7)
sigma = 4
burn = 1000
chain.len = 20000
chain.num = 4
# x0 = rchisq(4, df=1)
x0 = c(2, 6, 10, 14)   # initial values
# 生成四条不同初始值的链
MCMC.chain = matrix(0, nrow = chain.num, ncol = chain.len)
for (i in 1:chain.num){
  MCMC.chain[i, ] = Rchain(sigma, chain.len, x0[i])
}

psi0 = t(apply(MCMC.chain, 1, cumsum))
psi = psi0
for (i in 1:nrow(psi0)){
  psi[i,] = psi0[i,] / (1:ncol(psi0))
}

## ---- echo=FALSE--------------------------------------------------------------
# plot
for (i in 1:chain.num){
  if(i==1){
    plot((burn+1):chain.len, psi[i, (burn+1):chain.len], ylim=c(4,6), type="l",
         xlab='Index', ylab=bquote(phi))
  }else{
    lines(psi[i, (burn+1):chain.len], col=i)
  }
}

rhat = rep(0, chain.len)
for (j in (burn+1):chain.len){
  rhat[j] = Gelman.Rubin(psi[,1:j])
}
plot(rhat[(burn+1):chain.len], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## ----echo=FALSE---------------------------------------------------------------
library(coda)

## -----------------------------------------------------------------------------
obj = mcmc.list(mcmc(MCMC.chain[1, ]), mcmc(MCMC.chain[2, ]), 
                 mcmc(MCMC.chain[3, ]), mcmc(MCMC.chain[4, ]))
gelman.diag(obj)
gelman.plot(obj)

## -----------------------------------------------------------------------------
ui = c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
vi = c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3)
partial.l = function(lambda){
  equation = (-ui*exp(-lambda*(ui-vi))+vi) / (exp(-lambda*(ui-vi))-1)
  return(sum(equation))
}
cat("Lambda estimated by MLE:", uniroot(partial.l, lower = 0.01, upper = 3)$root)

## -----------------------------------------------------------------------------
# E-Step 的条件期望的被积函数
integrand = function(x, lambda, uk, vk) {
  numerator = x * lambda * exp(-lambda * x)
  denominator = exp(-lambda * uk) - exp(-lambda * vk)
  return(numerator / denominator)
}

# EM算法
em_algorithm = function(lambda_initial, num_iterations, eps, u, v) {
  n = length(u)
  lambda = lambda_initial
  pre_lambda = lambda_initial
  # 迭代
  for (iteration in 1:num_iterations) {
    E = numeric(n)
    for (k in 1:n) {
      E[k] = integrate(integrand, lower = u[k], upper = v[k], 
                       lambda = lambda, uk = u[k], vk = v[k])$value
    }
    lambda = n/sum(E)
    # whether to stop
    if (abs(lambda - pre_lambda) < eps) {
      break
    }
    pre_lambda = lambda
  }
  return(lambda)
}

lambda_initial = 0.01
num_iterations = 1e4
eps = 1e-5
em_lambda = em_algorithm(lambda_initial, num_iterations, eps, u = ui, v = vi)
cat("Lambda estimated by EM algorithm:", em_lambda)

## -----------------------------------------------------------------------------
library(boot)
solve.game = function(payoff.matrix){
  # preprocess
  min.A = min(A); A = A - min.A; max.A = max(A); A = A / max(A)
  m = nrow(A)
  n = ncol(A)
  iter = n^3
  
  # max v
  # objective function
  a = c(rep(0, m), 1)
  # constraint 1 >=，但simplex函数中的A1是<=，所以这里A1前有负号
  A1 = -cbind(t(A), rep(-1, n)) 
  b1 = rep(0, n)
  # constraint 2 sum(x)=1
  A3 = t(as.matrix(c(rep(1, m), 0)))
  b3 = 1
  sx = simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3, maxi=TRUE, n.iter=iter)
  
  # min v
  # objective function
  a = c(rep(0, n), 1)
  # constraint 1 <=
  A1 = cbind(A, rep(-1, m))
  b1 = rep(0, m)
  # constraint 2 sum(y)=1
  A3 = t(as.matrix(c(rep(1, n), 0)))
  b3 = 1
  sy = simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=FALSE, n.iter=iter)
  solution.list = list("A" = A * max.A + min.A,
                "x" = sx$soln[1:m],
                "y" = sy$soln[1:n],
                "v" = sx$soln[m+1] * max.A + min.A)
  solution.list
}

## -----------------------------------------------------------------------------
# payoff matrix
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)
B <- A + 2
# 计算A的解和B的解
round(cbind(solve.game(A)$x, solve.game(A)$y),7)
round(cbind(solve.game(B)$x, solve.game(B)$y),7)
# 验证game B的解是否为game A的极值点之一
epointA = round(c(0, 0, 25/61, 0, 20/61, 0, 16/61, 0, 0), 7)
epointA == round(solve.game(B)$x, 7)
epointA == round(solve.game(B)$y, 7)

## -----------------------------------------------------------------------------
solve.game(A)$v
solve.game(B)$v 

## -----------------------------------------------------------------------------
list1 = list('a', 2, 5, 'd')
unlist(list1)
as.vector(list1)

## -----------------------------------------------------------------------------
v = list(1, 2, 5, 7)
dim(v)
length(v)

## -----------------------------------------------------------------------------
m = matrix(c(1, 2, 5, 7), 2, 2)
is.matrix(m)
is.array(m)

## -----------------------------------------------------------------------------
df <- data.frame(
  col1 = c("a", "b", "c"),
  col2 = c(1, 2, 8),
  col3 = c(TRUE, TRUE, FALSE),
  stringsAsFactors = FALSE
)
cat('type of column 1:', typeof(df$col1), '\n')
cat('type of column 2:', typeof(df$col2), '\n')
cat('type of column 3:', typeof(df$col3), '\n')

typeof(as.matrix(df))

## -----------------------------------------------------------------------------
df = data.frame()
ncol(df)
nrow(df)

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
df1 = data.frame(
  col1 = c(1, 2, NA, 4),
  col2 = c(5, NA, 7, 8),
  col3 = c(NA, 10, 11, 12)
)

lapply(df1, scale01)

## -----------------------------------------------------------------------------
df2 = data.frame(
  col1 = seq(1,4,1),
  col2 = rep(TRUE, 4),
  col3 = seq(9,15,2),
  col4 = rep('a', 4),
  stringsAsFactors = FALSE
)
numeric.idx = vapply(df2, is.numeric, FUN.VALUE=logical(1))
lapply(df2[, numeric.idx], scale01)

## -----------------------------------------------------------------------------
df1 = data.frame(
  col1 = seq(1,4,1),
  col2 = seq(5,8,1),
  col3 = seq(9,15,2)
)
vapply(df1, sd, FUN.VALUE=0)

## -----------------------------------------------------------------------------
df2 = data.frame(
  col1 = seq(1,4,1),
  col2 = rep(TRUE, 4),
  col3 = seq(9,15,2),
  col4 = rep('a', 4),
  stringsAsFactors = FALSE
)
numeric.idx = vapply(df2, is.numeric, FUN.VALUE=logical(1))
vapply(df2[, numeric.idx], sd, FUN.VALUE=0)

## -----------------------------------------------------------------------------
gibbsR = function(N, burn, a, b, n, seed){
  # N: length of chain
  # burn: burn-in length
  # a, b, n: parameter of density
  # seed: random seed
  set.seed(seed)
  data = matrix(0, N, 2)
  data[1, ] = c(n-1, 1/2)
  for (i in 2:N) {
    y = data[i-1, 2]
    data[i, 1] = rbinom(1, n, y)
    x = data[i, 1]
    data[i, 2] = rbeta(1, x+a, n-x+b)
  }
  chain = data[(burn+1):N,]
  return(chain)
}

## -----------------------------------------------------------------------------
library(microbenchmark)
library(Rcpp)
time = microbenchmark(timeR=gibbsR(5000, 1000, 2, 3, 9, 12345), 
                      timeCpp = gibbsCpp(5000, 1000, 2, 3, 9))
summary(time)[,c(1,3,5,6)]

