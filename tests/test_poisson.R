set.seed(19970224)
nt = 100
n = 200
d = 100
p = 3
m = d * p
t = 1

library(splines)
library(SAM)

X = 0.5*matrix(runif(n*d),n,d) + matrix(rep(0.5*runif(n),d),n,d)
u = exp(-2*sin(X[,1]) + X[,2]^2-1/3 + X[,3]-1/2 + exp(-X[,4])+exp(-1)-1+1)
y = rep(0,n)
for(i in 1:n) y[i] = rpois(1,u[i])

Xt = 0.5*matrix(runif(nt*d),nt,d) + matrix(rep(0.5*runif(nt),d),nt,d)
ut = exp(-2*sin(Xt[,1]) + Xt[,2]^2-1/3 + Xt[,3]-1/2 + exp(-Xt[,4])+exp(-1)-1+1)
yt = rep(0,nt)
for(i in 1:nt) yt[i] = rpois(1,ut[i])

total_t = 0
total_l = 0
nlamb = 20
for (i in 1:t) {
  t0 = proc.time()
  out.trn = samEL(X, y, nlambda=nlamb)
  total_t = total_t + proc.time() - t0
  out.tst = predict(out.trn, Xt)
  total_l = total_l + mean((out.tst$expectations[,nlamb]-yt)^2)
}
print("sam Poisson regression:")
print(total_t / t)
cat("MSE:", total_l / t, "\n")
stopifnot(is.finite(total_l / t))
