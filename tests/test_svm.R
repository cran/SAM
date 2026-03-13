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
y = sign(((X[,1]-0.5)^2 + (X[,2]-0.5)^2)-0.06)
## flipping about 5 percent of y
y = y*sign(runif(n)-0.05)

Xt = 0.5*matrix(runif(nt*d),nt,d) + matrix(rep(0.5*runif(nt),d),nt,d)
yt = sign(((Xt[,1]-0.5)^2 + (Xt[,2]-0.5)^2)-0.06)
yt = yt*sign(runif(nt)-0.05)

total_t = 0
total_l = 0
nlamb = 20
for (i in 1:t) {
  t0 = proc.time()
  out.trn = samHL(X, y, nlambda=nlamb)
  total_t = total_t + proc.time() - t0
  out.tst = predict(out.trn, Xt)
  total_l = total_l + mean(out.tst$labels[,nlamb]==yt)
}
print("sam SVM (hinge loss):")
print(total_t / t)
cat("accuracy:", total_l / t, "\n")
stopifnot(total_l / t > 0.5)
