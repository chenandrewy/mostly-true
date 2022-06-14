# exponential test ====

lambda = 2

n = 1e5

x = rexp(n, 1/lambda)


q = 0.2

quantile(x,q)

lambda*log(1/(1-q))


h = 2
y = x[x>h]

quantile(y,q)
lambda*log(1/(1-q)) + h