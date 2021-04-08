# Precision of the estimates

rm(list=ls())
library(PLNmodels); library(mvtnorm)
par(pch=20)

# General parms
xtol = 1e-4
ftol  = 1e-6
lbvar = 1e-7
B = 1e2;
H = 1e3
n.list = c(100, 500, 30)
p.list = c(2, 3, 5)
d.list = c(1, 2, 3)

# for (n in n.list){
   # for (dp in 1:length(p.list)){
      # p = p.list[dp]; d = d.list[dp]
      n = 1e2; p = 3; d = 2
      cat('\nn=', n, ' p=', p, ' d=', d, '\n')

      # Simul parms
      parms.file = paste('SimulParms_n', n, '_p', p, '_d', d, '.Rdata', sep='')
      if(file.exists(parms.file)){load(parms.file)}else{
         beta = matrix((1:(d*p))-((d*p+1)/2), d, p) / sqrt(d*p)
         o = 1; O = matrix(o, n, p)
         X = matrix(rnorm(n*d), n, d)
         Sigma = exp(-as.matrix(dist(matrix(rnorm(2*p), p, 2))))
         parms = list(n=n, p=p, d=d, X=X, beta=beta, Sigma=Sigma, o=o, O=O)
         save(parms, file=parms.file)
      }
      X = parms$X; beta = parms$beta; Sigma = parms$Sigma; O = parms$O; o = parms$o
      Omega = solve(Sigma); oXbeta = o+X%*%beta
      Sigma.vec = as.vector(Sigma[upper.tri(Sigma, diag=T)])
      Omega.vec = as.vector(Omega[upper.tri(Omega, diag=T)])

      # Simul
      beta.pln = beta.glm = beta.pcaopt = beta.plnglm = beta.pcamax = array(dim=c(B, d, p));
      P = p*(p+1)/2; Sigma.pln= Sigma.pcaopt= Sigma.pcamax = Omega.pln = Omega.pcamax = matrix(NA, B, P)
      rank.pcaopt = rep(NA, B)
      par(mfrow=c(5, 5), mex=.5, pch=20)
      for (b in 1:B){
         cat(b, '')
         # Simul
         Z = rmvnorm(n, rep(0, p), Sigma)
         Y = matrix(rpois(n*p, exp(oXbeta + Z)), n, p)
         # Fit GLM
         invisible(sapply(1:p, function(j){
            GLM = glm(Y[, j] ~ -1 + X + offset(O[, j]), family='poisson')
            beta.glm[b, , j] <<- GLM$coef
         }))
         # Fit PLN
         PLN = PLN(Y, X=X, O=O, control=list(trace=0, lbvar=0.1,  xtol=xtol, method="LBFGS"))
         beta.pln[b, , ] = t(PLN$model.par$Theta)
         Sigma.pln[b, ] = as.vector(PLN$model.par$Sigma[upper.tri(PLN$model.par$Sigma, diag=T)])
         Omega.pln[b, ] = as.vector(PLN$model.par$Omega[upper.tri(PLN$model.par$Omega, diag=T)])
         # Fit of PLN
         logApln = (O + tcrossprod(X,PLN$model.par$Theta) + PLN$variational.par$M + PLN$variational.par$S/2)
         plot(log(1+Y), logApln, xlab='', ylab='', main=paste(n, p, d, b)); abline(a=0, b=1)
         # # Fit GLM on PLN residuals
         # Aoffset = exp(O + PLN$variational.par$M + PLN$variational.par$S/2)
         # invisible(sapply(1:p, function(j){
         #    GLM = glm(Y[, j] ~ -1 + X + offset(Aoffset[, j]), family='poisson')
         #    beta.plnglm[b, , j] <<- GLM$coef
         # }))
         # Fit PLNpcamax
         PCA = PLNPCA(Y, X=X, O=O, ranks=p, control=list(trace=0, xtol=xtol, ftol=ftol, lbvar=lbvar))
         PCAmax = PCA$getModel(p)
         beta.pcamax[b, , ] = t(PCAmax$model.par$Theta)
         Sigma.pcamax[b, ] = as.vector(PCAmax$model.par$Sigma[upper.tri(PCAmax$model.par$Sigma, diag=T)])
         if(min(eigen(PCAmax$model.par$Sigma)$values) > 1e-10){
            Omega.tmp = solve(PCAmax$model.par$Sigma)
            Omega.pcamax[b, ] = as.vector(Omega.tmp[upper.tri(Omega.tmp, diag=T)])
         }
         # # Fit PLNpca
         # PCA = PLNPCA(Y, X=X, O=O, ranks=1:p, control=list(trace=0, xtol=xtol, ftol=ftol))
         # PCA = PCA$getBestModel()
         # rank.pcaopt[b] = PCA$rank
         # beta.pcaopt[b, , ] = t(PCA$model.par$Theta)
         # Sigma.pcaopt[b, ] = as.vector(PCA$model.par$Sigma[upper.tri(PCA$model.par$Sigma, diag=T)])
         # Fit of PCAmax
         logAmax = (O + X%*%t(PCAmax$model.par$Theta) + PCAmax$variational.par$M + PCAmax$variational.par$S/2)
         points(log(1+Y), logAmax, col=4, pch=21)
      }

      # # Reshape beta.pln
      # beta.plnres = array(dim=c(B, d, p))
      # invisible(sapply(1:B, function(b){beta.plnres[b, , ] <<- t(matrix(as.vector(beta.pln[b, , ]), p, d))}))

      # Results for beta
      pdf(paste('PLNsim-beta_n', n, '_p', p, '_d', d, '_xtol', xtol, '.pdf', sep=''))
      par(mfcol=c(3, d*p), mex=.75); beta.lim = c(-1, 1)
      invisible(sapply(1:d, function(k){sapply(1:p, function(j){
         hist(beta.glm[, k, j]-beta[k, j], xlab='', ylab='', main='beta: glm', breaks=sqrt(B), xlim=beta.lim);
         abline(v=0, lwd=2, col=2)
         hist(beta.pln[, k, j]-beta[k, j], xlab='', ylab='', main='beta: pln', breaks=sqrt(B), xlim=beta.lim);
         abline(v=0, lwd=2, col=2)
         # hist(beta.plnres[, k, j]-beta[k, j], xlab='', ylab='', main='beta: pln-r', breaks=sqrt(B), xlim=beta.lim);
         # abline(v=0, lwd=2, col=2)
         # hist(beta.pln2[, k, j]-beta[k, j], xlab='', ylab='', main='', breaks=sqrt(B), xlim=beta.lim); abline(v=0, lwd=2, col=2)
         # hist(beta.plnglm[, k, j]-beta[k, j], xlab='', ylab='', main='', breaks=sqrt(B), xlim=beta.lim); abline(v=0, lwd=2, col=2)
         # hist(beta.pcaopt[, k, j]-beta[k, j], xlab='', ylab='', main='beta: pca-o', breaks=sqrt(B), xlim=beta.lim);
         # abline(v=0, lwd=2, col=2)
         hist(beta.pcamax[, k, j]-beta[k, j], xlab='', ylab='', main='beta: pca-m', breaks=sqrt(B), xlim=beta.lim);
         abline(v=0, lwd=2, col=2)
      })}))
      dev.off()

      # Results for Sigma
      pdf(paste('PLNsim-Sigma_n', n, '_p', p, '_d', d, '_xtol', xtol, '.pdf', sep=''))
      par(mfcol=c(4, P), mex=.75); Sigma.lim = c(-1, 1)
      invisible(sapply(1:P, function(p){
         hist(Sigma.pln[, p]-Sigma.vec[p], xlab='', ylab='', main='Sigma: pln', breaks=sqrt(B), xlim=Sigma.lim);
         abline(v=0, lwd=2, col=2)
         # hist(Sigma.pcaopt[, p]-Sigma.vec[p], xlab='', ylab='', main='Sigma: pca-o', breaks=sqrt(B), xlim=Sigma.lim);
         # abline(v=0, lwd=2, col=2)
         hist(Sigma.pcamax[, p]-Sigma.vec[p], xlab='', ylab='', main='Sigma: pca-m', breaks=sqrt(B), xlim=Sigma.lim);
         abline(v=0, lwd=2, col=2)
         hist(Omega.pln[, p]-Omega.vec[p], xlab='', ylab='', main='Omega: pln', breaks=sqrt(B));
         abline(v=0, lwd=2, col=2)
         hist(Omega.pcamax[, p]-Omega.vec[p], xlab='', ylab='', main='Omega: pca-m', breaks=sqrt(B));
         abline(v=0, lwd=2, col=2)
      }))
      dev.off()

      summary(rank.pcaopt)
   # }
# }
