library("MCMCglmm")

ped<-read.table("ped_loop.txt",header=T,na.strings="NA")

library("nadiv")
library("INLA")

reco=dim(ped)[1]

tr=reco*3

A=makeA(ped)


var_add =data.frame(add1=1,add2=1,add3=1)

var_err =data.frame(err1=1,err2=1,err3=1)

cov_add=data.frame(cov_add12=1,cov_add13=1,cov_add23=1)

cov_err=data.frame(cov_err12=1,cov_err13=1,cov_err23=1)

bv=data.frame(breed=1)
bv_a=data.frame(breed=1)



A_ch=chol(A)

A_ch=t(A_ch)


#G=c(1,0.4,0.6,0.4,1.5,0.9,0.6,0.9,2.5) # these are the simulated additive variance and covariance for three traits

G=c(5,-2,3,-2,7,4,3,4,10)#high heritability (0.7)
G=matrix(G,nrow=3)


R=c(20,-5,1,-5,28,3,1,3,35)#low heritabilty 0.2
R=matrix(R,nrow=3)

 G_ch=chol(G)
G_ch=t(G_ch)

R_ch=chol(R)
R_ch=t(R_ch)

K_A=kronecker(G_ch,A_ch)

K_E=kronecker(R_ch,diag(reco))

 




a=K_A%*%rnorm(tr,0,1)


e=K_E%*%rnorm(tr,0,1)



y=a+e
Tra1=y[1:reco] # first trait
Tra2=y[(reco+1):(2*reco)]# second trait
Tra3=y[(2*reco+1):(3*reco)]# third trait



dat=data.frame(Line=1:reco,Tra1=Tra1,Tra2=Tra2,Tra3=Tra3)

dat$Line=as.factor(dat$Line)


A=inverseA(ped)


n=reco
gaussian.1=c(dat$Tra1,rep(NA,reco),rep(NA,reco))
gaussian.2=c(rep(NA,reco),dat$Tra2,rep(NA,reco))
gaussian.3=c(rep(NA,reco),rep(NA,reco),dat$Tra3)
joint.response=list(gaussian.1, gaussian.2, gaussian.3)
##examples of priors settings
prec.A = list(param = c(0.5, 0.5), fixed = FALSE)
prec.rho =list(param = c(0,0.1), fixed =FALSE)
prec.e = list(param = c(0.5, 0.5), fixed = FALSE)

random.covariates=list(
y1a1=c(1:n,rep(NA,n),rep(NA,n)),
y2a2=c(rep(NA,n),1:n,rep(NA,n)),
y3a3=c(rep(NA,n),rep(NA,n),1:n),
y1e1=c(1:n,rep(NA,n),rep(NA,n)),
y2e2=c(rep(NA,n),1:n,rep(NA,n)),
y3e3=c(rep(NA,n),rep(NA,n),1:n),
y2a1rho=c(rep(NA,n),1:n,rep(NA,n)),
y3a1rho=c(rep(NA,n),rep(NA,n),1:n),
y3a2rho=c(rep(NA,n),rep(NA,n),1:n),
y2e1rho=c(rep(NA,n),1:n,rep(NA,n)),
y3e1rho=c(rep(NA,n),rep(NA,n),1:n),
y3e2rho=c(rep(NA,n),rep(NA,n),1:n))

##Make a data list
joint.data= random.covariates
joint.data$yyy=joint.response
formula.tri = yyy~f(y1a1,model = "generic0", hyper= list(theta = prec.A),constr=FALSE, Cmatrix=A$Ainv)+
f(y2a2,model = "generic0", hyper= list(theta = prec.A),constr=FALSE,Cmatrix=A$Ainv)+
f(y2a1rho,copy = "y1a1", hyper =list(theta = prec.rho))+f(y3a3,model = "generic0", hyper= list(theta = prec.A),constr=FALSE, Cmatrix=A$Ainv)+
f(y3a1rho,copy = "y1a1", hyper =list(theta = prec.rho))+f(y3a2rho,copy = "y2a2", hyper =list(theta = prec.rho))+
f(y1e1,model ="iid",hyper = list(theta = prec.e),constr=TRUE)+f(y2e2,model ="iid", hyper = list(theta = prec.e),constr=TRUE)+
f(y2e1rho,copy= "y1e1", hyper = list(theta = prec.rho))+f(y3e3,model ="iid", hyper = list(theta = prec.e),constr=TRUE)+
f(y3e1rho,copy= "y1e1", hyper = list(theta = prec.rho))+f(y3e2rho,copy= "y2e2", hyper = list(theta = prec.rho))-1

model = inla(formula.tri, family = c("gaussian", "gaussian","gaussian"), data=joint.data,verbose=F,
control.compute=list(dic=T),control.family=list(list(hyper = list(prec = list(initial=10 ,fixed = TRUE))),list(hyper = list(prec = list(initial=10 ,fixed = TRUE))),
list(hyper = list(prec = list(initial=10 ,fixed = TRUE)))),only.hyperparam = F)

g=inla.hyperpar.sample(n=10000,model,intern=FALSE)

a1=mean(1/g[,1])
a2=mean(g[,7]*g[,7]*1/g[,1]+1/g[,3])
a3=mean(g[,8]*g[,8]*1/g[,1]+g[,9]*g[,9]*1/g[,3]+1/g[,4])
e1=mean(1/g[,2])
e2=mean(g[,11]*g[,11]*1/g[,2]+1/g[,5])
e3=mean(g[,10]*g[,10]*1/g[,2]+g[,12]*g[,12]*1/g[,5]+1/g[,6])
a12=mean(g[,7]*1/g[,1])
a13=mean(g[,8]*1/g[,1])
a23=mean(g[,7]*g[,8]*1/g[,1]+g[,9]*1/g[,3])
e12=mean(g[,11]*1/g[,2])
e13=mean(g[,10]*1/g[,2])
e23=mean(g[,10]*1/g[,2]*g[,11]+g[,12]*1/g[,5])

var_add=rbind(var_add,c(a1,a2,a3))
var_err=rbind(var_err,c(e1,e2,e3))

cov_add=rbind(cov_add,c(a12,a13,a23))
cov_err=rbind(cov_err,c(e12,e13,e23))


breed=c(model$summary.random$y1a1$mean,model$summary.random$y2a2$mean,model$summary.random$y3a3$mean)

co_a=cor(breed,as.vector(a))
bv_a=rbind(bv,c(co_a))


co=cor(breed,as.vector(y))
bv=rbind(bv,c(co))




