#===============================================================================#

# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#===============================================================================#
# Notes :
# when Eror in FUN(X[[i]],...) : shows, run devtools::load_all(".")
#
# 1. use the current R Development Version (that will eventually become 3.4)
# 2. Run the tools::package_native_routine_registration_skeleton(".") and copy and
#    paste the full output in a packagename_init.c file to be put in src/
# 3. update NAMESPACE, verifying that useDynLib(packagename, .registration = TRUE)
# 4. If necessary, replace the exportPattern with export( list of object to be exported )
#===============================================================================#

#' Linear regression via coordinate descent with covariate clustering
#'
#' Covariate assignment to k clusters. The es
#'
#' @param Y vector of outcome variable
#' @param X matrix of covariates
#' @param k number of clusters
#' @param coefs vector of coefficients as starting values
#' @param clus vector of covariate cluster assignments as starting values
#' @param clusmns vector k cluster parameter centers
#' @param nC first nC-1 covariates in X not to cluster. Must be at least 1 for the intercept
#' @return clus cluster assignments
#' @return coefs vector of coefficients as starting values
#' @return clusmns vector of cluster means
#'
#'
#' @examples
#' set.seed(14) #Generate data
#' N = 1000; (bets = rep(-2:2,4)); p = length(bets); X = matrix(rnorm(N*p),N,p)
#' Y = cbind(1,X)%*%matrix(c(0.5,bets),ncol = 1)
#' begin_v<- rep(NA,p)
#' for (j in 1:p) {
#'  begin_v[j] = coef(lm(Y~X[,j]))[2]
#' }
#' set.seed(12); klus_obj<- kmeans(begin_v,centers = 5)
#' linrclus(Y,X,k=5,coefs=c(0,begin_v),clus=klus_obj$cluster,clusmns=klus_obj$centers)


#' @useDynLib cluscov linreg_coord_clus
#' @export

linrclus<- function(Y,X,k,coefs,clus,clusmns,nC=1){
  X = cbind(1,X)
  nrX = as.integer(nrow(X)); ncX = as.integer(ncol(X)); k = as.integer(k)
  Xdot = apply(X, 2, function(x)sum(x^2)); clus=as.integer(c(clus-1)); nC=as.integer(nC)
  ans=.C("linreg_coord_clus",as.double(Y),as.double(X),coefs=as.double(coefs),clus=as.integer(clus),
         klus=as.integer(c(0,clus)),double(nrX),as.double(Xdot),clusmns=as.double(clusmns),
         nrX,ncX,k,nC,PACKAGE = "cluscov")
  list(clus=(ans$clus+1),coefs=ans$coefs,clusmns=ans$clusmns)
}



#===================================================================================>
#' Integer Golden Search Minimisation
#'
#' This function conducts an integer golden search minimisation of a univariate function.
#'
#' @param fn  function to be minimised
#' @param interval a vector of length two containing the minimum and maximum interger
#' values within which to search for the minimiser.
#' @param  tol the tolerance level. Defaults at 1
#'
#' @return k minimiser of \code{fn()}
#' @return crit the minimum
#' @return iter total number of iterations
#' @return iterfn total number of function evaluations of \code{fn()}
#' @return fobj an object of the function minimisation
#' @return key a logical for warning if \code{fobj} may not correspond to \code{k}
#'
#' @examples
#' ## Not run:...
#' @export

# Golden Search Algorithm for the outer loop ()
goldopt<- function(fn,interval,tol=1){
  a=min(interval); b = max(interval); key=0
  xvals<- c(0); xvals[1]<- a; xvals[2]<- b
  fvals<- c(0)
  faobj<-fn(a); fa = faobj$BIC; fbobj<-fn(b); fb = fbobj$BIC;
  fvals[1]<- fa; fvals[2]<- fb
  cnt=2; # set counter to 2 for function evaluations
  phi<- (1+sqrt(5))/2

  c = ceiling(b - (b-a)/phi); d = floor(a + (b-a)/phi);
  if(any(xvals==c)){
    c=xvals[which(xvals==c)]; fc=fvals[which(xvals==c)]
    warning("Function object may not correspond to optimal number of clusters")
    key=1
  }else{
    cnt=cnt+1
    fcobj<- fn(c); fc=fcobj$BIC
  }
  if(any(xvals==d)){
    d=xvals[which(xvals==d)]; fd=fvals[which(xvals==d)];
    warning("Function object may not correspond to optimal number of clusters")
    key=1
  }else{
    cnt=cnt+1
    fdobj<- fn(d); fd=fdobj$BIC
  }


  l = 1; # set counter for iterations
  cat("iter = ",l,"\n");arreter = 0
  while(abs(c-d)>tol & arreter==0){# while c d
    if(fc<fd){
      b=d; fb=fd; fbobj=fdobj
      d=c;fd=fc; fdobj=fcobj
      c = ceiling(b-(b-a)/phi);
      if(any(xvals==c)){
        c=xvals[which(xvals==c)]; fc=fvals[which(xvals==c)]
        warning("Function object may not correspond to optimal number of clusters")
        key=1
      }else{
        cnt=cnt+1
        fcobj<- fn(c); fc=fcobj$BIC
      }

    }else if(fc>fd){
      a=c; fa = fc; faobj=fcobj;
      c=d; fc=fd; fcobj = fdobj
      d = floor(a + (b-a)/phi);
      if(any(xvals==d)){
        d=xvals[which(xvals==d)]; fd=fvals[which(xvals==d)]
        warning("Function object may not correspond to optimal number of clusters")
        key=1
      }else{
        cnt=cnt+1
        fdobj<- fn(d); fd=fdobj$BIC
      }
    }else{
      arreter=1
    }
    l=l+1; cat("iter = ",l,"\n")
  }

  optx=ifelse(fc>fd,d,c);
  if(fc>fd){
    fobj=fdobj
  }else{
    fobj=fcobj
  }
  res=list(k=optx,crit=min(fd,fc),iter=l,iterfn=cnt,fobj=fobj,key=key)
  return(res)
}
#===================================================================================>
#' Choice of model
#'
#' This function is used to choose th appropriate model
#'
#' @param  Y the dependent or response variable
#' @param  X design matrix (without intercept)
#' @param  model a string for the desire model name. The following are supported
#'           "lm", "logit", "probit", "gammainverse", "gammaidentity", "gammalog",
#'           "poissonlog", "poissonsqrt", "negbin", "quantreg"
#' @param  ... additional model arguments
#'
#' @return mobj fitted model object
#'
#' @examples
#' ## Not run: summary(ch.model(cars$speed,cars$dist,"gammaidentity"))
#' @export

ch.model<- function(Y,X,model="lm",...){
  if(model=="lm"){
    mobj<- stats::glm(Y~.,family = gaussian(link = "identity"),data = data.frame(Y,X),...)
  }else if(model=="logit"){
    mobj<- stats::glm(Y~.,family = binomial(link = "logit"),data = data.frame(Y,X),...)
  }else if(model=="quantreg"){
    mobj<- quantreg::rq(Y~.,data = data.frame(Y,X),...)
  }else if(model=="probit"){
    mobj<- stats::glm(Y~.,family = binomial(link = "probit"),data = data.frame(Y,X),...)
  }else if(model=="gammainverse"){
    mobj<- stats::glm(Y~.,family = Gamma(link = "inverse"),data = data.frame(Y,X),...)
  }else if(model=="gammaidentity"){
    mobj<- stats::glm(Y~.,family = Gamma(link = "identity"),data = data.frame(Y,X),...)
  }else if(model=="gammalog"){
    mobj<- stats::glm(Y~.,family = Gamma(link = "log"),data = data.frame(Y,X),...)
  }else if(model=="poissonlog"){
    mobj<- stats::glm(Y~.,family = poisson(link = "log"),data = data.frame(Y,X),...)
  }else if(model=="poissonidentity"){
    mobj<- stats::glm(Y~.,family = poisson(link = "identity"),data = data.frame(Y,X),...)
  }else if(model=="poissonsqrt"){
    mobj<- stats::glm(Y~.,family = poisson(link = "sqrt"),data = data.frame(Y,X),...)
  }else if(model=="negbin"){
    mobj<- MASS::glm.nb(Y~.,data = data.frame(Y,X),...)
  }
  return(mobj)
}

#===================================================================================>
#' Construct a network design matrix
#'
#' This function creates the design matrix for a latent network structure using a balanced
#' panel
#'
#' @param datf  the entire data frame of balanced panel with NT rows of unit-time
#' observations
#' @param  Y  dependent variable in the data frame datf
#' @param  X  the covariate(s) generating spillovers
#' @param  Wi  other unit-varying (can be time-invariant) control variables
#' @param  W  global variables. these are only time varying but are common to all units.
#' eg. GDP
#' for individual/state-level data. Note that W has to be a vector of length T so cannot be
#' in the data frame datf
#' @param  panvar  the panel variable eg. unique person/firm identifiers
#' @param  tvar  time variable, eg. years
#' @param  factors  a vector of characters of factors in the data
#' @param scaling  a logical indicating whether non-discrete covariates should be scaled by their standard deviations
#' @param  unicons  a logical indicating whether to include unit-specific constant term

#' @return  Y  vector of dependent variables
#' @return  X  a block matrix of spillover matrix (TN x N^2)
#' @return  Wm  a matrix corresponding to covariate Wi
#' @return  Wf  a matrix of dummies corresponding to factors
#' @export
#'
netdat<- function(datf,Y,X,Wi,W=NULL,panvar,tvar,factors,scaling=TRUE,unicons=TRUE){
  if(any(is.na(data.frame(datf[panvar],datf[tvar])))){stop("NA in panel or time variables unallowed")}
  dat<- datf[order(datf[panvar],datf[tvar]),] #sort data into a panel by units and time
  Nd = nrow(unique(datf[panvar])) # extract number of units
  Td = nrow(unique(datf[tvar])) # extract the number of time periods
  fmat<- dat[factors]
  # check if panel is balanced
  if(dim(dat)[1]!=(Nd*Td)){stop("Unbalanced panel")}
  Xj = matrix(unlist(dat[X]),ncol = Nd) # the covariate(s) generating spillovers or externalities
  # write code to scale continuous covariates and store scales.
  if(scaling){
    fn<- function(x) {stats::sd(stats::na.omit(x))*sqrt(length(which(!is.na(x)))/length(x))}
    vX<- apply(Xj,2,fn);
    if(any(vX<(10^-6))){cat("some units lack variation in x. consider removing them.")}


    if(unicons){
      Xm = kronecker(diag(1,Nd),cbind(1,(Xj/vX)))
    }else{
      Xm = kronecker(diag(1,Nd),(Xj/vX))
    } # a (TdxNd) x (Nd^2) block matrix
  }else{
    if(unicons){
      Xm = kronecker(diag(1,Nd),cbind(1,Xj))
    }else{
      Xm = kronecker(diag(1,Nd),Xj)
    } # a (TdxNd) x (Nd^2) block matrix
  }

  if(all(dim(Xm)!=c((Td*Nd),(Nd^2)))){stop("Network data creation failed.")}

  if(scaling){
    vWm = apply(dat[Wi],2,fn);
    if(any(vWm<(10^-6))){cat("some units lack variation in covariates. consider removing them.")}
    Wm = dat[Wi]/vWm
  }else{
    Wm = dat[Wi]
  }

  Ym = unlist(dat[Y]); zY=rep(1,length(Ym)); zY[!is.na(Ym)]<- Ym[!is.na(Ym)]
  fmod<- apply(fmat, 2, as.factor) #convert the columns of fmat into factors
  Wf=stats::model.matrix(zY~.,data.frame(zY,fmod))[,-1] # remove the constant term
  # combine terms
  #datDM<-data.frame(Ym,Xm,Wm,Wf); names(datDM)[1]<- "Y"
  if(scaling){
    res=list(Y=Ym,X=Xm,Wm=Wm,Wf=Wf,sdX=vX,sdWm=vWm)
  }else{
    res=list(Y=Ym,X=Xm,Wm=Wm,Wf=Wf)
  }
  return(res)
}


#=============================================================================================>
#' Clustering of vector elements
#'
#' A deterministic clustering device of vector elements into k clusters
#'
#' @param k  number of clusters
#' @param vec  the vector of real valued elements
#'
#' @return clus  integer assignment of corresponding elements in vec in up to k clusters
#' @examples
#' ## Not run: set.seed(2); (v=c(rnorm(4,0,0.5),rnorm(3,3,0.5))[sample(1:7)])
#' ## dcluspar(k=2,vec = v)
#' @export
#'
dcluspar<- function(k,vec){
  svec=sort(vec)
  dsvec=diff(svec)
  idub = which(rank(-dsvec)<k)
  idlb = 1+idub
  lb=svec[idlb]
  ub=svec[idub]
  lb = c(min(vec),lb); ub = c(ub,max(vec))
  clus = rep(1,length(vec))
  for (j in 1:k) {
    clus[which(vec>=lb[j] & vec<=ub[j])]=j
  }
  clus
}
#=============================================================================================>
#' Regression via sequential optimisation
#'
#' \code{CCRls} runs regressions with potentially more covariates than observations.
#' See \code{ch.model()} for the list of models supported.
#'
#' @param Y vector of dependent variable Y
#' @param X design matrix (without intercept)
#' @param kap maximum number of parameters to estimate in each active sequential step,
#' as a fraction of the less of total number of observations n or number of covariates p.
#' @param model a string denoting the desired model
#' @param tol level of tolerance for convergence; default \code{tol=1e-06}
#' @param reltol a logical for relative tolerance instead of level. Defaults at TRUE
#' @param rndcov seed for randomising assignment of covariates to partitions; default \code{NULL}
#' @param report number of iterations after which to report progress; default \code{NULL}
#' @param ... additional arguments to be passed to the model
#'
#' @return \code{betas}  parameter estimates (intercept first),
#' @return \code{iter}  number of iterations,
#' @return \code{dev}  increment in the objective function value at convergence
#' @return \code{fval} objective function value at convergence
#'
#' @examples ##Not run:
#' @export
## list(betas=c(bet0,bet_vec),iter=l,dev=dev,fval=val0)
CCRls<- function(Y,X,kap=0.1,model="lm",tol=1e-6,reltol=TRUE,rndcov=NULL,report=NULL,...){
  p = ncol(X)
  n = nrow(X)
  bet_vec<- rep(NA,p) # vector to store parameters, excludes intercept
  asz = 1e-20
  slc<- floor(kap*min(c(n,p))) # maximum size of a local covariate cluster
  lcls<- ceiling((1:p)/slc) # partition covariates into clusters
  nlcls<- max(lcls) #number of local clusters of covariates
  if(!is.null(rndcov)){set.seed(rndcov);lcls<- sample(lcls,p)}

  # initialise parameters
  bet0<- 0 #initialise intercept
  for (j in 1:p) {
    coefs<- coef(ch.model(Y,X[,j],model = model,...))
    bet_vec[j]<- coefs[2]
  }

  val0<- Inf; val1<- 0; l<- 0; dev=(val0-val1)

  while((dev>tol)){
    if(l>0){
      val0<- val1; obj0<- obj1;
      coefs<- coef(obj1)
      bet_vec[IDls]<- coefs[c(2:(1+nIDls))] # local parameter updating
      bet0<- coefs[1]                       # local parameter updating, intercept
      bet_vec[-IDls]<- tail(coefs,n=1)*bet_vec[-IDls] #non-local parameter updating
    }

    l<- l + 1; IND<- l - (ceiling(l/nlcls)-1)*nlcls;

    IDls<- which(lcls==IND); nIDls<- length(IDls)
    # construct (local) design matrix
    XB_ = X[,-IDls]%*%matrix(bet_vec[-IDls],ncol = 1)
    Xl = X[,IDls]
    XX = data.frame(Xl,XB_)
    obj1 = ch.model(Y,as.matrix(XX),model = model,...)
    val1<- -logLik(obj1)
    if(!is.null(report)){if(l%%report==0){ cat("Iter =",l,"fval =",val1,"\n")}}
    if(reltol){dev=(val0-val1)/(asz+abs(val0))}else{dev=(val0-val1)}
    if(l==1){dev=1}
    # 1 in denominator to avoid dividing by zero
  }

  list(betas=c(bet0,bet_vec),iter=l,dev=dev,fval=val0)
}

#=============================================================================================>

clfun2<- function(Y,X,Xnc=NULL,clus,k,model=model,...){
  nrX <- nrow(X) #number of rows of X
  uniClus <- unique(clus)
  X1 <- matrix(NA,nrX,k)
  for(j in uniClus) X1[,j] <- apply(as.matrix(X[,(which(clus == j))]),1,sum)
  if(is.null(Xnc)){
    model1 <-ch.model(Y,X1,model = model,...)
  }else{
    model1 <-ch.model(Y,as.matrix(cbind(Xnc,X1)),model = model,...)
  }
  model1
}

# Next step: build a clustered covariate estimation procedure
# compute BIC at a given k clustering
CCRk<- function(k,betas,Y,X,model="lm",...){
  clus<- dcluspar(k,betas)
  regobj<- clfun2(Y,X,clus,k,model=model,...)
  stats::BIC(regobj)
}

# regression object and clustering at a given number of clusters k
CCRk2<- function(k,betas,Y,X,Xnc=NULL,model="lm",...){
  clus<- dcluspar(k,betas)
  regobj<- clfun2(Y,X,Xnc=Xnc,clus,k,model=model,...)
  list(BIC=stats::BIC(regobj),regobj=regobj,clus=clus)
}




# Golden Search Algorithm for the outer loop ()
goldensearch<- function(fn,interval,tol=1){
  a=min(interval); b = max(interval)
  xvals<- c(0); xvals[1]<- a; xvals[2]<- b
  fvals<- c(0)
  fa = fn(a); fb = fn(b);
  fvals[1]<- fa; fvals[2]<- fb
  cnt=2; # set counter to 2 for function evaluations
  phi<- (1+sqrt(5))/2

  c = ceiling(b - (b-a)/phi); d = floor(a + (b-a)/phi);
  if(any(xvals==c)){
    c=xvals[which(xvals==c)]; fc=fvals[which(xvals==c)]
  }else{
    cnt=cnt+1
    fc=fn(c)
  }
  if(any(xvals==d)){
    d=xvals[which(xvals==d)]; fd=fvals[which(xvals==d)]
  }else{
    cnt=cnt+1
    fd=fn(d)
  }


  l = 1; # set counter for iterations
  cat("iter = ",l,"\n");arreter = 0
  while(abs(c-d)>tol & arreter==0){# while c d
    if(fc<fd){
      b=d; fb=fd;d=c;fd=fc;
      c = ceiling(b-(b-a)/phi);
      if(any(xvals==c)){
        c=xvals[which(xvals==c)]; fc=fvals[which(xvals==c)]
      }else{
        cnt=cnt+1
        fc=fn(c)
      }

    }else if(fc>fd){
      a=c; fa = fc; c=d; fc=fd
      d = floor(a + (b-a)/phi);
      if(any(xvals==d)){
        d=xvals[which(xvals==d)]; fd=fvals[which(xvals==d)];
      }else{
        cnt=cnt+1
        fd=fn(d)
      }
    }else{
      arreter=1
    }
    l=l+1; cat("iter = ",l,"\n")
  }

  optx=ifelse(fc>fd,d,c)
  res=list(k=optx,value=min(fd,fc),iter=l,iterfn=cnt)
  return(res)
}
