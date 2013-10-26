\name{PartialCox}
\alias{PartialCox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function trains a partial Cox model on given training data matrix and survival object.
}
\description{
This function implements the partial Cox model training as proposed by Li and Gui (2004). The model order is controlled by ncomp, which by default is set to 3 latent variables. Increasing the number of latent variables increases runtime radically.
}
\usage{
PartialCox(Xtrain, SurvObject, ncomp = 3, univbetas = c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xtrain}{
Training data matrix with samples in rows and variables in columns.
}
  \item{SurvObject}{
This is an object of type Surv and describes the follow-up times and event indicators.
}
  \item{ncomp}{
Number of latent variables to evaluate (analogous to number of principal components). By default 3. Incresing ncomp radically increases runtime for high dimensional data.
}
  \item{univbetas}{
For certain feature selection scenarios this variable can be useful for passing the univariate Cox coefficients that do not change across iterations.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{beta}{The partial Cox model coefficients}
\item{univbetas}{The univariate Cox model coefficients}
\item{loadings}{Variable loadings}
\item{predcoef}{A list containing beta and means from the training data}
}
\references{
Li, H. and J. Gui (2004), “Partial Cox regression analysis for high-dimensional microarray gene expression data,” Bioinformatics, 20(Suppl 1), i208–i215.
}
\author{
Miika Ahdesmaki
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# generate artificial training data
N=50
P=100
Xtrain = matrix(rnorm(N*P),N,P)
# true coefficients, based on which survival times are generated
betas = matrix(0,P,1)
betas[c(1,10,20,100)] = c(-5,5,-3,3)
mySurvTimes = Xtrain %*% betas
mySurvTimes = mySurvTimes - min(mySurvTimes)
# event indicator
myEvents = matrix(1,N,1)
# censor samples 20-40
myEvents[20:40] = 0
# generate a survival object
mySurvObject = Surv(time=mySurvTimes, event=myEvents)
# start partial Cox inference
PCoxModel = PartialCox(Xtrain=Xtrain, SurvObject=mySurvObject)
# generate a test matrix
Xtest = matrix(rnorm(N*P),N,P)
testSurvTimes = Xtest %*% betas
# predict risk scores (high risk score implies a low survival time)
preds = PartialCox.predict(Xtest,PCoxModel)
require('Hmisc')
# evaluate concordance between predicted risk and true survival times
print("Concordance index between predictions and survival times for artificial test data:")
print(rcorr.cens(-preds, testSurvTimes)["C Index"])

## The function is currently defined as
PartialCox = function(Xtrain, SurvObject, ncomp=3, univbetas=c())
{
  #require("survival")
  if (!is.matrix(Xtrain)) stop("Training data must be given as matrix!")
  if (!is.Surv(SurvObject)) stop("SurvObject must be given as an object of Surv!")
  #if (missing(SurvTimes)) stop("Survival times are missing!")
  #if (missing(EventIndicator)) stop("Event indicator vector missing!")
  Xdim = dim(Xtrain)
  if (min(Xdim)<ncomp){
    padding = TRUE
    ncompeval = min(Xdim) # to avoid numerical trouble 
  }else{
    padding = FALSE
    ncompeval = ncomp
  }
  #SurvObj = Surv(time=SurvTimes, event=EventIndicator) # form a Surv object for further use in coxph
  
  # fit partial Cox
  means = colMeans(Xtrain)
  Xcen = sweep(Xtrain,2,means,'-')
  rownames(Xcen) = rownames(Xcen, do.NULL = FALSE, prefix = "Obs.")
  control = coxph.control()
  if (is.null(univbetas)){
    univbetas = matrix(0,Xdim[2])
    for (ii in 1:Xdim[2]){
      #test1 = list(Xcen=Xcen[,ii])
      #model = coxph(SurvObj ~ Xcen, test1)
      #univbetas[ii] = model$coefficients
      modelfit = suppressWarnings(coxph.fit(x=Xcen[,ii, drop = F],y=cbind(SurvTimes,EventIndicator),strata = c(), offset = c(), init = c(), control = control, weights = c(), method = "efron", rownames = rownames(Xcen)))
      univbetas[ii] = modelfit$coefficients
    }
  }
  # initialise
  loadings = matrix(0,Xdim[2],ncompeval) # initialise loadings matrix
  ProjMat = array(0,c(Xdim[2],Xdim[2],ncompeval)) # memory hog
  ProjMat[,,1] = diag(1,Xdim[2],Xdim[2]) # initialise as identity matrix
  b = matrix(0,Xdim[2],ncompeval) # these are "pseudo betas"
  b[,1] = univbetas
  w = matrix(0,Xdim[2],ncompeval)
  w[,1] = as.matrix(apply(Xcen,2,var),nrow=Xdim[2])
  w[,1] = w[,1]/sum(w[,1]) # weights (p-by-1) proportional to variances, sum to one
  
  # reserve space for T:s (N-by-1 vectors, NumComp of them), initialise
  T = matrix(0,Xdim[1],ncompeval);
  T[,1] = crossprod(t(Xcen),w[,1]*b[,1]);

  # iterate
  if(ncompeval>=2){
    for (jjj in 2:ncompeval){
      wb = matrix(w[,jjj-1]*b[,jjj-1],nrow=Xdim[2]) # variances times betas from previous round
      tempRes = crossprod(wb, crossprod(crossprod(t(Xcen),ProjMat[,,jjj-1]))) # ProjMat{jjj-1}'*VTV*ProjMat{jjj-1}; % store temporary result (for next line of code)
      
      H = crossprod(t(ProjMat[,,jjj-1]), crossprod(t(wb), tempRes)) / as.vector(crossprod(t(tempRes),wb))
      ProjMat[,,jjj] = ProjMat[,,jjj-1] - H # projection matrix
      
      Vcur = crossprod(t(Xcen) , ProjMat[,,jjj])
      w[,jjj] = as.matrix(apply(Vcur,2,var),nrow=Xdim[2])
      w[,jjj] = w[,jjj]/sum(w[,jjj]) # update weights
      for (kkk in 1:Xdim[2]){
	  temp = suppressWarnings(coxph.fit(x=cbind(T[,1:(jjj-1)],Vcur[,kkk]), y=cbind(SurvTimes,EventIndicator),strata = c(), offset = c(), init = c(), control = control, weights = c(), method = "efron", rownames = rownames(Xcen)))
	  b[kkk,jjj] = temp$coefficients[length(temp$coefficients)] # pick the last one, corresponding to the left-over data
      }
      T[,jjj] = crossprod(t(Vcur), w[,jjj,drop=F]*b[,jjj,drop=F] )       
    }
  }
  
  # fit T:s to the data and reverse engineer weights for the original Xcen
  Tb = suppressWarnings(coxph.fit(x=T, y=cbind(SurvTimes,EventIndicator),strata = c(), offset = c(), init = c(), control = control, weights = c(), method = "efron", rownames = rownames(Xcen)))
  loadings[,1] = as.vector(Tb$coefficients[1]) * w[,1,drop=F]*b[,1,drop=F]
  if(ncompeval>=2){
    for (lll in 2:ncompeval){
      loadings[,lll] = Tb$coefficients[lll] * crossprod(t(ProjMat[,,lll]), w[,lll,drop=F]*b[,lll,drop=F])
    }
  }

  
  beta = as.matrix(rowSums(loadings),ncol=1)
  out = list(beta=beta, univbetas=univbetas, loadings=loadings, 
             predcoef=cbind(means, beta))
  class(out)="PartialCox"
  return (out)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
