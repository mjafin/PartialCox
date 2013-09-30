### PartialCox.R  (2012-04-27)
###
###    Partial Cox PH model training
###
### Copyright 2012-04 Miika Ahdesmaki
###
###
### This file is part of the `PartialCox' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

PartialCox = function(Xtrain, SurvTimes, EventIndicator, ncomp=3, univbetas=c())
{
  require("survival")
  #require("Hmisc")
  if (!is.matrix(Xtrain)) stop("Training data must be given as matrix!")
  if (missing(SurvTimes)) stop("Survival times are missing!")
  if (missing(EventIndicator)) stop("Event indicator vector missing!")
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
