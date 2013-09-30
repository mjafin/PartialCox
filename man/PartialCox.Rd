\name{PartialCox}
\alias{PartialCox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PartialCox(Xtrain, SurvTimes, EventIndicator, ncomp = 3, univbetas = c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xtrain}{
%%     ~~Describe \code{Xtrain} here~~
}
  \item{SurvTimes}{
%%     ~~Describe \code{SurvTimes} here~~
}
  \item{EventIndicator}{
%%     ~~Describe \code{EventIndicator} here~~
}
  \item{ncomp}{
%%     ~~Describe \code{ncomp} here~~
}
  \item{univbetas}{
%%     ~~Describe \code{univbetas} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Xtrain, SurvTimes, EventIndicator, ncomp = 3, univbetas = c()) 
{
    require("survival")
    if (!is.matrix(Xtrain)) 
        stop("Training data must be given as matrix!")
    if (missing(SurvTimes)) 
        stop("Survival times are missing!")
    if (missing(EventIndicator)) 
        stop("Event indicator vector missing!")
    Xdim = dim(Xtrain)
    if (min(Xdim) < ncomp) {
        padding = TRUE
        ncompeval = min(Xdim)
    }
    else {
        padding = FALSE
        ncompeval = ncomp
    }
    means = colMeans(Xtrain)
    Xcen = sweep(Xtrain, 2, means, "-")
    rownames(Xcen) = rownames(Xcen, do.NULL = FALSE, prefix = "Obs.")
    control = coxph.control()
    if (is.null(univbetas)) {
        univbetas = matrix(0, Xdim[2])
        for (ii in 1:Xdim[2]) {
            modelfit = suppressWarnings(coxph.fit(x = Xcen[, 
                ii, drop = F], y = cbind(SurvTimes, EventIndicator), 
                strata = c(), offset = c(), init = c(), control = control, 
                weights = c(), method = "efron", rownames = rownames(Xcen)))
            univbetas[ii] = modelfit$coefficients
        }
    }
    loadings = matrix(0, Xdim[2], ncompeval)
    ProjMat = array(0, c(Xdim[2], Xdim[2], ncompeval))
    ProjMat[, , 1] = diag(1, Xdim[2], Xdim[2])
    b = matrix(0, Xdim[2], ncompeval)
    b[, 1] = univbetas
    w = matrix(0, Xdim[2], ncompeval)
    w[, 1] = as.matrix(apply(Xcen, 2, var), nrow = Xdim[2])
    w[, 1] = w[, 1]/sum(w[, 1])
    T = matrix(0, Xdim[1], ncompeval)
    T[, 1] = crossprod(t(Xcen), w[, 1] * b[, 1])
    if (ncompeval >= 2) {
        for (jjj in 2:ncompeval) {
            wb = matrix(w[, jjj - 1] * b[, jjj - 1], nrow = Xdim[2])
            tempRes = crossprod(wb, crossprod(crossprod(t(Xcen), 
                ProjMat[, , jjj - 1])))
            H = crossprod(t(ProjMat[, , jjj - 1]), crossprod(t(wb), 
                tempRes))/as.vector(crossprod(t(tempRes), wb))
            ProjMat[, , jjj] = ProjMat[, , jjj - 1] - H
            Vcur = crossprod(t(Xcen), ProjMat[, , jjj])
            w[, jjj] = as.matrix(apply(Vcur, 2, var), nrow = Xdim[2])
            w[, jjj] = w[, jjj]/sum(w[, jjj])
            for (kkk in 1:Xdim[2]) {
                temp = suppressWarnings(coxph.fit(x = cbind(T[, 
                  1:(jjj - 1)], Vcur[, kkk]), y = cbind(SurvTimes, 
                  EventIndicator), strata = c(), offset = c(), 
                  init = c(), control = control, weights = c(), 
                  method = "efron", rownames = rownames(Xcen)))
                b[kkk, jjj] = temp$coefficients[length(temp$coefficients)]
            }
            T[, jjj] = crossprod(t(Vcur), w[, jjj, drop = F] * 
                b[, jjj, drop = F])
        }
    }
    Tb = suppressWarnings(coxph.fit(x = T, y = cbind(SurvTimes, 
        EventIndicator), strata = c(), offset = c(), init = c(), 
        control = control, weights = c(), method = "efron", rownames = rownames(Xcen)))
    loadings[, 1] = as.vector(Tb$coefficients[1]) * w[, 1, drop = F] * 
        b[, 1, drop = F]
    if (ncompeval >= 2) {
        for (lll in 2:ncompeval) {
            loadings[, lll] = Tb$coefficients[lll] * crossprod(t(ProjMat[, 
                , lll]), w[, lll, drop = F] * b[, lll, drop = F])
        }
    }
    beta = as.matrix(rowSums(loadings), ncol = 1)
    out = list(beta = beta, univbetas = univbetas, loadings = loadings, 
        predcoef = cbind(means, beta))
    class(out) = "PartialCox"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
