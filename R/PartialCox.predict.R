### PartialCox.predict.R  (2012-04-27)
###
###    Partial Cox PH model predictions
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

# Given a trained model, apply the model to test data and return risk scores (negative values represent low risk and positive values high risk)
PartialCox.predict = function(Xtest, PartialCoxModel){
  if (!is.matrix(Xtest)) stop("Test data must be given as a matrix with samples in rows!")
  out = sweep(Xtest, 2, PartialCoxModel$predcoef[,1]) %*% PartialCoxModel$predcoef[,2]
  return (out)
}
