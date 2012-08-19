\name{l1svm}
\alias{l1svm}
\title{
Training function of L1 norm SVM
}
\description{
The classifier is learned using training data.
}
\usage{
l1svm(x, y, lambda = NULL, rho = 1, thol = 1e-04, maxIter = 1e4, rescale = T)
}
\arguments{
  \item{x}{
The training dataset represented in a \code{n} by \code{d} matrix, where \code{n} is sample size and \code{d} is dimension.
}
  \item{y}{
The labels of training dataset represented in a \code{n} by \code{d} matrix, where \code{n} is sample size and \code{d} is dimension. The label MUST be encoded as either 1 or -1.
}
  \item{lambda}{
A sequence of decresing positive numbers to control the regularization. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence \code{lambda = log(d)*c(20:1)/100/n}. Users can also specify a sequence to override this. Please use with care - a "bad" sequence may generate an full-zero path.
}
  \item{rho}{
The penalty parameter used in the optimization algorithm. Please use it with care - setting \code{rho} too large or too small may make the algorithm early stop or diverge.
}
  \item{thol}{
Stopping Criterion: Maximum relative change of primal and dual parameters.
}
  \item{maxIter}{
The number of maximum iterations.
}
  \item{rescale}{
Rescale all variables in the training dataset to 0-1 scale if \code{rescale=T}. We strongly recommond users NOT to turn it off.
}
}
\details{
We adopt the linearized alternative direction method of multipliers. The computation is further accelerated by "warm-start" and "active-set" tricks. Based on our experience, the optimal regularization parameters may dramatically change under different settings. The default lambda sequence cannot secure the performance.
}
\value{
  \item{rescale}{
The indicator of the rescaling. 
}
  \item{x.min}{
A vector with each entry corresponding to the minimum of each input variable. (Used for rescaling in testing)
}
  \item{x.max}{
A vector with each entry corresponding to the maximum of each input variable. (Used for rescaling in testing)
}
  \item{lambda}{
A sequence of regularization parameter used in training.
}
  \item{w}{
The solution path matrix (\code{d} by length of \code{lambda}) with each column corresponding to a regularization parameter.
}
  \item{b}{
The solution path of the intercept.
}
  \item{df}{
The degree of freedom of the solution path (The number of non-zero parameters)
}
  \item{lab}{
Predicted labels represented in a \code{n} by \code{d} matrix, where \code{n} is sample size and \code{d} is dimension. 
}
  \item{dec}{
Decision values represented in a \code{n} by \code{d} matrix, where \code{n} is sample size and \code{d} is dimension. 
}
}
\references{
T. Zhao and H.Liu. "Sparse Additive Machine", \emph{International Conference on Artificial Intelligence and Statistics}, 2012.\cr
P. Bradley and O. Mangasarian. ""Feature selection via concaveminimization and support vector machines"", \emph{International Conference on Machine Learing}, 1998.
}
\author{
Tuo Zhao, Xingguo Li, Han Liu, Lie Wang and Kathryn Roeder\cr
Maintainers: Tuo Zhao<tourzhao@gmail.com>
}
\seealso{
\code{\link{SAM}},\code{\link{spam}},\code{\link{plot.l1svm},\link{print.l1svm},\link{predict.l1svm}}
}

\examples{
## generating training data
x = rbind(0.5+matrix(rnorm(100),50,2),-0.5+matrix(rnorm(100),50,2))
x = cbind(x,matrix(rnorm(800),100,8))

## generating labels
y = c(rep(1,50),rep(-1,50))

## Training
fit = l1svm(x,y)
fit

## plotting solution path
plot(fit)

## generating testing data
xt = rbind(0.5+matrix(rnorm(200),50,2),-0.5+matrix(rnorm(200),50,2))
xt = cbind(xt,matrix(rnorm(1600),100,8))

## predicting labels
out = predict(fit,newdata=xt)
}