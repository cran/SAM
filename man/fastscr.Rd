\name{fastscr}
\alias{fastscr}
\title{
Screening variables before the traning stage
}
\description{
Variables are preselected by t test ot wilcoxon test
}
\usage{
fastscr(x, y, nscr = NULL, method = "t.test")
}
\arguments{
  \item{x}{
The training dataset represented in a \code{n} by \code{d} matrix, where \code{n} is fastscrple size and \code{d} is dimension.
}
  \item{y}{
The labels of training dataset represented in a \code{n} by \code{d} matrix, where \code{n} is fastscrple size and \code{d} is dimension. The label must be encoded as either 1 or -1.
}
  \item{nscr}{
The number of variables retained. The default value is \code{min(n,d)}.
}
  \item{method}{
The screening will be based on t.test if \code{method = "t.test")}. An alternative is wilcoxon test if \code{method = "wilcox")}.
}
}
\details{
The fast screening can greatly reduce the compuational burden, and often achieve equally (or better) performance in practice.
}
\value{
  \item{idx}{
A vector of length \code{nscr} corresponding to the indices of retained variables. 
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
\code{\link{SAM}},\code{\link{l1svm}},\code{\link{spam}}
}

\examples{
## generating training data
x = rbind(0.5+matrix(rnorm(100),50,2),-0.5+matrix(rnorm(100),50,2))
x = cbind(x,matrix(rnorm(7800),100,78))

## generating labels
y = c(rep(1,50),rep(-1,50))

## fast screening
idx = fastscr(x,y,nscr = 10)

## Training
## fit = spam(x[,idx],y,k=4)
}