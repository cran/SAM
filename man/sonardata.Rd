\name{sonardata}
\alias{sonardata}
\docType{data}
\title{Sonar Dataset}
\description{
This data set consists of 208 sonar signals (111 from metal and 97 from rock).
}
\usage{data(sonardata)}
\format{
  The format is a list containing conatins a matrix and a vector.
  1. fea - 208x60: 208 signals from 60 different angles.
  2. lab - 60x1: "1" denotes "rock" and "-1" denotes "metal".
}
\details{
This data set can be used to perform the sparse additive machine to demonstrate the supriority of its non-linearity over L1 norm SVM.
}
\source{
It is publicly available at http://archive.ics.uci.edu/ml/datasets/Connectionist+Bench+(Sonar,+Mines+vs.+Rocks)
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
\code{\link{spam}},\code{\link{l1svm}}
}
\examples{
data(sonardata)
dim(sonardata$fea)
sum(sonardata$lab==1)
}
