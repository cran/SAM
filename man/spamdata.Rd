\name{spamdata}
\alias{spamdata}
\docType{data}
\title{SPAM Detection Dataset}
\description{
This data set consists of 4601 emails (1813 spam and 2788 regular).
}
\usage{data(spamdata)}
\format{
  The format is a list containing conatins a matrix and a vector.
  1. fea - 4601x57: 4601 emails with 57 word counts.
  2. lab - 4601x1: "1" denotes "spam" and "-1" denotes "regular".
}
\details{
This data set can be used to perform the sparse additive machine to demonstrate the supriority of its non-linearity over L1 norm SVM.
}
\source{
It is publicly available at http://archive.ics.uci.edu/ml/datasets/Spambase
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
data(spamdata)
dim(spamdata$fea)
sum(spamdata$lab==1)
}
