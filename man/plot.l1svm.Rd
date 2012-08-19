\name{plot.l1svm}
\alias{plot.l1svm}

\title{
Plot function for S3 class \code{"l1svm"}
}
\description{
This function plots the regularization path (regularization parameter versus coefficient)
}
\usage{
\method{plot}{l1svm}(x,...)
}
\arguments{
	\item{x}{
An object with S3 class \code{"l1svm"} 
}
	\item{\dots}{
System reserved (No specific usage)
}
}
\details{
The horizontal axis is for the regularization parameter index. The vertical axis is for the value of each coefficient.
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
\code{\link{l1svm}}
}
