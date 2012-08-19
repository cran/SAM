\name{predict.l1svm}
\alias{predict.l1svm}
\title{
Prediction function for S3 class \code{"l1svm"}
}
\description{
Predict the class labels for testing data.
}
\usage{
\method{predict}{l1svm}(object, newdata, ...)
}
\arguments{
  \item{object}{
An object with S3 class \code{"l1svm"}.
}
  \item{newdata}{
The testing dataset represented in a \code{n} by \code{d} matrix, where \code{n} is sample size and \code{d} is dimension.
}
	\item{\dots}{
System reserved. (No specific usage)
}
}
\details{
The testing dataset is rescale to the l1svme range, and expanded by the l1svme spline basis funcions as the training data.
}
\value{
  \item{lab}{
Predicted labels also represented in a \code{n} by \code{d} matrix, where \code{n} is sample size and \code{d} is dimension. 
}
  \item{dec}{
Decision values also represented in a \code{n} by \code{d} matrix, where \code{n} is sample size and \code{d} is dimension. 
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
\code{\link{l1svm}}
}