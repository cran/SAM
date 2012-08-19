\name{predict.spam}
\alias{predict.spam}
\title{
Prediction function for S3 class \code{"spam"}
}
\description{
Predict the class labels for testing data.
}
\usage{
\method{predict}{spam}(object, newdata, ...)
}
\arguments{
  \item{object}{
An object with S3 class \code{"spam"}.
}
  \item{newdata}{
The testing dataset represented in a \code{n} by \code{d} matrix, where \code{n} is spample size and \code{d} is dimension.
}
	\item{\dots}{
System reserved. (No specific usage)
}
}
\details{
The testing dataset is rescale to the spame range, and expanded by the spame spline basis funcions as the training data.
}
\value{
  \item{lab}{
Predicted labels also represented in a \code{n} by \code{d} matrix, where \code{n} is spample size and \code{d} is dimension. 
}
  \item{dec}{
Decision values also represented in a \code{n} by \code{d} matrix, where \code{n} is spample size and \code{d} is dimension. 
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
\code{\link{spam}}
}