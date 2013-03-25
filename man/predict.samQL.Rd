\name{predict.samQL}
\alias{predict.samQL}
\title{
Prediction function for S3 class \code{"samQL"}
}
\description{
Predict the responses for testing data.
}
\usage{
\method{predict}{samQL}(object, newdata, ...)
}
\arguments{
  \item{object}{
An object with S3 class \code{"samQL"}.
}
  \item{newdata}{
The testing dataset represented in a \code{n} by \code{d} matrix, where \code{n} is testing sample size and \code{d} is dimension.
}
	\item{\dots}{
System reserved. (No specific usage)
}
}
\details{
The testing dataset is rescale to the samQLe range, and expanded by the samQLe spline basis funcions as the training data.
}
\value{
  \item{values}{
Predicted values also represented in a \code{n} by the length of \code{lambda} matrix, where \code{n} is testing sample size. 
}
}
\references{
P. Ravikumar, J. Lafferty, H.Liu and L. Wasserman. "Sparse Additive Models", \emph{Journal of Royal Statistical Society: Series B}, 2009.\cr
T. Zhao and H.Liu. "Sparse Additive Machine", \emph{International Conference on Artificial Intelligence and Statistics}, 2012.\cr
}
\author{
Tuo Zhao, Xingguo Li, Han Liu, Kathryn Roeder\cr
Maintainers: Tuo Zhao<tourzhao@gmail.com>
}
\seealso{
\code{\link{samQL}}
}