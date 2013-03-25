\name{plot.samLL}
\alias{plot.samLL}

\title{
Plot function for S3 class \code{"samLL"}
}
\description{
This function plots the regularization path (regularization parameter versus functional norm)
}
\usage{
\method{plot}{samLL}(x,...)
}
\arguments{
	\item{x}{
An object with S3 class \code{"samLL"} 
}
	\item{\dots}{
System reserved (No specific usage)
}
}
\details{
The horizontal axis is for the regularization parameters in log scale. The vertical axis is for the functional norm of each component.
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
\code{\link{samLL}}
}
