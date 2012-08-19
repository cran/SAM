\name{plot.spam}
\alias{plot.spam}

\title{
Plot function for S3 class \code{"spam"}
}
\description{
This function plots the regularization path (regularization parameter versus functional norm)
}
\usage{
\method{plot}{spam}(x,...)
}
\arguments{
	\item{x}{
An object with S3 class \code{"spam"} 
}
	\item{\dots}{
System reserved (No specific usage)
}
}
\details{
The horizontal axis is for the regularization parameter index. The vertical axis is for the functional norm of each component.
}
\references{
T. Zhao and H.Liu. "Sparse Additive Machine", \emph{International Conference on Artificial Intelligence and Statistics}, 2012.\cr
P. Bradley and O. Mangasarian. ""Feature selection via concaveminimization and support vector machines"", \emph{Intermational Conference on Machine Learing}, 1998.
}
\author{
Tuo Zhao, Xingguo Li, Han Liu, Lie Wang and Kathryn Roeder\cr
Maintainers: Tuo Zhao<tourzhao@gmail.com>
}
\seealso{
\code{\link{spam}}
}
