\name{print.spam}
\alias{print.spam}
\title{
Printing function for S3 class \code{spam}
}
\description{
Summarize the information of the object with S3 class \code{spam}.
}
\usage{
\method{print}{spam}(x, ...)
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
The output includes length and d.f. of the regularization path, and whether the input data is rescaled.
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
