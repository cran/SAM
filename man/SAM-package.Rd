\name{SAM-package}
\alias{SAM-package}
\alias{SAM}
\docType{package}
\title{
Sparse Additive Modelling
}
\description{
The package SAM targets at high dimensional predictive modeling (regression and classification) for complex data analysis. SAM is short for sparse additive modeling, and adopts the computationally efficient basis spline technique. We solve the optimization problems by various computational algorithms including the block coordinate descent algorithm, fast iterative soft-thresholding algorithm, and newton method. The computation is further accelerated by warm-start and active-set tricks.
}
\details{
\tabular{ll}{
Package: \tab SAM\cr
Type: \tab Package\cr
Version: \tab 1.0.4\cr
Date: \tab 2013-09-17\cr
License: \tab GPL-2 \cr
}
}
\references{
P. Ravikumar, J. Lafferty, H.Liu and L. Wasserman. "Sparse Additive Models", \emph{Journal of Royal Statistical Society: Series B}, 2009.\cr
T. Zhao and H.Liu. "Sparse Additive Machine", \emph{International Conference on Artificial Intelligence and Statistics}, 2012.\cr
}
\author{
Tuo Zhao, Xingguo Li, Han Liu, and Kathryn Roeder\cr
Maintainers: Tuo Zhao<tourzhao@gmail.com>
}
\seealso{
\code{\link{samQL}},\code{\link{samHL}},\code{\link{samLL}},\code{\link{samEL}}
}