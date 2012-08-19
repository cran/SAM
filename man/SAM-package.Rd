\name{SAM-package}
\alias{SAM-package}
\alias{SAM}
\docType{package}
\title{
Sparse Additive Machine
}
\description{
The package "SAM" targets at high dimensional classification problem for complex data analysis. SAM is short for "Sparse Additive Machine", which is a explicit combination of high dimensional sparse additive modeling and support vector machine. Different from existing non-linear classification methods, which usually use , SAM adopts the computationally efficient basis spline technique. The optimization is solved by the Linearized Alternative Direction Method of Multipliers (L-ADMM). The computation is further accelerated by "warm-start" and "active-set" tricks. For users who are interested in large-scale problems, we also provide an implementation of L1 norm SVM for computational convenience.
}
\details{
\tabular{ll}{
Package: \tab SAM\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2012-08-20\cr
License: \tab GPL-2 \cr
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
\code{\link{spam},\link{l1svm}}
}