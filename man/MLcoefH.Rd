\name{MLcoefH}
\alias{MLcoefH}

\title{
    Two-level scalability coefficents H
}

\description{
    Computes all types of two-level scalability coefficients (Snijders, 2001; Crisan, Van de Pol, & Van der Ark, 2016),
    that is, between- and within-rater item-pair scalability coefficents and their ratio (HBij, HWij, and BWij, respectively),
    between- and within-rater item scalability coefficients and their ratio (HBi, HWi, and BWi, respectively),
    and between- and within-rater total scale scalability coefficients and their ratio (HB, HW, and BW, respectively).
    In addition, standard errors are estimated (Koopman, Zijlstra, & Van der Ark, 2019). Note that this version is an 
    adaptation of the estimation methods described in Snijders and in Koopman et al., as the group proportions are now by 
    default weighted for group size (Koopman, Zijlstra, & Van der Ark, 2020). As a result, the estimates for the 
    autonomySupport data differs from the printed table in Koopman et al. (2019).
}

\usage{
    MLcoefH(X, se = TRUE, nice.output = TRUE, subject = 1, fixed.itemstep.order = NULL, weigh.props = TRUE)
}

\arguments{
  \item{X}{matrix or data frame of numeric data containing a subject indicator column and the
           responses of \code{nrow(X)} raters to \code{ncol(X) - 1} items. 
           Missing values are not allowed}
  \item{se}{Logical: If \code{TRUE}, the standard errors are printed alongside the scalability coefficients}
  \item{nice.output}{Logical: If \code{TRUE}, scalability coefficients and standard errors are combined 
           in an a single object of class \code{noquote}. Item-pair ratios BWij are only printed if \code{FALSE}}
  \item{subject}{Represents the subject column. Default is column 1.}
  \item{fixed.itemstep.order}{matrix with number of rows equal to the number of item steps (m) and number of columns equal to the number of items (J). The matrix should consis the integers 1 : (m * J), indicating a predefined order of the items steps with respect to popularity. Value 1 indicates the easiest (most popular) item step, value (m * J) indicates the most difficult item step.}
  \item{weight.props}{If \code{TRUE}: Use weighted proportions across groups to estimate coefficients and standard errors, if \code{FALSE}: Use averaged proportions across groups to estimate coefficients and standard errors.}
}

\details{
 
 If \code{se = TRUE} and \code{nice.output = TRUE}, the result is a list of 3 objects of class \code{noquote};  
 if \code{se = TRUE} and \code{nice.output = FALSE}, the result is a list of 3 matrices, one per set of coefficients; and 
 if \code{se = FALSE}, the result is a list of 3 matrices containing only the scalability coefficients.
 Computation of standard errors can be slow for a combination of many subjects and a large number of items.
}

\value{
 \item{Hij}{Scalability coefficients of the item pairs, upper triangle are the between-rater coefficients, and the
            lower triangle the within-rater coefficients (possibly with standard errors; see details). 
            If \code{nice.output = FALSE} this returns a matrix with the type of coefficients per column}  
 \item{Hij}{Within-rater scalability coefficients of the item pairs (possibly with standard errors; see details)}
 \item{Hi}{Between- and within-rater scalability coefficients of the items and their ratio (possibly with standard errors; see details)}  
 \item{H}{between- and within-rater scalability coefficients of the entire scale and their ratio (possibly with standard error; see details)}  
}

\references{
  Crisan, D. R., Van de Pol, J. E., & Van der Ark, L. A. (2016). Scalability coefficients for two-level polytomous item scores: 
  An introduction and an application. In L. A. Van der Ark, D. M. Bolt, W.-C. Wang, & M. Wiberg (Eds.), 
  \emph{Quantitative psychology research}. New York, NY: Springer. doi: 10.1007/978-3-319-38759-8_11
  
  Koopman, L., Zijlstra, B. J. H. & Van der Ark, L. A., (2019). 
  Standard errors of two-level scalability coefficients. \emph{British Journal of Statistical and Mathematical Psychology}.
  Advance online publication. doi: 10.1111/bmsp.12174 

  Koopman, L. Zijlstra, B. J. H, & Van der Ark, L. A. (2020). 
  \emph{A two-step procedure for scaling multilevel data using Mokken's scalability coefficients}. Manuscript submitted 
     for publication.

  Snijders, T. A. B. (2001). Two-level non-parametric scaling for dichotomous data. In A. Boomsma, 
  M. A. J. van Duijn, & T. A. B. Snijders (Eds.), \emph{Essays on item response theory} (p. 319-338). 
  New York, NY: Springer. doi: 10.1007/978-1-4613-0169-1_17
}

\author{
   Letty Koopman \email{V.E.C.Koopman@uva.nl}
   L. Andries van der Ark \email{L.A.vanderArk@uva.nl} 
}

\seealso{
\code{\link{coefH}}, 
}


\examples{

# Small data example (Koopman et al., 2019)
smallData <- data.frame(Subs = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3),
                        Xa   = c(0, 0, 1, 0, 1, 1, 1, 2, 1, 0, 1, 2, 0, 0, 0), 
                        Xb   = c(0, 0, 1, 0, 2, 2, 2, 1, 2, 1, 2, 2, 1, 1, 0))

MLcoefH(smallData)

# Load real data example. Note that due to an estimation adaptation (Koopman et al., 2020)
# the results differ from the table in Koopman et al. (2019).
data(autonomySupport)

# Compute scalability coefficients with or without standard errors and nice output
\donttest{          
H.se.nice <- MLcoefH(autonomySupport)
H.se.nice 
}
\donttest{          
H.se.not  <- MLcoefH(autonomySupport, nice.output = FALSE)
H.se.not  
}

\donttest{          
H.nice    <- MLcoefH(autonomySupport, se = FALSE)
H.nice    
}

\donttest{          
H.not     <- MLcoefH(autonomySupport, se = FALSE, nice.output = FALSE)
H.not     
}

}
\keyword{two-level scalability coefficients}
