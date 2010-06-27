\name{LRRbreak.test}
\alias{LRRbreak.test}
\title{Test of cointegration with breaks}
\description{Test of no cointegration against cointegration wih breaks of gregory and Hansen (1996)
}
\usage{
LRRbreak.test(y,x,model=c("C","CT","CS"), AdfLagChoice=c("User", "AIC", "BIC", "P-Val"),ADFlagmax=6, trim=0.15, hpc=c("none", "foreach")) }

\arguments{
\item{y}{A vector. The left hand-side variable in the cointegration relationship}
\item{x}{A vector or matrix. The right hand-side variable(s) in the cointegration relationship}
\item{model}{The type of break to test for. See details. }
\item{AdfLagChoice}{The procedure to use to determine the lag order for the ADF-like test. }
\item{ADFlagmax}{The lag order for ADF in case AdfLagChoice='user'}
\item{trim}{trimming parameter indicating the minimal percentage of observations in each regime}
\item{hpc}{Experimental. Whether algorithms should be run in parallel. }
}

\value{A list of containing may outputs, among them the critical values and break dates. 
}

\details{Gregory and hansen suggest sup-statistics for the break, and provide critical values, for three different breaks structure:

\tabular{rl}{
\code{C}{Level shift: only shift in the constant}\cr
\code{CT}{Level shift, trend present: only shift in the constant, but a trend is added in the equation}\cr
\code{CS}{Level and slope shift: The shift is both in the constant and coinetgarting parameter. No trend is present in the model}\cr
}
}

\author{ Matthieu Stigler}
\examples{

data(zeroyld)
data<-zeroyld

LRRbreak.test(zeroyld[,1],x=zeroyld[,2],model="CS",AdfLagChoice="AIC")
}

\keyword{ ts }

\references{
Hansen (1999) Testing for linearity, Journal of Economic Surveys, Volume 13, Number 5, December 1999 , pp. 551-576(26)
avalaible at: \url{http://www.ssc.wisc.edu/~bhansen/papers/cv.htm}

Lo and Zivot (2001) "Threshold Cointegration and Nonlinear Adjustment to the Law of One Price," Macroeconomic Dynamics, Cambridge University Press, vol. 5(4), pages 533-76, September.


}
