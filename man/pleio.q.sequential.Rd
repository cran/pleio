\name{pleio.q.sequential}
\alias{pleio.q.sequential}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform sequential tests of pleiotropy
}
\description{
  Perform sequential tests of pleiotropy to determine the number of
  traits, and which traits, are associatd with a vector of genotypes.
}
\usage{
pleio.q.sequential(obj.fit, pval.threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.fit}{ result of pleio.q.fit }
  \item{pval.threshold}{p-value for rejecting the null hypothesis of the
    specified number of non-zero beta's.
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  Perform sequential tests of pleiotropy, starting at the usual
  multivarite null hypothesis that all beta's = 0. If this test
  rejects becuase the p-value < pval.threshold, then allow one beta to
  be non-zero in order to test whether the remaining beta's = 0.
  If the test of one non-zero beta rejects, then allow two
  non-zero beta's, and continue this sequential testing until the
  p-value for a test is greater than the specified pval.threshold.
}
\value{
%%  ~Describe the value returned
A list containing:
\item{pval }{p-value of the final test from the sequential testing}
\item{index.beta}{index of trait(s) in y that have non-zero beta's.
  These indices indicate which traits are associated with the genotype,
  accounting for the correlations among the traits.}
}
\references{
  %% ~put references to the literature/web site here ~
   Schaid DJ, Tong X, Larrabee B, Kennedy RB, Poland GA, Sinnwell JP.
  Statistical Methods for Testing Genetic Pleiotropy. To appear in Genetics, 2016.
}
\author{
  %%  ~~who you are~~
  Dan Schaid and Jason Sinnwell
}
\examples{
data(pleio.qsim)
fit <- pleio.q.fit(y, geno)
test.seq <- pleio.q.sequential(fit, pval.threshold=.05)
test.seq
}

