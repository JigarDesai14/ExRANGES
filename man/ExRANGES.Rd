\name{ExRANGES}
\alias{calc.slope}
\alias{sample.pval.calc}
\title{RANGES Functions}
\usage{
calc.slopes(time.series, cycle=F, last.time.step)

sample.pval.calc(slopes, sample.size=10000)
}
\description{
The following two functions are used to calcuates the RANGES values to be used to calculate ExRANGES.

}
\arguments{
calc.slopes()

time.series - A matrix with rows as genes and columns as sample names. Samples must be labeled as XX_samplename, XX=Numeric number for time.

cycle - If data is cyclical then define the last time step form the last column to the first.

last.time.step - Time step from last column to the fist column

#########################################################################

sample.pval.calc()

slopes - Output of calc.slopes(). Should be transposed pvalues of calculates slopes between time points.

sample.size - How many time should the slopes be sampled for each gene to calculate a pvalue.

}
\examples{

}
