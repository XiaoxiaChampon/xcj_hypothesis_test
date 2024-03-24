# install.packages("optparse")
# Load the optparse library
library(optparse)

# Define options
option_list <- list(
    make_option(c("-n", "--jobid"), type="integer", default=123,
                help="Job Index", metavar="JOBID"),
    make_option(c("-a", "--numcpus"), type="integer", default=32,
                help="Num CPUSs", metavar="NUMCPUS")
)

# Create parser and parse options
parser <- OptionParser(option_list=option_list)
options <- parse_args(parser)

# Use the options
cat("Job Idx:", options$jobid, "\n")
cat("Num CPUs:", options$numcpus, "\n")
