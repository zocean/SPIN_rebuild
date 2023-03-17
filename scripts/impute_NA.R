#!Rscript

library(argparse)
library(missForest)

parser <- ArgumentParser(description="")
parser$add_argument('input',help='input file')
parser$add_argument('--output',dest='output',help='output file name after imputation')
parser$add_argument('--na_string',dest='na_string',help='string used to represent NA value')
parser$add_argument('--seed',type='integer',dest='seed',help='random seed number')

# parse argument
args <- parser$parse_args()
input = args$input
output = args$output
na_string = as.character(args$na_string)
seed_number = args$seed

# set random seed
set.seed(seed_number)

# load
data <- read.table(file = input, header = T, sep = '\t', na.strings = na_string, stringsAsFactors = F)

print(summary(data))

# impute
data.imp <- missForest(data, verbose = TRUE)

# check QC
print(data.imp$OOBerror)

# report to output
result = data.imp$ximp
colnames(result) <- colnames(data)
write.table(result, file = output, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
