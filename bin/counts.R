#!/usr/bin/env Rscript

suppressPackageStartupMessages(suppressWarnings(library(GenomicRanges)))
suppressPackageStartupMessages(suppressWarnings(library(BiocParallel)))
suppressPackageStartupMessages(suppressWarnings(library(Rsamtools)))
suppressPackageStartupMessages(suppressWarnings(library(bamsignals)))
suppressPackageStartupMessages(suppressWarnings(library(preprocessCore)))
suppressPackageStartupMessages(suppressWarnings(library(matrixStats)))

script.dir <- dirname(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
library(Rcpp)
sourceCpp(file.path(script.dir, "counts.cpp"))

# Reading all the arguments
args = commandArgs(trailingOnly=TRUE)
# print(args) # Uncomment to check for provided arguments
tab_files <- read.table(args[1],header=F)
ut.regions <- read.table(args[2],header = F)
bin_size <- as.integer(args[3])
shift_bp <- as.integer(args[4])
output_path <- as.character(args[5])
cores <- as.integer(args[6])
paired_end <- as.character(args[7])

print("All inputs are read...")

tab_files$V2 <- as.character(tab_files$V2)

# Generating refined regions 
regions <- GRanges(ut.regions[[1]], IRanges(start=ut.regions[[2]]+1, end=ut.regions[[3]]))
starts <- start(regions)-1
ends <- end(regions)
newstarts <- bin_size*(ceiling(starts/bin_size))
newends <- bin_size*(floor(ends/bin_size))
valid <- newends > newstarts
regions <- regions[valid]
newstarts <- newstarts[valid]
newends <- newends[valid]
start(regions) <- newstarts+1
end(regions) <- newends


# List of BAM files
bam_files <- c(tab_files$V2)

# Initialize BiocParallel
register(MulticoreParam(workers = cores))  # Adjust the number of workers as needed

# Process each BAM file in parallel
results <- bplapply(bam_files, function(bam_path) {
  bprof <- bamProfile(bampath = bam_path, regions, binsize = bin_size, shift = shift_bp, paired.end=paired_end)
  unlist(as.list(bprof))
})


# Combine the results into a list
counts_list <- as.data.frame(results)
colnames(counts_list) <- c(as.character(tab_files$V1))

# Find unique column name patterns
unique_patterns <- unique(sub("\\.\\d+$", "", names(counts_list)))

# Use lapply to group and sum rows based on patterns
summed_data_list <- lapply(unique_patterns, function(pattern) {
  print("Starting normalization...")
    if(ncol(counts_list[grep(pattern, names(counts_list))]) > 1) {
        input_df <- counts_list[grep(pattern, names(counts_list))]
        #outlier removal
        for(i in 1:ncol(input_df)){
          input_df <- subset(input_df, input_df[,i]<=quantile(input_df[,i],prob=.999))
        }
        normalized_counts <- quantileNormalization(input_df)
        normalized_counts <- as.matrix(normalized_counts)
        print("calculating rowmedians...")
        summed_cols <- floor(rowMedians(normalized_counts))
    } else {
        summed_cols <- counts_list[grep(pattern, names(counts_list))]
    }
    data.frame(summed_cols)
}
)


# Combine the list of dataframes into a single dataframe
summed_data <- do.call(cbind, summed_data_list)
colnames(summed_data) <- unique_patterns

refined_regions <- data.frame(chr= as.character(seqnames(regions)), start=start(regions)-1, end=end(regions))

# write the count matrix and refined regions to output file 
print("writing tables...")
write.table(summed_data,file=paste(output_path,"counts.txt",sep="_"),sep="\t",quote=FALSE,row.names=F)
write.table(refined_regions,file=paste(output_path,"refine_chr_regions.bed",sep="_"),sep="\t",quote=FALSE,row.names=F, col.names=FALSE)