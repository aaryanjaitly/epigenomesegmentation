#include <Rcpp.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <ctime> 

// need to be compiled like
// g++ -o [outfile] [file] -fopenmp

using namespace Rcpp;

// pretty basic solution for now, in future may look into Boost MultiArray!
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> data, int threads = 1) {
    // this assumes that all inner vectors have the same size and
    // allocates space for the complete result in advance
    std::vector<std::vector<double>> result(data[0].size(), std::vector<double>(data.size()));
    #pragma omp parallel num_threads(threads)
    {
        #pragma omp for
        for (std::vector<double>::size_type i = 0; i < data[0].size(); i++)
        { 
            for (std::vector<double>::size_type j = 0; j < data.size(); j++) 
            {
                result[i][j] = data[j][i];
            }
        }
    }
    return result;
}

std::vector<double> calculate_rowMeans(std::vector<std::vector<double>> vec, int threads = 1)
{
    std::vector<double> res;

    #pragma omp parallel num_threads(threads)
    {
        #pragma omp for
        for(int i = 0; i < vec.size(); i++)
        {
            double value = std::accumulate(vec[i].begin(), vec[i].end(), 0.0) / (double)vec[i].size();
            res.push_back(value);
                
        }
    }

    return res;
}

// [[Rcpp::export]]
IntegerVector quantileNormalization(DataFrame df, int threads = 1)
{
    auto start = std::chrono::system_clock::now();
    //CharacterVector column_names = df.names();
    int num_cols = df.length();

    if(num_cols < 2)
    {
        stop("DataFrame must have more than 1 column!");
    }

    Rcout << "Starting the cpp script..." << std::endl;

    // matrix to store the ranks of the column values
    // initial_ranks: first = original, second = rank
    std::vector<std::vector<std::pair<double,int>>> initial_ranks;
    std::vector<std::vector<std::pair<double,int>>> sorted_ranks;
    std::vector<std::vector<double>> sorted_values;
    std::vector<std::vector<std::pair<int, int>>> mean_ranks;
    std::vector<int> output;

    // enable parallel computing
    #pragma omp parallel num_threads(threads)
    {
        //iterate over all columns of the dataframe
        #pragma omp for
        for(int i = 0; i < num_cols; i++)
        {
            // convert each column to c++ vector and rank it
            // STEP 1
            Rcout << "Assigning first ranks to column " << i << "..." << std::endl;
            std::vector<double> column = as<std::vector<double>>(df[i]);
            std::vector<std::pair<double,int>> ranks;
            std::vector<double> values;
            for(int j = 1; j <= column.size(); j++)
            {
                ranks.push_back(std::pair<double,int>(column[j-1],j));
            }
            initial_ranks.push_back(ranks);
            
            // sort the vector of pairs
            // STEP 2
            Rcout << "Sorting column " << i << "..." << std::endl;
            std::vector<double> vals = column;
            std::sort(vals.begin(), vals.end());
            std::sort(ranks.begin(), ranks.end());
            sorted_ranks.push_back(ranks);
            sorted_values.push_back(vals);
        }

        // compute row means
        // STEP 3
        Rcout << "Computing row means..." << std::endl;
        std::vector<std::vector<double>> transposed_values = transpose(sorted_values);

        std::vector<double> rowmeans = calculate_rowMeans(transposed_values, threads);

        // assign the rowmeans with the sorted ranks
        // STEP 4
        #pragma omp for
        for(int i = 0; i < sorted_ranks.size(); i++)
        {
            std::vector<std::pair<int, int>> new_col;
            for(int j = 0; j < sorted_ranks[i].size(); j++)
            {
                std::pair<int, int> p{sorted_ranks[i][j].second, (int)rowmeans[j]};
                new_col.push_back(p);
            }
            std::sort(new_col.begin(), new_col.end());

            for(auto it : new_col)
            {
                output.push_back(it.second);
            }
            mean_ranks.push_back(new_col);
        }
    }

    Rcout << "Converting cpp objects to R objects..." << std::endl;
    IntegerVector m = Rcpp::wrap(output);
    m.attr("dim") = Dimension(initial_ranks[0].size(), num_cols);

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    Rcout << "finished computation at " << std::ctime(&end_time) << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

    return m;
}