#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <algorithm>

#include <omp.h>
#include "boost/program_options.hpp"
#include "TROOT.h"
#include "TError.h"

#include "adjustableDurationHMM.h"
#include "matrix_reader.h"

namespace bpo = boost::program_options;

bool train, adjustTopology;
size_t maxIteration, iterationsPretrain, numAdjustments, threads;
double epsilon;
std::string modelInput, modelOutput, viterbiPath, posteriorDecoding, countMatrix, methylationMatrix, regions;

bool parseArguments(int argc, char* argv[])
{
    bpo::variables_map vm;
    bpo::options_description desc;

    desc.add_options()
        ("help,h", "Output help message.")
        ("modelInput,m", bpo::value<std::string>(&modelInput)->required(), "Serialized HMM with starting parameters.")
        ("regions,r", bpo::value<std::string>(&regions)->default_value(""), "Start indices for independent regions (starting at zero and whitespace separated), if not set no splitting is performed.")
        ("countMatrix,c", bpo::value<std::string>(&countMatrix)->required(), "Count matrix for training or decoding.")
        ("methylationMatrix,x", bpo::value<std::string>(&methylationMatrix)->default_value(""), "Optional matrix for DNA methylation, only required if model has set methylation to true.")
        ("threads,p", bpo::value<size_t>(&threads)->default_value(1), "Number of threads for training the HMM model.");

    try
    {
        bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),vm);
        if (vm.count("help")) 
        {
            std::cout << "Command line options: " << std::endl;
            std::cout << desc << "\n";
            return false;
        }
        bpo::notify(vm);
    }
    catch(bpo::error& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        desc.print(std::cerr);
        return false;
    }
    return true;
}

int main(int argc, char* argv[])
{
    ROOT::EnableThreadSafety();
    gErrorIgnoreLevel = kFatal;
    if(!parseArguments(argc, argv))
    {
        return -1;
    }
    
    omp_set_num_threads(threads);

    std::ifstream counts (countMatrix, std::ifstream::in);
    if (!counts) 
    {
        std::cerr << "Cannot open input file: " + countMatrix << std::endl;
        return -1;
    }
    
    Reader r;
    HMM::const_matrix_ptr<int> observation = std::make_shared<Matrix<int>>(r.parse_matrix(counts));
    counts.close();
    if (!r.get_message().empty())
    {
        std::cerr << r.get_message() << std::endl;
        return -1;
    }

    std::ifstream startModel (modelInput, std::ifstream::in);
    if (!startModel)
    {
        std::cerr << "Cannot open input file: " + modelInput << std::endl;
        return  -1;
    } 

    try
    {
        AdjustableDurationHMM hmm;
        startModel >> hmm;
        startModel.close();
        
        HMM::const_matrix_ptr<int> nObservation = std::make_shared<Matrix<int>>(Matrix<int>());
        if (hmm.has_methylation())
        {
            std::ifstream methylation (methylationMatrix, std::ifstream::in);
            if (!methylation) 
            {
                std::cerr << "Cannot open input file: " + methylationMatrix << std::endl;
                return -1;
            }
                
            nObservation = std::make_shared<Matrix<int>> (r.parse_methylation_matrix(methylation));
            methylation.close();
            if (!r.get_message().empty())
            {
                std::cerr << r.get_message() << std::endl;
                return -1;
            }

            if (nObservation->nrows() != observation->nrows())
            {
                std::cerr << "Number of rows of count matrix and methylation matrix must be the same." << std::endl;
                return -1;
            }
        }

        std::vector<size_t> startIndex {0};
        if (!regions.empty())
        {
            std::ifstream index (regions, std::ifstream::in);
            if (!index) 
            {
                std::cerr << "Cannot open input file: " + regions << std::endl;
                return -1;
            }
            startIndex = r.parse_regions(index, observation->nrows());
            index.close();
            if (!r.get_message().empty())
            {
                std::cerr << r.get_message() << std::endl;
                return -1;
            }
        }
        
        double logL = 0;
        for (size_t k = 0; k < startIndex.size(); ++k)
        {
            size_t end = k < startIndex.size()-1 ? startIndex[k+1] : observation->nrows();
            logL += hmm.log_likelihood(observation, nObservation, std::make_pair(startIndex[k], end));
        }

        std::cout << logL;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return -1;
    }
    return 0;
}