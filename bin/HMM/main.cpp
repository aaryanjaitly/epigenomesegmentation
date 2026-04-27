#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <omp.h>

#include "TROOT.h"
#include "TError.h"

#include "boost/program_options.hpp"

#include "HMM.h"
#include "matrix_reader.h"

namespace bpo = boost::program_options;

bool train;
size_t maxIteration, threads;
double epsilon;
std::string modelInput, modelOutput, viterbiPath, posteriorDecoding, countMatrix, methylationMatrix, regions;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()
        ("help,h", "Output help message.")
        ("modelInput,m", bpo::value<std::string>(&modelInput)->required(), "Serialized HMM with starting parameters.")
        ("countMatrix,c", bpo::value<std::string>(&countMatrix)->default_value(""), "Count matrix for training or decoding.")
        ("methylationMatrix,x", bpo::value<std::string>(&methylationMatrix)->default_value(""), "Optional matrix for DNA methylation, only required if model has set methylation to true.")
        ("regions,r", bpo::value<std::string>(&regions)->default_value(""), "Start indices for independent regions (starting at zero and whitespace separated), if not set no splitting is performed.")
		("train,t", bpo::bool_switch(&train), "Should the HMM model be trained first?")
		("maxIteration,i", bpo::value<size_t>(&maxIteration)->default_value(300), "Maximum number of iteration during training.")
		("epsilon,e", bpo::value<double>(&epsilon)->default_value(0.1), "Convergence of likelihoods as termination criterion during training.")
		("modelOutput,o", bpo::value<std::string>(&modelOutput)->default_value(""), "Optional output file for parameters of new HMM (only required if HMM is trained first).")
        ("viterbiPath,v", bpo::value<std::string>(&viterbiPath)->default_value(""), "Optional output file for state sequence using Viterbi decoding.")
        ("posteriorDecoding,d", bpo::value<std::string>(&posteriorDecoding)->default_value(""), "Optional output file for state sequence using posterior decoding.")
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

    std::ifstream startModel (modelInput, std::ifstream::in);
    if (!startModel)
    {
        std::cerr << "Cannot open input file: " + modelInput << std::endl;
        return  -1;
    } 

    try
    {
        HMM model;
        startModel >> model;
        startModel.close();

        Reader r;
        HMM::const_matrix_ptr<int> observation = std::make_shared<Matrix<int>>(Matrix<int>());
        if (countMatrix != "")
        {
            std::ifstream counts (countMatrix, std::ifstream::in);
            if (!counts) 
            {
                std::cerr << "Cannot open input file: " + countMatrix << std::endl;
                return -1;
            }
            observation = std::make_shared<Matrix<int>>(r.parse_matrix(counts));
            counts.close();
            if (!r.get_message().empty())
            {
                std::cerr << r.get_message() << std::endl;
                return -1;
            }
        }
        
        HMM::const_matrix_ptr<int> nObservation = std::make_shared<Matrix<int>>(Matrix<int>());
        if (model.has_methylation())
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

            if (observation->nrows() > 0 && nObservation->nrows() != observation->nrows())
            {
                std::cerr << "Number of rows of count matrix and methylation matrix must be the same." << std::endl;
                return -1;
            }
        }

        size_t rows = observation->nrows() > 0 ? observation->nrows() : nObservation->nrows();

        std::vector<size_t> startIndex {0};
        if (!regions.empty())
        {
            std::ifstream index (regions, std::ifstream::in);
            if (!index) 
            {
                std::cerr << "Cannot open input file: " + regions << std::endl;
                return -1;
            }
            startIndex = r.parse_regions(index, rows);
            index.close();
            if (!r.get_message().empty())
            {
                std::cerr << r.get_message() << std::endl;
                return -1;
            }
        }

        if (train)
        {
            model.train(observation, nObservation, startIndex, epsilon, maxIteration, false);

            if (!modelOutput.empty())
            {
                std::ofstream finalModel (modelOutput, std::ifstream::out);
                if (!finalModel) 
                {
                    std::cerr << "Cannot open output file: " << modelOutput << std::endl;
                    return -1;
                } 
                finalModel << model;
                finalModel.close();
            }
        }   

        if (!viterbiPath.empty())
        {
            std::ofstream outputDecoding (viterbiPath, std::ifstream::out);
            if (!outputDecoding)
            {
                std::cerr << "Cannot open output file: " << viterbiPath << std::endl;
                return -1;
            }

            for (size_t k = 0; k < startIndex.size(); ++k)
            {
                size_t end = k < startIndex.size()-1 ? startIndex[k+1] : rows;
                std::pair<double, std::vector<int>> decoding = model.viterbi_decoding(std::make_pair(startIndex[k], end), observation, nObservation);
                for (const auto s : decoding.second)
                {
                    outputDecoding << s << '\n';
                }
            }
            outputDecoding.close();
        }

        if (!posteriorDecoding.empty())
        {
            std::ofstream outputDecoding (posteriorDecoding, std::ifstream::out);
            if (!outputDecoding)
            {
                std::cerr << "Cannot open output file: " << posteriorDecoding << std::endl;
                return -1;
            }
            for (size_t k = 0; k < startIndex.size(); ++k)
            {
                size_t end = k < startIndex.size()-1 ? startIndex[k+1] : rows;
                std::vector<std::pair<int, double>> decoding = model.posterior_decoding(std::make_pair(startIndex[k], end), observation, nObservation);
                for (const auto s : decoding)
                {
                    outputDecoding << s.first << '\t' << s.second << '\n';
                }
            }
            outputDecoding.close();
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return -1;
    }
    return 0;
}