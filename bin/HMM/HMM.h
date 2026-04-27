#ifndef _HMM_H
#define _HMM_H

#include <vector>
#include <string>
#include <tuple>
#include <memory>

#include "matrix.h"
#include "distribution.h"

/**
 * @brief Multivariate Hidden Markov Model allowing for different distributions for all variables.
 * Variables are assumed to be independent. 
 * Calculations are performed in log-space for better numerical stability.
 */
class HMM
{
    public:
        template<class T>
        using matrix_ptr = std::shared_ptr<Matrix<T>>;
        template<class T>
        using const_matrix_ptr = std::shared_ptr<const Matrix<T>>;

        /**
         * @brief Constructs an uninitialized HMM, which should in the next step be initialized from a file via an input stream.
         */
        HMM(){};

        /**
         * @brief Constructs a new HMM object with given number of states and specified emission probabilities.
         * Transition matrix is uniformly initialized.
         * 
         * @param states 
         * @param emission (matrix states x variables)
         */
        explicit HMM(size_t states, Matrix<std::shared_ptr<DiscreteDistribution>>& emission, bool methylation = false);

        /**
         * @brief Returns the complete log likelihood of an observation in the current model.
         * 
         * @param observation
         * @param nObservation: optional, for two column input
         * @return double 
         */
        double log_likelihood(const_matrix_ptr<int> observation, const_matrix_ptr<int> nObservation = std::make_shared<Matrix<int>>(Matrix<int>()), std::pair<size_t, size_t> range = std::make_pair<size_t, size_t>(0,0)) const;

        /**
         * @brief Implementation of the Baum-Welch algorithm to fit the parameters of the HMM on multiple independent obsrvations.
         * 
         * @param training observations
         * @param nObservation two column input (trials, successes)
         * @param starting indices of observations
         * @param eps (stopping criterion for log-likelihood)
         * @param max_iteration 
         * @param update_init (change initial state distribution)
         */
        void train(const_matrix_ptr<int> observations, const_matrix_ptr<int> nObservations, std::vector<size_t>& index, double eps = 0.1, size_t maxIteration = 0, bool updateInit = true);

        /**
         * @brief Implementation of the Viterbi algorihtm to decode the observation sequence into a state sequence (find maximum likelihood path).
         * Returns path probability and state sequence.
         * 
         * @param observation
         * @param nObservation two column input (trials, successes)
         * @return std::pair<double, std::vector<int>> 
         */
        std::pair<double, std::vector<int>> viterbi_decoding(std::pair<size_t, size_t> range, const_matrix_ptr<int> observation, const_matrix_ptr<int> nObservation = std::make_shared<Matrix<int>>(Matrix<int>())) const;

        /**
         * @brief Posterior decoding to get sequence of most likely states (can lead to an invalid path if not all transitions are possible).
         * Return for each observation the most likely state and its posterior probability.
         * 
         * @param observation
         * @param nObservation two column input (trials, successes)
         * @return std::vector<std::pair<int, double>>
         */
        std::vector<std::pair<int, double>> posterior_decoding(std::pair<size_t, size_t> range, const_matrix_ptr<int> observation, const_matrix_ptr<int> nObservation = std::make_shared<Matrix<int>>(Matrix<int>())) const;

        /**
         * @brief calculate membership coefficients
         * 
         * @param observation
         * @param nObservation two column input (trials, successes)
         * @return Matrix<double>
         */
        std::vector<matrix_ptr<double>> membership_coefficients(std::vector<size_t>& index, const_matrix_ptr<int> observation, const_matrix_ptr<int> nObservation = std::make_shared<Matrix<int>>(Matrix<int>())) const;

        /**
         * @brief Returns the transition matrix.
         * 
         * @return matrix<double> 
         */
        Matrix<double> get_transition_matrix() const;

        /**
         * @brief Returns the emission matrix.
         * 
         * @return matrix<std::shared_ptr<DiscreteDistribution>> 
         */
        const Matrix<std::shared_ptr<DiscreteDistribution>>& get_emission_matrix() const {return emission;};

        /**
         * @brief Set a different transition matrix.
         * 
         * @param transition matrix
         */
        void set_transitions(const Matrix<double>&);
        
        /**
         * @brief Set a different initial state distribution.
         * 
         * @param initial state distribution
         */
        void set_initial(const std::vector<double>&);

        /**
         * @brief Returns if the HMM model excepts methylation data
         * 
         * @return methylation
         */
        bool has_methylation() {return methylation;};

        /**
         * @brief Outputs the initial state distribution, transition matrix and parameters of the distribution in the different states.
         * 
         * @return std::ostream& 
         */
        friend std::ostream& operator<<(std::ostream&, const HMM&);

        /**
         * @brief Outputs the initial state distribution, transition matrix and parameters of the distribution in the different states.
         * 
         * @return std::ostream& 
         */
        friend std::istream& operator>>(std::istream&, HMM&);

    protected:

        /**
         * @brief Calculates the forward step of the Baum-Welch algorithm and return the log-likelihood of the observation
         * 
         * @param logEmission
         * @param logAlpha
         * 
         * @return log-likelihood 
         */
        void forward(const_matrix_ptr<double>, matrix_ptr<double>, double&) const;

        /**
         * @brief Calculates the backward step of the Baum-Welch algorihtm.
         * 
         * @param logEmission
         * @param logBeta
         */
        void backward(const_matrix_ptr<double>, matrix_ptr<double>) const;

        /**
         * @brief Performs the E-step of the Baum-Welch algorithm.
         * 
         * @param logEmission
         * @param logAlpha
         * @param logBeta
         * @param logGamma
         */
        void E_step(const_matrix_ptr<double>, matrix_ptr<double>, matrix_ptr<double>, matrix_ptr<double>) const;

        /**
         * @brief Performs M step for multiple observation
         * 
         * @param observations
         * @param nObservation two column input (trials, successes)
         * @param logEmission
         * @param logGamma
         * @param logBeta
         * @param update initial state distribution
         */
        virtual void M_step(const_matrix_ptr<int>, const_matrix_ptr<int>, std::vector<matrix_ptr<double>>&,  std::vector<matrix_ptr<double>>&,  std::vector<matrix_ptr<double>>&, bool);
    
        /**
         * @brief Calculates the emission probability in state s at time t for the given observation.
         * 
         * @param log emission probability
         * @param observation
         * @param nObservation
         */
        virtual void calculate_log_emission(matrix_ptr<double>, const_matrix_ptr<int>, const_matrix_ptr<int>, size_t i = 0) const;

        /**
         * @brief Uniform initialization of the transition matrix.
         */
        virtual void init_transitions();

        /**
         * @brief Uniform initialization of the initial state distribution.
         */
        virtual void init_initial();
        
        /**
         * @brief Normalize initial state distribution and transition matrix to sum up to one 
         * (required after setting an entry to zero).
         */
        void normalize_transition_matrix();

        size_t N;
        size_t m;
        std::vector<std::string> marker;
        bool methylation;
        Matrix<double> logA;
        std::vector<double> logPi;
        Matrix<std::shared_ptr<DiscreteDistribution>> emission;
};

#endif