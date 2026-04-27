#ifndef _ADJUSTABLE_DURATION_HMM_H
#define _ADJUSTABLE_DURATION_HMM_H

#include <vector>
#include <tuple>

#include "HMM.h"

/**
 * @brief Multivariate Hidden Markov Model with duration modelling allowing for different distributions for all variables.
 * Each state can have a different sub-HMM topology (different number of replicate states with same emission probabilities and self-transition probabilities)
 * Variables are assumed to be independent. 
 * Calculations are performed in log-space for better numerical stability.
 */
class AdjustableDurationHMM: public HMM
{
    public:

    /**
     * @brief Constructs an uninitialized HMM, which should in the next step be initialized from a file via an input stream.
     */
    explicit AdjustableDurationHMM() {};

    /**
     * @brief Post-process state sequence to assign same label to all states belonging to the same sub-HMM
     **/
    void process_state_sequence(std::vector<int>& stateSequence) const;

    /**
     * @brief Calculate state sequence and adjust HMM topolgy to better model state durations
     * 
     * @param observations
     * @param nObservation two column input (trials, successes)
     * @param starting indices of observations
     **/
    void adjust_topology(const_matrix_ptr<int>, const_matrix_ptr<int>, std::vector<size_t>&);

    /**
     * @brief Outputs the initial state distribution, transition matrix, states in sub-HMMs and parameters of the distribution in the different states.
     * 
     * @return std::ostream& 
     */
    friend std::ostream& operator<<(std::ostream&, const AdjustableDurationHMM&);

    /**
     * @brief Reads HMM from input file.
     * 
     * @return std::ostream& 
     */
    friend std::istream& operator>>(std::istream&, AdjustableDurationHMM&);

    private:

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
    virtual void M_step(const_matrix_ptr<int>, const_matrix_ptr<int>, std::vector<matrix_ptr<double>>&,  std::vector<matrix_ptr<double>>&,  std::vector<matrix_ptr<double>>&, bool) override;

    /**
     * @brief Updates the transition matrix of states i
     * 
     * @param index
     * @param logEmission
     * @param logGamma
     * @param logBeta
     * @param update initial state distribution
     * @param logA
     */
    void update_transitions(size_t, std::vector<HMM::matrix_ptr<double>>&, std::vector<HMM::matrix_ptr<double>>&, std::vector<HMM::matrix_ptr<double>>&, bool, Matrix<double>&);

    /**
     * @brief Updates the emission probabilities
     * 
     * @param observations
     * @param nObservation two column input (trials, successes)
     * @param logGamma
     */
    void update_emission(size_t, HMM::const_matrix_ptr<int>, HMM::const_matrix_ptr<int>, std::vector<HMM::matrix_ptr<double>>&);

    /**
     * @brief Calculates the emission probability in state s at time t for the given observation.
     * 
     * @param log emission probability
     * @param observation
     * @param nObservation
     */
    virtual void calculate_log_emission(HMM::matrix_ptr<double>, HMM::const_matrix_ptr<int>, HMM::const_matrix_ptr<int>, size_t i = 0) const override;

    /**
     * @brief Initialization of the transition matrix.
     */
    virtual void init_transitions() override;

    /**
     * @brief Initialization of the initial state distribution.
     */
    virtual void init_initial() override;

    /**
     * @brief Initialization of the transition matrix.
     * 
     * @param self transition probabilities
     */
    void init_transitions(const std::vector<double>&);

    /**
     * @brief Calculates the segment lengths of a state decoding with n different states
     * 
     * @param n states
     * @param decoding
    */
    std::vector<std::vector<size_t>> calculate_segement_lengths(size_t, const std::vector<std::vector<int>>&);

    bool updateEmission;
    bool updateTransition;
    std::vector<std::vector<size_t>> stateIndices;
    std::vector<size_t> stateAssignment;
};

#endif
