#include "adjustableDurationHMM.h"

#include <thread>
#include <algorithm>
#include "log_prob.h"

#include <fstream>


void AdjustableDurationHMM::calculate_log_emission(HMM::matrix_ptr<double> logEmission, HMM::const_matrix_ptr<int> observation, HMM::const_matrix_ptr<int> nObservation, size_t start) const
{
    size_t states = stateIndices.size();
    size_t T = logEmission->nrows();
    #pragma omp parallel for collapse(2) schedule(static)
    for (size_t t = 0; t < T; ++t)
    {
        for (size_t i = 0; i < N; ++i)
        {
            size_t disIndex = stateAssignment[i];
            double p = 0.0;

            if (m > 0)
            {
                auto it = observation->row_begin(start+t);
                //marks are assumed to be independent
                for (size_t k = 0; k < m; ++k)
                {
                    p = lp::log_mul(p, emission(disIndex, k)->log_pmf(*it));
                    ++it;
                }
            }
            if (methylation)
            {
                std::shared_ptr<TwoValueDiscreteDistribution> dis = std::dynamic_pointer_cast<TwoValueDiscreteDistribution>(emission(disIndex, m));
                p = lp::log_mul(p, dis->log_pmf((*nObservation)(start+t, 0), (*nObservation)(start+t, 1)));
            }
            (*logEmission)(t, i) = p;
        }
    }
}

void AdjustableDurationHMM::update_transitions(size_t index, std::vector<HMM::matrix_ptr<double>>& logEmission, std::vector<HMM::matrix_ptr<double>>& logGamma, std::vector<HMM::matrix_ptr<double>>& logBeta, bool updateInit, Matrix<double>& new_logA)
{
    const std::vector<size_t>& substates = stateIndices[index];
    size_t numStates = substates.size();
    size_t startIndex = substates[0];
    size_t endIndex = substates[numStates-1];
    
    size_t obs = logGamma.size();
    double lowerProbBound = lp::ext_log(pow(10.0, -10.0));
    double xi_tij = 0;

    //update initial state distribution if required (only for first state of sub-HMM)
    if (updateInit)
    {
        logPi[startIndex] = (*logGamma[0])(0, startIndex);
        for (size_t k = 1; k < obs; ++k)
        {
            logPi[startIndex] = lp::log_add(logPi[startIndex], (*logGamma[k])(0, startIndex));
        }
        logPi[startIndex] -= log(obs);
    }

    // if sub-HMM contains a single state update in the same way as ergodic HMMs
    if (numStates == 1)
    {
        double denom = NAN;
        for (size_t k = 0; k < obs; ++k)
        {
            size_t T = logGamma[k]->nrows();
            for (size_t t = 0; t < T-1; ++t)
            {
                double gamma = (*logGamma[k])(t, startIndex);
                denom = lp::log_add(denom, gamma);
            }
        }

        for (size_t j = 0; j < N; ++j)
        {
            double num = NAN;
            for (size_t k = 0; k < obs; ++k)
            {
                size_t T = logGamma[k]->nrows();
                for (size_t t = 0; t < T-1; ++t)
                {
                    xi_tij = (*logGamma[k])(t, startIndex) + logA(startIndex, j) + (*logEmission[k])(t+1, j) + (*logBeta[k])(t+1, j) - (*logBeta[k])(t, startIndex);
                    num = lp::log_add(num, xi_tij);
                }
            }
            new_logA(startIndex, j) = lp::log_mul(num, -denom);
        }
    }
    else
    {
        // if sub-HMM has > 1 states, update transitions such that all states have the same self transition probability

        double denom = NAN;
        for (size_t k = 0; k < obs; ++k)
        {
            size_t T = logGamma[k]->nrows();
            for (size_t t = 0; t < T-1; ++t)
            {
                for (size_t s : substates)
                {
                    denom = lp::log_add(denom, (*logGamma[k])(t, s));
                }
            }
        }

        double num = NAN;
        for (size_t k = 0; k < obs; ++k)
        {
            size_t T = logGamma[k]->nrows();
            for (size_t t = 0; t < T-1; ++t)
            {
                for (size_t s : substates)
                {
                    xi_tij = (*logGamma[k])(t, s) + logA(s, s) + (*logEmission[k])(t+1, s) + (*logBeta[k])(t+1, s) - (*logBeta[k])(t, s);
                    num = lp::log_add(num, xi_tij);
                }
            }
        }

        for (size_t i = 0; i < numStates-1; ++i)
        {
            new_logA(substates[i], substates[i]) = lp::log_mul(num, -denom);
            new_logA(substates[i], substates[i+1]) = lp::log1x(-lp::ext_exp(new_logA(substates[i], substates[i])));
        }
        new_logA(endIndex, endIndex) = lp::log_mul(num, -denom);


        // update outgoing transitions of last state 
        denom = NAN;
        for (size_t k = 0; k < obs; ++k)
        {
            size_t T = logGamma[k]->nrows();
            for (size_t t = 0; t < T-1; ++t)
            {
                denom = lp::log_add(denom, (*logGamma[k])(t, endIndex));
            }
        }
        
        double sum = NAN;
        for (size_t j = 0; j < N; ++j)
        {
            if (j != endIndex)
            {
                double num = NAN;
                for (size_t k = 0; k < obs; ++k)
                {
                    size_t T = logGamma[k]->nrows();
                    for (size_t t = 0; t < T-1; ++t)
                    {
                        xi_tij = (*logGamma[k])(t, endIndex) + logA(endIndex, j) + (*logEmission[k])(t+1, j) + (*logBeta[k])(t+1, j) - (*logBeta[k])(t, endIndex);
                        num = lp::log_add(num, xi_tij);
                    }
                }
                new_logA(endIndex, j) = lp::log_mul(num, -denom);

                sum = lp::log_add(sum, new_logA(endIndex, j));
            }
        }


        // normalize outgoing transitions to one 
        double ratio = lp::log_mul(new_logA(substates[0], substates[1]), -sum);
        for (size_t j = 0; j < N; ++j)
        {
            if (j != endIndex)
            {
                new_logA(endIndex, j) = lp::log_mul(new_logA(endIndex, j), ratio);  
            }
        }
    }
}

void AdjustableDurationHMM::update_emission(size_t state, HMM::const_matrix_ptr<int> observation, HMM::const_matrix_ptr<int> nObservation, std::vector<HMM::matrix_ptr<double>>& logGamma)
{
    size_t obs = logGamma.size();
    size_t totalT = observation->nrows() > 0 ? observation->nrows() : nObservation->nrows();
    std::vector<double> completeGamma (totalT, 0.0);

    size_t curr = 0;
    for (size_t k = 0; k < obs; ++k)
    {
        size_t T = logGamma[k]->nrows();
        for (size_t t = 0; t < T; ++t)
        {
            for (auto s : stateIndices[state])
            {
                completeGamma[curr] += lp::ext_exp((*logGamma[k])(t, s));
            }
            ++curr;
        }
    }

    //only update parameters if membership coefficients are not all zero
    if (std::accumulate(completeGamma.begin(), completeGamma.end(), 0.0) == 0.0)
    {
        #pragma omp critical
        std::clog << "Membership coefficient are zero for state " << state << '\n';
    }
    else
    {
        for (size_t k = 0; k < m; ++k)
        {
            emission(state, k)->update(completeGamma.begin(), completeGamma.end(), observation->col_begin(k), observation->col_end(k));
        }
        if (methylation)
        {
            std::shared_ptr<TwoValueDiscreteDistribution> dis = std::dynamic_pointer_cast<TwoValueDiscreteDistribution>(emission(state, m));
            dis->update_methylation(completeGamma.begin(), completeGamma.end(), nObservation->col_begin(0), nObservation->col_begin(1));
        }
    }
}

void AdjustableDurationHMM::M_step(HMM::const_matrix_ptr<int> observation, HMM::const_matrix_ptr<int> nObservation, std::vector<matrix_ptr<double>>& logEmission, std::vector<matrix_ptr<double>>& logGamma, std::vector<matrix_ptr<double>>& logBeta, bool updateInit)
{
    size_t states = stateIndices.size();
    Matrix<double> new_logA = Matrix<double> (N, N, NAN);
    
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < states; ++i)
    {
        update_transitions(i, logEmission, logGamma, logBeta, updateInit, new_logA);
        update_emission(i, observation, nObservation, logGamma);
    }
    logA = new_logA;
}

void AdjustableDurationHMM::init_initial()
{
    logPi = std::vector<double> (N, NAN);
    double logP = lp::ext_log(1.0 / stateIndices.size());
    for (auto s : stateIndices)
    {
        logPi[s[0]] = logP;
    }
}

void AdjustableDurationHMM::init_transitions(const std::vector<double>& selfP)
{
    size_t states = stateIndices.size();
    
    logA = Matrix<double> (N, N, NAN);

    // outgoing transition probabilities are initialized uniformly (considering fixed self transition probability)
    for (size_t i = 0; i < states; ++i)
    {
        double logSelfP = lp::ext_log(selfP[i]);
        double logNextP = lp::ext_log((1.0 - selfP[i]));
        
        // set transition probabilities for inner states of  the sub-HMM
        const std::vector<size_t>& subHMM = stateIndices[i];
        for (size_t j = 0; j < subHMM.size()-1; ++j)
        {
            logA(subHMM[j], subHMM[j+1]) = logNextP;
            logA(subHMM[j], subHMM[j]) = logSelfP;
        }
        
        // set transition probabilities for last state of the sub-HMM
        double logSwicthP = lp::ext_log((1.0 - selfP[i]) * (1.0 / (states-1)));
        size_t endState = subHMM[subHMM.size()-1];
        logA(endState, endState) = logSelfP;
        for (size_t j = 0; j < states; ++j)
        {
            if (i != j)
            {
                size_t nextState = stateIndices[j][0];
                logA(endState, nextState) = logSwicthP;
            }
        } 
    }    
}

void AdjustableDurationHMM::init_transitions()
{
    init_transitions(std::vector<double>(stateIndices.size(), 0.8));   
}

void AdjustableDurationHMM::process_state_sequence(std::vector<int>& stateSequence) const
{
    // all states of the sub-HMM have the same state label
    for (size_t i = 0; i < stateSequence.size(); ++i)
    {
        stateSequence[i] = stateAssignment[stateSequence[i]-1]+1;
    }
}

std::vector<std::vector<size_t>> AdjustableDurationHMM::calculate_segement_lengths(size_t states, const std::vector<std::vector<int>>& decoding)
{
    std::vector<std::vector<size_t>> segmentLengths (states);
    for (size_t i = 0; i < decoding.size(); ++i)
    {
        size_t length = 0;
        size_t index = decoding[i][0];
        for (size_t j = 0; j < decoding[i].size(); ++j)
        {
            if (decoding[i][j] == index)
            {
                length += 1;
            }
            else
            {
                segmentLengths[index-1].push_back(length);
                index = decoding[i][j];
                length = 1;
            }
        }
        segmentLengths[index-1].push_back(length);
    }
    return std::move(segmentLengths);
}

void AdjustableDurationHMM::adjust_topology(const_matrix_ptr<int> observation, const_matrix_ptr<int> nObservation, std::vector<size_t>& startIndex)
{
    size_t states = stateIndices.size();
    size_t T = observation->nrows() > 0 ? observation->nrows() : nObservation->nrows();

    std::vector<std::vector<int>> decoding;
    for (size_t k = 0; k < startIndex.size(); ++k)
    {
        size_t end = k < startIndex.size()-1 ? startIndex[k+1] : T;
        decoding.push_back(viterbi_decoding(std::make_pair(startIndex[k], end), observation, nObservation).second);
        process_state_sequence(decoding[k]);
    }

    std::vector<std::vector<size_t>> segmentLengths = calculate_segement_lengths(states, decoding);

    std::vector<double> selfP;
    std::vector<size_t> subStates;
    double p = 0.9;

    for (size_t i = 0; i < states; ++i)
    {
        double mean = std::accumulate(segmentLengths[i].begin(), segmentLengths[i].end(), 0.0) / segmentLengths[i].size();
        double r = mean * (1.0 - p) / p;
        subStates.push_back(std::min<size_t>(5, std::ceil(r)));
        selfP.push_back(r >= 1.0 ? p : lp::ext_exp(logA(stateIndices[i][0], stateIndices[i][0])));
    }

    N = std::accumulate(subStates.begin(), subStates.end(), 0);
    stateAssignment = std::vector<size_t>(N);
    size_t index = 0;
    for (size_t i = 0; i < states; ++i)
    {
        stateIndices[i] = std::vector<size_t>(subStates[i]);
        for (auto it = stateIndices[i].begin(); it != stateIndices[i].end(); ++it)
        {
            *it = index;
            stateAssignment[index] = i;
            ++index;
        }
    }

    // update transition probabilities
    init_initial();
    init_transitions(selfP);
}