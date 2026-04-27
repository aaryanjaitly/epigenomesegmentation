#include <limits>
#include <cmath>
#include <math.h>
#include <stdexcept>
#include <omp.h>
#include <thread>

#include "HMM.h"
#include "log_prob.h"

HMM::HMM(size_t states, Matrix<std::shared_ptr<DiscreteDistribution>>& emission, bool methylation):
N(states), emission(emission), methylation(methylation)
{
    if (states != emission.nrows())
        throw std::invalid_argument("Number of rows in the emmission matrix must be equal to number of states.");

    m = methylation ? emission.ncols()-1 : emission.ncols();
    init_transitions();
    init_initial();
}

void HMM::init_transitions()
{
    logA = Matrix<double> (N, N);
    double p = lp::ext_log(1.0 / N);
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            logA(i, j) = p;
        }
    }
}

void HMM::init_initial()
{
    logPi = std::vector<double> (N);
    double p = lp::ext_log(1.0 / N);
    for (size_t i = 0; i < N; ++i)
    {
        logPi[i] = p;
    }
}

void HMM::set_transitions(const Matrix<double>& A)
{
    if (A.nrows() != N || A.ncols() != N)
        throw std::invalid_argument("Transition matrix has wrong dimensions.");
    
    logA = Matrix<double> (N, N);
    for (size_t i = 0; i < A.nrows(); ++i)
    {
        double rowSum = 0.0;
        for (size_t j = 0; j < A.ncols(); ++j)
        {
            logA(i, j) = lp::ext_log(A(i, j));
            rowSum += A(i, j);
        }
    }
    normalize_transition_matrix();
}

void HMM::set_initial(const std::vector<double>& pi)
{
    if (pi.size() != N)
        throw std::invalid_argument("Initial state distribution has wrong dimensions.");
    
    logPi = std::vector<double> (N);
    for (size_t i = 0; i < pi.size(); ++i)
    {
        logPi[i] = lp::ext_log(pi[i]);
    }
}

Matrix<double> HMM::get_transition_matrix() const
{
    Matrix<double> A (N, N);

    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            A(i, j) = lp::ext_exp(logA(i, j));
        }
    }

    return A;
}

void HMM::train(const_matrix_ptr<int> observations, const_matrix_ptr<int> nObservation, std::vector<size_t>& index, double eps, size_t maxIteration, bool updateInit)
{
    if (observations->ncols() + nObservation->ncols() / 2 != emission.ncols())
        throw std::invalid_argument("Number of columns in both observation matrices must together be equal to the number of columns in the emission matrix.");

    int obs = index.size();
    size_t totalT = observations->nrows() > 0 ? observations->nrows() : nObservation->nrows();
    std::vector<matrix_ptr<double>> logAlpha (obs);
    std::vector<matrix_ptr<double>> logBeta (obs);
    std::vector<matrix_ptr<double>> logGamma (obs);
    std::vector<matrix_ptr<double>> logEmission (obs);

    for (size_t k = 0; k < obs; ++k)
    {
        size_t T = k < obs-1 ? index[k+1] - index[k] : totalT - index[k];
        logAlpha[k] = std::make_shared<Matrix<double>>(Matrix<double> (T, N));
        logBeta[k] = std::make_shared<Matrix<double>>(Matrix<double> (T, N));
        logGamma[k] = std::make_shared<Matrix<double>>(Matrix<double> (T, N));
        logEmission[k] = std::make_shared<Matrix<double>>(Matrix<double> (T, N));
    }

    size_t iteration = 0;
    double oldLikelihood = std::numeric_limits<double>::lowest();
    double newLikelihood = std::numeric_limits<double>::lowest();

    int threads = std::min(obs, omp_get_max_threads());
    do
    {
        ++iteration;
        oldLikelihood = newLikelihood;
        newLikelihood = 0.0;

        // perform forward, backward and E-step for each independent observation
        #pragma omp parallel for num_threads(threads) schedule(static)
        for (size_t k = 0; k < obs; ++k)
        {
            calculate_log_emission(logEmission[k], observations, nObservation, index[k]);
            std::thread backwardPass(&HMM::backward, this, logEmission[k], logBeta[k]); 
            double logL = 0.0;
            forward(logEmission[k], logAlpha[k], logL);
            backwardPass.join();

            #pragma omp critical
            {
                newLikelihood += logL;
            }

            E_step(logEmission[k], logAlpha[k], logBeta[k], logGamma[k]);
        }

        // perform M step for all observations
        M_step(observations, nObservation, logEmission, logGamma, logBeta, updateInit);

        if (maxIteration > 0 && iteration == maxIteration)
            break;

        if (newLikelihood < oldLikelihood)
        {
            std::cout << "Training stopped: new likelihood smaller than old likelihood." << std::endl;
            break;
        }

        if (iteration % 10 == 0)
        {
            std::cout << "Performed iteration " << iteration << " with loglikelihood " << newLikelihood << std::endl;
        }
    }
    while (fabs(newLikelihood - oldLikelihood) > eps);
    std::cout << "Training successfully completed after " << iteration << " iterations with loglikelihood " << newLikelihood << std::endl;

    int s = emission.nrows();
    // number of emission parameters
    int k = 0;
    for (size_t i = 0; i < emission.ncols(); ++i)
    {
        k += emission(0, i)->get_parameters().size();
    }
    // add extra parameters of topology model
    k = N > s ? k + 1 : k;

    // total number of free parameters
    int p = pow(s, 2) + s * k - 1;

    double AIC = -2 * newLikelihood + 2 * p;
    double BIC = -2 * newLikelihood + p * lp::ext_log(totalT);

    std::cout << "AIC: " << AIC << std::endl;
    std::cout << "BIC: " << BIC << std::endl;
}

double HMM::log_likelihood(const_matrix_ptr<int> observation, const_matrix_ptr<int> nObservation, std::pair<size_t, size_t> range) const
{    
    if (observation->ncols() + nObservation->ncols() / 2 != emission.ncols())
        throw std::invalid_argument("Number of columns in the observation matrices must together be equal to the number of columns in the emission matrix.");

    size_t totalT = observation->nrows() > 0 ? observation->nrows() : nObservation->nrows();
    size_t T = range.second > 0 ? range.second - range.first : totalT;
    matrix_ptr<double> logAlpha = std::make_shared<Matrix<double>> (Matrix<double>(T, N));
    matrix_ptr<double> logEmission = std::make_shared<Matrix<double>> (Matrix<double>(T, N));
    calculate_log_emission(logEmission, observation, nObservation, range.first);

    double logL = 0.0;
    forward(logEmission, logAlpha, logL);
    return logL;
}

std::pair<double, std::vector<int>> HMM::viterbi_decoding(std::pair<size_t, size_t> range, const_matrix_ptr<int> observation, const_matrix_ptr<int> nObservation) const
{
    if (observation->ncols() + nObservation->ncols() / 2 != emission.ncols())
        throw std::invalid_argument("Number of columns in the observation matrices must together be equal to the number of columns in the emission matrix.");
    
    size_t T = range.second - range.first;
    matrix_ptr<double> logEmission = std::make_shared<Matrix<double>> (Matrix<double>(T, N));
    Matrix<double> trellis (T, N);
    Matrix<int> backtracking (T, N);

    calculate_log_emission(logEmission, observation, nObservation, range.first);

    for (size_t i = 0; i < N; ++i)
    {
        trellis(0, i) = lp::log_mul(logPi[i], (*logEmission)(0, i));
        backtracking(0, i) = -1;
    }

    for (size_t t = 1; t < T; ++t)
    {
        for (size_t i = 0; i < N; ++i)
        {
            double maxP = lp::log_mul(trellis(t-1, 0), logA(0, i));
            int maxS = 0;

            for (size_t j = 0; j < N; ++j)
            {
                double p = lp::log_mul(trellis(t-1, j), logA(j, i));
                if (lp::log_greater(p, maxP))
                {
                    maxP = p;
                    maxS = j;
                }
            }
            trellis(t, i) = lp::log_mul((*logEmission)(t, i), maxP);
            backtracking(t, i) = maxS;
        }
    }

    // backtracking to find state sequence
    std::vector<int> stateSequence (T);
    int next = 0;
    for (size_t i = 1; i < N; ++i)
    {
        if (lp::log_greater(trellis(T-1, i), trellis(T-1, next)))
            next = i;
    }

    double maxProb = trellis(T-1, next);
    stateSequence[T-1] = next+1;
    next = backtracking(T-1, next);

    for (int t = T-1; t > 0; --t)
    {
        stateSequence[t-1] = next+1;
        next = backtracking(t-1, next);
    }

    return std::make_pair(maxProb, stateSequence);
}

std::vector<std::pair<int, double>> HMM::posterior_decoding(std::pair<size_t, size_t> range, const_matrix_ptr<int> observation, const_matrix_ptr<int> nObservation) const
{
    if (observation->ncols() + nObservation->ncols() / 2 != emission.ncols())
        throw std::invalid_argument("Number of columns in the observation matrices must together be equal to the number of columns in the emission matrix.");
    
    size_t T = range.second - range.first;
    matrix_ptr<double> logAlpha = std::make_shared<Matrix<double>> (Matrix<double>(T, N));
    matrix_ptr<double> logBeta = std::make_shared<Matrix<double>> (Matrix<double>(T, N));
    matrix_ptr<double> logEmission = std::make_shared<Matrix<double>> (Matrix<double>(T, N));
    double logLikelihood = 0.0;

    calculate_log_emission(logEmission, observation, nObservation, range.first);
    std::thread backwardPass(&HMM::backward, this, logEmission, logBeta);
    forward(logEmission, logAlpha, logLikelihood);
    backwardPass.join();    

    std::vector<std::pair<int, double>> stateSequence (T);
    for (size_t t = 0; t < T; ++t)
    {
        int state = 1;
        double prob = lp::log_mul((*logAlpha)(t, 0), (*logBeta)(t, 0));
        for (size_t i = 1; i < N; ++i)
        {
            double p = lp::log_mul((*logAlpha)(t, i), (*logBeta)(t, i));
            if (lp::log_greater(p, prob))
            {
                state = i+1;
                prob = p;
            }
        }
        double logPosterior = lp::log_mul(prob, -logLikelihood);
        stateSequence[t] = std::make_pair(state, lp::ext_exp(logPosterior));
    }
    return stateSequence;
}

std::vector<HMM::matrix_ptr<double>> HMM::membership_coefficients(std::vector<size_t>& index, const_matrix_ptr<int> observation, const_matrix_ptr<int> nObservation) const
{
    size_t obs = index.size();
    size_t totalT = observation->nrows() > 0 ? observation->nrows() : nObservation->nrows();
    std::vector<matrix_ptr<double>> logAlpha (obs);
    std::vector<matrix_ptr<double>> logBeta (obs);
    std::vector<matrix_ptr<double>> logGamma (obs);
    std::vector<matrix_ptr<double>> logEmission (obs);

    for (size_t k = 0; k < obs; ++k)
    {
        size_t T = k < obs-1 ? index[k+1] - index[k] : totalT - index[k];
        logAlpha[k] = std::make_shared<Matrix<double>>(Matrix<double> (T, N));
        logBeta[k] = std::make_shared<Matrix<double>>(Matrix<double> (T, N));
        logGamma[k] = std::make_shared<Matrix<double>>(Matrix<double> (T, N));
        logEmission[k] = std::make_shared<Matrix<double>>(Matrix<double> (T, N));
    }

    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < obs; ++k)
    {
        calculate_log_emission(logEmission[k], observation, nObservation, index[k]);
        std::thread backwardPass(&HMM::backward, this, logEmission[k], logBeta[k]); 
        double logL = 0.0;
        forward(logEmission[k], logAlpha[k], logL);
        backwardPass.join();
        E_step(logEmission[k], logAlpha[k], logBeta[k], logGamma[k]);
    }
    return logGamma;
}

void HMM::calculate_log_emission(matrix_ptr<double> logEmission, const_matrix_ptr<int> observation, const_matrix_ptr<int> nObservation, size_t start) const
{    
    size_t T = logEmission->nrows();
    #pragma omp parallel for collapse(2) schedule(static)
    for (size_t t = 0; t < T; ++t)
    {
        for (size_t i = 0; i < N; ++i)
        {
            double p = 0.0;
            if (m > 0)
            {
                auto it = observation->row_begin(start+t);
                //marks are assumed to be independent
                for (size_t k = 0; k < m; ++k)
                {
                    p = lp::log_mul(p, emission(i, k)->log_pmf(*it));
                    ++it;
                }
            }
            if (methylation)
            {
                auto mIt = nObservation->row_begin(start+t);
                std::shared_ptr<TwoValueDiscreteDistribution> dis = std::dynamic_pointer_cast<TwoValueDiscreteDistribution>(emission(i, m));
                p = lp::log_mul(p, dis->log_pmf((*nObservation)(start+t, 0), (*nObservation)(start+t, 1)));
            }
            (*logEmission)(t, i) = p;
        }
    }
}

void HMM::forward(const_matrix_ptr<double> logEmission, matrix_ptr<double> logAlpha, double& likelihood) const
{
    size_t T = logEmission->nrows();
    
    for (size_t i = 0; i < N; ++i)
    {
        (*logAlpha)(0, i) = lp::log_mul(logPi[i], (*logEmission)(0, i));
    }

    for (size_t t = 1; t < T; ++t)
    {
        for (size_t i = 0; i < N; ++i)
        {
            double alpha = NAN;
            for (size_t j = 0; j < N; ++j)
            {
                alpha = lp::log_add(alpha, lp::log_mul((*logAlpha)(t-1, j), logA(j, i)));
            }
            (*logAlpha)(t, i) = lp::log_mul(alpha, (*logEmission)(t, i));
        }
    }
    likelihood = lp::log_sum(logAlpha->row_begin(T-1), logAlpha->row_end(T-1));
}

void HMM::backward(const_matrix_ptr<double> logEmission, matrix_ptr<double> logBeta) const
{
    size_t T = logEmission->nrows();
    for (size_t i = 0; i < N; ++i)
    {
        (*logBeta)(T-1, i) = 0;
    }

    for (size_t t = 2; t <= T; ++t)
    {
        for (size_t i = 0; i < N; ++i)
        {
            double beta = NAN;
            for (size_t j = 0; j < N; ++j)
            {
                beta = lp::log_add(beta, lp::log_mul(logA(i, j), lp::log_mul((*logEmission)(T-t+1, j), (*logBeta)(T-t+1, j))));
            }
            (*logBeta)(T-t, i) = beta;
        }
    }
}

void HMM::E_step(const_matrix_ptr<double> logEmission, matrix_ptr<double> logAlpha, matrix_ptr<double> logBeta, matrix_ptr<double> logGamma) const
{
    size_t T = logEmission->nrows();

    #pragma omp parallel for schedule(static)
    for (size_t t = 0; t < T; ++t)
    {
        double denom = NAN;
        for (size_t i = 0; i < N; ++i)
        {
            (*logGamma)(t, i) = lp::log_mul((*logAlpha)(t, i), (*logBeta)(t, i));
            denom = lp::log_add(denom, (*logGamma)(t, i));
        }

        for (size_t i = 0; i < N; ++i)
        {
            (*logGamma)(t, i) = lp::log_mul((*logGamma)(t, i), -denom);
        }
    }
}

void HMM::M_step(const_matrix_ptr<int> observations, const_matrix_ptr<int> nObservations, std::vector<matrix_ptr<double>>& logEmission, std::vector<matrix_ptr<double>>& logGamma, std::vector<matrix_ptr<double>>& logBeta, bool update_init)
{    
    size_t obs = logGamma.size();
    size_t totalT = observations->nrows() > 0 ?  observations->nrows() : nObservations->nrows();
    double lowerProbBound = lp::ext_log(pow(10.0, -10.0));
    bool normalize = false;
    Matrix<double> new_logA = Matrix<double> (N, N);
    
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < N; ++i)
    {
        //update initial state distribution if required
        if (update_init)
        {
            logPi[i] = (*logGamma[0])(0, i);
            for (size_t k = 1; k < obs; ++k)
            {
                logPi[i] = lp::log_add(logPi[i], (*logGamma[k])(0, i));
            }
            logPi[i] -= lp::ext_log(obs);

            if (lp::log_greater(lowerProbBound, logPi[i]))
            {
                logPi[i] = NAN;
                #pragma omp critical
                {
                    normalize = true;
                }
            }
        }

        double denom = NAN;
        std::vector<double> completeGamma (totalT);
        size_t curr = 0;
        for (size_t k = 0; k < obs; ++k)
        {
            size_t T = logGamma[k]->nrows();
            for (size_t t = 0; t < T-1; ++t)
            {
                double gamma = (*logGamma[k])(t, i);
                denom = lp::log_add(denom, gamma);
                completeGamma[curr] = lp::ext_exp(gamma);
                ++curr;
            }
            completeGamma[curr] = lp::ext_exp((*logGamma[k])(T-1, i));
            ++curr;
        }

        //update transition matrix
        for (size_t j = 0; j < N; ++j)
        {
            double num = NAN;
            for (size_t k = 0; k < obs; ++k)
            {
                size_t T = logGamma[k]->nrows();
                for (size_t t = 0; t < T-1; ++t)
                {
                    double xi_tij = (*logGamma[k])(t, i) + logA(i, j) + (*logEmission[k])(t+1, j) + (*logBeta[k])(t+1, j) - (*logBeta[k])(t, i);
                    num = lp::log_add(num, xi_tij);
                }
            }
            new_logA(i, j) = lp::log_mul(num, -denom);
            //if probability is < 10^-10 set it to zero
            if (lp::log_greater(lowerProbBound, new_logA(i, j)))
            {
                new_logA(i, j) = NAN;
                #pragma omp critical
                {
                    normalize = true;
                }
            }
        }

        //only update parameters if membership coefficients are not all zero
        if (std::accumulate(completeGamma.begin(), completeGamma.end(), 0.0) == 0.0)
        {
            #pragma omp critical
            std::cerr << "Membership coefficient are zero for state " << i+1 << std::endl;
        }
        else
        {
            for (size_t k = 0; k < m; ++k)
            {
                emission(i, k)->update(completeGamma.begin(), completeGamma.end(), observations->col_begin(k), observations->col_end(k));
            }
            if (methylation)
            {
                std::shared_ptr<TwoValueDiscreteDistribution> dis = std::dynamic_pointer_cast<TwoValueDiscreteDistribution>(emission(i, m));
                dis->update_methylation(completeGamma.begin(), completeGamma.end(), nObservations->col_begin(0), nObservations->col_begin(1));
            }
        }
    }

    logA = new_logA;

    //if a value of the transition matrix was set to zero, normalize such that all rows sum up to one again
    if (normalize)
    {
        normalize_transition_matrix();
    }
}

void HMM::normalize_transition_matrix()
{
    double norm = lp::log_sum(logPi.begin(), logPi.end());
    for (auto& p : logPi)
    {
        p = lp::log_mul(p, -norm);
    }

    for (size_t i = 0; i < N; ++i)
    {
        norm = lp::log_sum(logA.row_begin(i), logA.row_end(i));

        for (size_t j = 0; j < N; ++j)
        {
            logA(i, j) = lp::log_mul(logA(i, j), -norm);
        }
    }
}