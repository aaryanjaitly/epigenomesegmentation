#ifndef _LOGL_FCN_H
#define _LOGL_FCN_H

#include <vector>
#include "distribution.h"

/**
 * @brief Functor calculating log likelihood to pass to minimizer
 * 
 * @param T function to calculate log pmf (double*, std::vector<double>&) -> double
 * 
**/
template<typename T>
class LogLikelihoodFunction
{
    public:
    LogLikelihoodFunction(DiscreteDistribution::g_iterator gBegin, DiscreteDistribution::g_iterator gEnd, 
                            DiscreteDistribution::o_iterator oBegin, DiscreteDistribution::o_iterator oEnd, std::vector<double>& lp, T f):
                            gammaBegin(gBegin), gammaEnd(gEnd), obsBegin(oBegin), obsEnd(oEnd), calculate_pmf(f) {logpmf = lp;}
    
    double operator () (const double* x) 
    {
        // try-catch block needed for Sichel distribution (due to inprecison and overflow problems of Bessel function)
        try
        {
            calculate_pmf(x, logpmf);
            if (isnan(logpmf[0]) || isnan(logpmf[1]))
            {
                return std::numeric_limits<double>::max();
            }
        }
        catch(const std::exception& e)
        {
            return std::numeric_limits<double>::max();
        }

        double logL = 0.0;
        auto obsIt = obsBegin;
        for (auto it = gammaBegin; it != gammaEnd; ++it)
        {
            if (*it > 0)
            {
                logL += *it * logpmf[*obsIt];
            }
            ++obsIt;
        }
        return - logL;
    }

    private:
    DiscreteDistribution::g_iterator gammaBegin, gammaEnd;
    DiscreteDistribution::o_iterator obsBegin, obsEnd;
    std::vector<double> logpmf;
    T calculate_pmf;
};

/**
 * @brief Functor calculating log likelihood to pass to minimizer for beta binomial (fixed parameter n)
 * 
 * @param T function to calculate log pmf (double*, int, std::vector<double>&) -> double
 * 
**/
template<typename T>
class LogLikelihoodFunctionBB
{
    public:
    LogLikelihoodFunctionBB(DiscreteDistribution::g_iterator gBegin, DiscreteDistribution::g_iterator gEnd, 
                            DiscreteDistribution::o_iterator oBegin, DiscreteDistribution::o_iterator oEnd, std::vector<double>& lp, int n, T f):
                            gammaBegin(gBegin), gammaEnd(gEnd), obsBegin(oBegin), obsEnd(oEnd), calculate_pmf(f) {logpmf = lp; N = n;}
    
    double operator () (const double* x) 
    {
        calculate_pmf(x, N, logpmf);

        double logL = 0.0;
        auto obsIt = obsBegin;
        for (auto it = gammaBegin; it != gammaEnd; ++it)
        {
            if (*it > 0)
            {
                logL += *it * logpmf[*obsIt];
            }
            ++obsIt;
        }
        return - logL;
    }

    private:
    DiscreteDistribution::g_iterator gammaBegin, gammaEnd;
    DiscreteDistribution::o_iterator obsBegin, obsEnd;
    std::vector<double> logpmf;
    T calculate_pmf;
    int N;

};

/**
 * @brief Functor calculating log likelihood to pass to minimizer for DNA methylation data
 * 
 * @param T function to calculate log pmf (double*, int, int) -> double
 * 
**/
template<typename T>
class LogLikelihoodFunctionMeth
{
    public:
    LogLikelihoodFunctionMeth(TwoValueDiscreteDistribution::dna_g_iterator gBegin, TwoValueDiscreteDistribution::dna_g_iterator gEnd, 
                            TwoValueDiscreteDistribution::dna_o_iterator cBegin, TwoValueDiscreteDistribution::dna_o_iterator mBegin, T f):
                            gammaBegin(gBegin), gammaEnd(gEnd), covBegin(cBegin), methBegin(mBegin), calculate_pmf(f) {}
    
    double operator () (const double* x) 
    {
        double logL = 0.0;
        auto covIt = covBegin;
        auto methIt = methBegin;
        for (auto it = gammaBegin; it != gammaEnd; ++it)
        {
            if (*it > 0)
            {
                logL += (*it) * calculate_pmf(x, *covIt, *methIt);
            }
            ++covIt;
            ++methIt;
        }
        return - logL;
    }

    private:
    TwoValueDiscreteDistribution::dna_g_iterator gammaBegin, gammaEnd;
    TwoValueDiscreteDistribution::dna_o_iterator covBegin, methBegin;
    T calculate_pmf;
};


#endif