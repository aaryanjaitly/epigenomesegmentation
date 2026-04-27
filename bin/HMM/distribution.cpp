#include <cmath>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <functional>

#include "boost/math/special_functions/digamma.hpp"
#include "boost/math/special_functions/bessel.hpp"
#include "boost/math/special_functions/beta.hpp" 

#include "FCN.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "distribution.h"
#include "log_prob.h"
#include "sample_statistics.h"

namespace bmath = boost::math;

Poisson::Poisson(double l): lambda(l)
{
    if (lambda <= 0.0)
        throw std::invalid_argument("Lambda parameter of a Poisson distribution must be > 0.");
}

double Poisson::pmf(int x) const
{
    return lp::ext_exp(log_pmf(x));
}

double Poisson::log_pmf(int x) const
{
    return x * log(lambda) - lambda - lgamma(x+1);
}

std::vector<double> Poisson::get_parameters() const
{
    return std::vector<double>{lambda};
}

void Poisson::update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{    
    double sampleMean = stat::sample_mean(gammaBegin, gammaEnd, obsBegin, obsEnd);
    if (sampleMean > 0.0)
    {
        lambda = sampleMean;
    }
}

ZeroAdjustedPoisson::ZeroAdjustedPoisson(double l, double pi): lambda(l), pi(pi)
{
    if (lambda <= 0.0)
        throw std::invalid_argument("Lambda parameter of a Poisson distribution must be > 0.");

    if (pi < 0.0 || pi > 1.0)
        throw std::invalid_argument("Zero mixing parameter must be a probability in the interval [0, 1].");

    minimizer = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    minimizer->SetMaxFunctionCalls(300);
    minimizer->SetTolerance(1e-16);
}

double ZeroAdjustedPoisson::pmf(int x) const
{
    return lp::ext_exp(log_pmf(x));
}

double ZeroAdjustedPoisson::log_pmf(int x) const
{
    if (x == 0)
        return lp::ext_log(pi);
    
    if (pi > 1.0 - std::numeric_limits<double>::epsilon())
        return NAN;

    double prob0 = lp::ext_exp(-lambda);
    return lp::ext_log(1.0-pi) + x * log(lambda) - lambda - lgamma(x+1) - lp::log1x(-prob0);
}

std::vector<double> ZeroAdjustedPoisson::get_parameters() const
{
    return std::vector<double>{lambda, pi};
}

void ZeroAdjustedPoisson::update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{   
    pi = stat::zero_frequency(gammaBegin, gammaEnd, obsBegin, obsEnd);
    
    if (pi > 1.0 - std::numeric_limits<double>::epsilon())
    {
        pi = 1.0;
        lambda = 1.0;
    }
    else
    {
        using namespace ROOT::Math;

        int max = *std::max_element(obsBegin, obsEnd);
        std::vector<double> lpmf (max+1);
        LogLikelihoodFunction L (gammaBegin, gammaEnd, obsBegin, obsEnd, lpmf, precalculate_log_pmf);
        Functor f(L, 1);

        minimizer->SetFunction(f);
        minimizer->SetLowerLimitedVariable(1, "lambda", lambda, 0.01, std::numeric_limits<double>::min());
        minimizer->Minimize();

        double params [] = {lambda};
        double logL = L(&params[0]);

        if (logL > minimizer->MinValue() && minimizer->X()[0] >= 0.0)
        {
            lambda = minimizer->X()[0];
        }
    }
}

void ZeroAdjustedPoisson::precalculate_log_pmf(const double* params, std::vector<double>& logpmf)
{
    logpmf[0] = 0.0;

    double lambda = params[0];
    double prob0 = lp::ext_exp(-lambda);
  
    for (size_t i = 1; i < logpmf.size(); ++i)
    {
        logpmf[i] = i * log(lambda) - lambda - lgamma(i+1) - lp::log1x(-prob0);
    }
}

Binomial::Binomial(double p, int n): p(p), n(n)
{
    if (n <= 0)
        throw std::invalid_argument("N must be greater than 0.");
    
    if (p < 0.0 || p > 1.0)
        throw std::invalid_argument("P must be in the interval [0, 1].");

    lgammaN = lgamma(n+1);
}

Binomial::Binomial(double p): p(p)
{    
    if (p < 0.0 || p > 1.0)
        throw std::invalid_argument("P must be in the interval [0, 1].");
    n = 0;
}

double Binomial::pmf(int x) const
{
    return lp::ext_exp(log_pmf(x));
}

double Binomial::log_pmf(int x) const
{
    if (x > n)
        return NAN;
    return lgammaN - lgamma(x+1) - lgamma(n-x+1) + lp::xlogy(x, p) + lp::xlogy(n-x, 1-p); 
}

double Binomial::pmf(int n, int x) const
{
    return lp::ext_exp(log_pmf(n, x));
}

double Binomial::log_pmf(int n, int x) const
{
    if (n < 0)
        throw std::invalid_argument("Invalid parameter in Binomial distribution (n < 0)");

    if (x > n)
        return NAN;

    if (n == 0 & x == 0)
        return 0.0;
        
    return lgamma(n+1) - lgamma(x+1) - lgamma(n-x+1) + lp::xlogy(x, p) + lp::xlogy(n-x, 1-p); 
}

std::vector<double> Binomial::get_parameters() const
{
    return std::vector<double>{p, (double)n};
}

void Binomial::update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{
    double sampleMean = stat::sample_mean(gammaBegin, gammaEnd, obsBegin, obsEnd); 
    p = sampleMean / n;
}

void Binomial::update_methylation(dna_g_iterator gammaBegin, dna_g_iterator gammaEnd, dna_o_iterator covBegin, dna_o_iterator methBegin)
{
    auto covIt = covBegin;
    auto methIt = methBegin;
    double X = 0;
    double N = 0;
    for (auto it = gammaBegin; it != gammaEnd; ++it)
    {
        N += (*it) * (*covIt);
        X += (*it) * (*methIt);
        ++covIt;
        ++methIt;
    }
    p = X / N;
}

NegativeBinomial::NegativeBinomial(double p, double r): p(p), r(r)
{
    if (r <= 0.0)
        throw std::invalid_argument("Number of successes for the negative binomial distribution must be > 0.");

    if (p < 0.0 || p > 1.0)
        throw std::invalid_argument("Success probability for the negative binomial distribution must be in the interval [0, 1].");

    lgammaR = lgamma(r);

    minimizer = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    minimizer->SetMaxFunctionCalls(300);
    minimizer->SetTolerance(1e-16);
}

double NegativeBinomial::pmf(int x) const
{
    return lp::ext_exp(log_pmf(x));
}

double NegativeBinomial::log_pmf(int x) const
{
    return lp::xlogy(r, p) + lp::xlogy(x, 1.0-p) + lgamma(x+r) - lgamma(x+1) - lgammaR;
}

std::vector<double> NegativeBinomial::get_parameters() const
{
    return std::vector<double>{p, r};
}

void NegativeBinomial::update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{    
    double sampleMean = stat::sample_mean(gammaBegin, gammaEnd, obsBegin, obsEnd);
    double sampleVariance = stat::sample_variance(sampleMean, gammaBegin, gammaEnd, obsBegin, obsEnd);

    if (sampleVariance <= 0.0 || sampleVariance - sampleMean <= 0.0)
    {
        return;
    }

    double newR = std::pow(sampleMean, 2) / (sampleVariance - sampleMean);
    double newP = sampleMean / sampleVariance;
    if (!std::isfinite(newR))
    {
        return;
    }
    find_mle(newP, newR, gammaBegin, gammaEnd, obsBegin, obsEnd);
}

void NegativeBinomial::find_mle(double initP, double initR, g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{
    using namespace ROOT::Math;

    int max = *std::max_element(obsBegin, obsEnd);
    std::vector<double> lpmf (max+1);
    LogLikelihoodFunction L (gammaBegin, gammaEnd, obsBegin, obsEnd, lpmf, precalculate_log_pmf);
    Functor f(L, 2);

    minimizer->SetFunction(f);
    minimizer->SetLimitedVariable(0, "p", initP, 0.01, 0.0, 1.0);
    minimizer->SetLowerLimitedVariable(1, "r", initR, 0.01, std::numeric_limits<double>::min());
    minimizer->Minimize();

    double params [] = {p, r};
    double logL = L(&params[0]);

    if (logL > minimizer->MinValue() && minimizer->X()[0] >= 0.0 && minimizer->X()[0] <= 1.0 && minimizer->X()[1] > 0.0)
    {
        p = minimizer->X()[0];
        r = minimizer->X()[1];
        lgammaR = lgamma(r);
    }
}

void NegativeBinomial::precalculate_log_pmf(const double* params, std::vector<double>& logpmf)
{
    double p = params[0];
    double r = params[1];
    for (size_t i = 0; i < logpmf.size(); ++i)
    {
        logpmf[i] = lp::xlogy(r, p) + lp::xlogy(i, 1.0-p) + lgamma(i+r) - lgamma(i+1) - lgamma(r);
    }
}

ZeroAdjustedNegativeBinomial::ZeroAdjustedNegativeBinomial(double p, double r, double pi): p(p), r(r), pi(pi)
{
    if (r <= 0.0)
        throw std::invalid_argument("Number of successes for the negative binomial distribution must be > 0.");

    if (p < 0.0 || p > 1.0)
        throw std::invalid_argument("Success probability for the negative binomial distribution must be in the interval [0, 1].");

    if (pi < 0.0 || pi > 1.0)
        throw std::invalid_argument("Zero probability must be in the interval [0, 1].");

    lgammaR = lgamma(r);

    minimizer = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    minimizer->SetMaxFunctionCalls(300);
    minimizer->SetTolerance(1e-16);
}

double ZeroAdjustedNegativeBinomial::pmf(int x) const
{
    return lp::ext_exp(log_pmf(x));
}

double ZeroAdjustedNegativeBinomial::log_pmf(int x) const
{
    if (x == 0)
        return lp::ext_log(pi);

    if (pi > 1.0 - std::numeric_limits<double>::epsilon())
        return NAN;

    double link = lp::log1x(-lp::ext_exp(lp::xlogy(r, p) + lp::xlogy(0, 1.0-p) + lgamma(0+r) - lgamma(0+1) - lgammaR));
    return lp::ext_log(1.0-pi) + lp::xlogy(r, p) + lp::xlogy(x, 1.0-p) + lgamma(x+r) - lgamma(x+1) - lgammaR - link;
}

std::vector<double> ZeroAdjustedNegativeBinomial::get_parameters() const
{
    return std::vector<double>{p, r, pi};
}

void ZeroAdjustedNegativeBinomial::update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{    
    pi = stat::zero_frequency(gammaBegin, gammaEnd, obsBegin, obsEnd);
    if (pi > 1.0 - std::numeric_limits<double>::epsilon())
    {
        pi = 1.0;
        p = 1.0;
    }
    else
    {
        using namespace ROOT::Math;

        int max = *std::max_element(obsBegin, obsEnd);
        std::vector<double> lpmf (max+1);
        LogLikelihoodFunction L (gammaBegin, gammaEnd, obsBegin, obsEnd, lpmf, precalculate_log_pmf);
        Functor f(L, 2);

        minimizer->SetFunction(f);
        minimizer->SetLimitedVariable(0, "p", p, 0.01, 0.0, 1.0);
        minimizer->SetLowerLimitedVariable(1, "r", r, 0.01, std::numeric_limits<double>::min());
        minimizer->Minimize();

        double params [] = {p, r};
        double logL = L(&params[0]);

        if (logL > minimizer->MinValue() && minimizer->X()[0] >= 0.0 && minimizer->X()[0] <= 1.0 && minimizer->X()[1] > 0.0)
        {
            p = minimizer->X()[0];
            r = minimizer->X()[1];
            lgammaR = lgamma(r);
        }
    }
}

void ZeroAdjustedNegativeBinomial::precalculate_log_pmf(const double* params, std::vector<double>& logpmf)
{
    NegativeBinomial::precalculate_log_pmf(params, logpmf);

    double link = lp::log1x(-lp::ext_exp(logpmf[0]));
    std::for_each(logpmf.begin(), logpmf.end(), [link](double& p){p = p - link;});
    logpmf[0] = 0.0;
}

BetaBinomial::BetaBinomial(double a, double b, int n): alpha(a), beta(b), n(n)
{
    if (n <= 0)
        throw std::invalid_argument("N must be greater than 1.");
    
    if (alpha <= 0.0 || beta <= 0.0)
        throw std::invalid_argument("Both alpha and beta of the beta binomial distribution must be > 0.");

    lgammaN1 = lgamma(n+1);
    lgammaNAB = lgamma(n+alpha+beta);
    lgammaAB = lgamma(alpha+beta);
    lgammaA = lgamma(alpha);
    lgammaB = lgamma(beta);

    minimizer = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    minimizer->SetMaxFunctionCalls(300);
    minimizer->SetTolerance(1e-16);
}

BetaBinomial::BetaBinomial(double a, double b): alpha(a), beta(b)
{   
    if (alpha <= 0.0 || beta <= 0.0)
        throw std::invalid_argument("Both alpha and beta of the beta binomial distribution must be > 0.");

    lgammaAB = lgamma(alpha+beta);
    lgammaA = lgamma(alpha);
    lgammaB = lgamma(beta);
    lgammaN1 = 0.0;
    lgammaNAB = 0.0;
    n = 0;

    minimizer = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    minimizer->SetMaxFunctionCalls(300);
    minimizer->SetTolerance(1e-16);
}

double BetaBinomial::pmf(int x) const
{
    return lp::ext_exp(log_pmf(x));
}

double BetaBinomial::log_pmf(int x) const
{
    if (n == 0)
        throw std::invalid_argument("Invalid parameter in Beta Binomial distribution (n = 0)");

    if (x > n)
        return NAN;

    return lgammaN1 + lgamma(x+alpha) + lgamma(n-x+beta) + lgammaAB - (lgamma(x+1) + lgamma(n-x+1) + lgammaNAB + lgammaA + lgammaB);
}

double BetaBinomial::pmf(int n, int x) const
{
    return lp::ext_exp(log_pmf(n, x));
}

double BetaBinomial::log_pmf(int n, int x) const
{
    if (n < 0)
        throw std::invalid_argument("Number of trials must be > 0 for the beta binomial distribution.");

    if (x > n)
        return NAN;

    if (x == 0 && n == 0)
        return 0.0;

    return lgamma(n+1) + lgamma(x+alpha) + lgamma(n-x+beta) + lgammaAB - (lgamma(x+1) + lgamma(n-x+1) + lgamma(n+alpha+beta) + lgammaA + lgammaB);
}

std::vector<double> BetaBinomial::get_parameters() const
{
    return std::vector<double>{alpha, beta, (double)n};
}

void BetaBinomial::update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{
    double mean = stat::sample_mean(gammaBegin, gammaEnd, obsBegin, obsEnd);
    double var = stat::sample_variance(mean, gammaBegin, gammaEnd, obsBegin, obsEnd);

    //method of moment estimators
    double pi = mean / n;
    if (pi > 1.0)
    {
        std::clog << "The mean value of the binomial parameter p should be in [0, 1] (wrong parameter n chosen?) \n No update perfomed!\n";
        return;
    }
    double theta = (var - n * pi * (1.0 - pi)) / (std::pow(n, 2) * pi * (1.0 - pi) - var);

    if (pi > 0.0 && theta > 0.0)
    {
        find_mle(gammaBegin, gammaEnd, obsBegin, obsEnd, pi / theta, (1 - pi) / theta);

        //update helper variables
        lgammaNAB = lgamma(n+alpha+beta);
        lgammaAB = lgamma(alpha+beta);
        lgammaA = lgamma(alpha);
        lgammaB = lgamma(beta);
    }
}

void BetaBinomial::update_methylation(dna_g_iterator gammaBegin, dna_g_iterator gammaEnd, dna_o_iterator covBegin, dna_o_iterator methBegin)
{
    using namespace ROOT::Math;

    LogLikelihoodFunctionMeth L (gammaBegin, gammaEnd, covBegin, methBegin, calculate_log_pmf);
    Functor f(L, 2);

    minimizer->SetFunction(f);
    minimizer->SetLowerLimitedVariable(0, "alpha", alpha, 0.01, std::numeric_limits<double>::min());
    minimizer->SetLowerLimitedVariable(1, "beta", beta, 0.01, std::numeric_limits<double>::min());
    minimizer->Minimize();

    double p [] = {alpha, beta};
    double logL = L(&p[0]);

    if (logL > minimizer->MinValue() && minimizer->X()[0] > 0 && minimizer->X()[1] > 0)
    {
        alpha = minimizer->X()[0];
        beta = minimizer->X()[1];
    }

    //update helper variables
    lgammaAB = lgamma(alpha+beta);
    lgammaA = lgamma(alpha);
    lgammaB = lgamma(beta);
}

double BetaBinomial::calculate_log_pmf(const double* params, int n, int x)
{
    if (x > n)
        return NAN;

    if (x == 0 && n == 0)
        return 0.0;

    double alpha = params[0];
    double beta = params[1];
    return (lgamma(n+1) + lgamma(x+alpha) + lgamma(n-x+beta) + lgamma(alpha+beta) - (lgamma(x+1) + lgamma(n-x+1) + lgamma(n+alpha+beta) + lgamma(alpha) + lgamma(beta)));
}

void BetaBinomial::find_mle(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd, double initA, double initB)
{
    using namespace ROOT::Math;

    int max = *std::max_element(obsBegin, obsEnd);
    std::vector<double> lpmf (max+1);
    LogLikelihoodFunctionBB L (gammaBegin, gammaEnd, obsBegin, obsEnd, lpmf, n, precalculate_log_pmf);
    Functor f(L, 2);

    minimizer->SetFunction(f);
    minimizer->SetLowerLimitedVariable(0, "alpha", initA, 0.01, std::numeric_limits<double>::min());
    minimizer->SetLowerLimitedVariable(1, "beta", initB, 0.01, std::numeric_limits<double>::min());
    minimizer->Minimize();

    double p [] = {alpha, beta};
    double logL = L(&p[0]);

    if (logL > minimizer->MinValue() && minimizer->X()[0] > 0 && minimizer->X()[1] > 0)
    {
        alpha = minimizer->X()[0];
        beta = minimizer->X()[1];
    }
}

void BetaBinomial::precalculate_log_pmf(const double* params, int n, std::vector<double>& logpmf)
{
    double alpha = params[0];
    double beta = params[1];

    for (size_t i = 0; i < logpmf.size(); ++i)
    {
        logpmf[i] = lgamma(n+1) + lgamma((i) + alpha) + lgamma(n-i+beta) + lgamma(alpha+beta) 
                    - (lgamma(i+1) + lgamma(n-i+1) + lgamma(n + alpha + beta) + lgamma(alpha) + lgamma(beta));
    }
}

BetaNegativeBinomial::BetaNegativeBinomial(double a, double b, double r): alpha(a), beta(b), r(r)
{
    if (alpha <= 1.0 || beta <= 0.0 || r < 1.0)
        throw std::invalid_argument("Alpha and r must be > 1 and beta > 0.");

    lgammaR = lgamma(r);
    lgammaAR = lgamma(alpha + r);
    lgammaAB = lgamma(alpha + beta);
    lgammaA = lgamma(alpha);
    lgammaB = lgamma(beta); 

    minimizer = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    minimizer->SetMaxFunctionCalls(300);
    minimizer->SetTolerance(1e-16);
}

double BetaNegativeBinomial::pmf(int x) const
{
    return lp::ext_exp(log_pmf(x));
}

double BetaNegativeBinomial::log_pmf(int x) const
{
    return lgamma(r + x) + lgammaAR + lgamma(beta + x) + lgammaAB - (lgamma(x+1) + lgammaR + lgamma(alpha + r + beta + x) + lgammaA + lgammaB); 
}

std::vector<double> BetaNegativeBinomial::get_parameters() const
{
    return std::vector<double>{alpha, beta, r};
}

void BetaNegativeBinomial::update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{
    find_mle(gammaBegin, gammaEnd, obsBegin, obsEnd);

    lgammaR = lgamma(r);
    lgammaAR = lgamma(alpha + r);
    lgammaAB = lgamma(alpha + beta);
    lgammaA = lgamma(alpha);
    lgammaB = lgamma(beta); 
}

void BetaNegativeBinomial::find_mle(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{
    using namespace ROOT::Math;

    int max = *std::max_element(obsBegin, obsEnd);
    std::vector<double> lpmf (max+1);
    LogLikelihoodFunction L (gammaBegin, gammaEnd, obsBegin, obsEnd, lpmf, precalculate_log_pmf);
    Functor f(L, 3);

    minimizer->SetFunction(f);
    minimizer->SetLowerLimitedVariable(0, "alpha", alpha, 0.01, 1.0);
    minimizer->SetLowerLimitedVariable(1, "beta", beta, 0.01, std::numeric_limits<double>::min());
    minimizer->SetLowerLimitedVariable(2, "r", r, 0.01, 1.0);
    minimizer->Minimize();

    double p [] = {alpha, beta, r};
    double logL = L(&p[0]);

    if (logL > minimizer->MinValue() && minimizer->X()[0] > 1 && minimizer->X()[1] > 0 && minimizer->X()[2] >= 1)
    {
        alpha = minimizer->X()[0];
        beta = minimizer->X()[1];
        r = minimizer->X()[2];
    }
}

void BetaNegativeBinomial::precalculate_log_pmf(const double* params, std::vector<double>& pmf)
{
    double alpha = params[0];
    double beta = params[1];
    double r = params[2];
  
    for (size_t i = 0; i < pmf.size(); ++i)
    {
        pmf[i] = lgamma(r + (i)) + lgamma(alpha + r) + lgamma(beta + i) + lgamma(alpha+beta) - (lgamma(i+1) + lgamma(r) + lgamma(alpha + r + beta + i) + lgamma(alpha) + lgamma(beta));
    }
}


ZeroAdjustedBetaNegativeBinomial::ZeroAdjustedBetaNegativeBinomial(double a, double b, double r, double pi): alpha(a), beta(b), r(r), pi(pi)
{
    if (alpha <= 1.0 || beta <= 0.0 || r < 1.0)
        throw std::invalid_argument("Alpha and r must be > 1 and beta > 0.");

    if (pi < 0.0 || pi > 1.0)
        throw std::invalid_argument("Zero probability must be in the interval [0, 1].");

    lgammaR = lgamma(r);
    lgammaAR = lgamma(alpha + r);
    lgammaAB = lgamma(alpha + beta);
    lgammaA = lgamma(alpha);
    lgammaB = lgamma(beta); 

    minimizer = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    minimizer->SetMaxFunctionCalls(300);
    minimizer->SetTolerance(1e-16);
}

double ZeroAdjustedBetaNegativeBinomial::pmf(int x) const
{
    return lp::ext_exp(log_pmf(x));
}

double ZeroAdjustedBetaNegativeBinomial::log_pmf(int x) const
{
    if (x == 0)
        return lp::ext_log(pi);

    if (pi > 1.0 - std::numeric_limits<double>::epsilon())
        return NAN;

    double link = lp::log1x(-lp::ext_exp(lgamma(r + 0) + lgammaAR + lgamma(beta + 0) + lgammaAB - (lgamma(0+1) + lgammaR + lgamma(alpha + r + beta + 0) + lgammaA + lgammaB)));
    return lp::ext_log(1.0-pi) + lgamma(r + x) + lgammaAR + lgamma(beta + x) + lgammaAB - (lgamma(x+1) + lgammaR + lgamma(alpha + r + beta + x) + lgammaA + lgammaB) - link; 
}

std::vector<double> ZeroAdjustedBetaNegativeBinomial::get_parameters() const
{
    return std::vector<double>{alpha, beta, r, pi};
}

void ZeroAdjustedBetaNegativeBinomial::update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{
    find_mle(gammaBegin, gammaEnd, obsBegin, obsEnd);

    lgammaR = lgamma(r);
    lgammaAR = lgamma(alpha + r);
    lgammaAB = lgamma(alpha + beta);
    lgammaA = lgamma(alpha);
    lgammaB = lgamma(beta); 
}

void ZeroAdjustedBetaNegativeBinomial::find_mle(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{
    pi = stat::zero_frequency(gammaBegin, gammaEnd, obsBegin, obsEnd);
    if (pi > 1.0 - std::numeric_limits<double>::epsilon())
    {
        pi = 1.0;
    }
    else
    {
        using namespace ROOT::Math;

        int max = *std::max_element(obsBegin, obsEnd);
        std::vector<double> lpmf (max+1);

        LogLikelihoodFunction L (gammaBegin, gammaEnd, obsBegin, obsEnd, lpmf, precalculate_log_pmf);
        Functor f(L, 3);

        minimizer->SetFunction(f);
        minimizer->SetLowerLimitedVariable(0, "alpha", alpha, 0.01, 1.0);
        minimizer->SetLowerLimitedVariable(1, "beta", beta, 0.01, std::numeric_limits<double>::min());
        minimizer->SetLowerLimitedVariable(2, "r", r, 0.01, 1.0);
        minimizer->Minimize();

        double p [] = {alpha, beta, r};
        double logL = L(&p[0]);

        if (logL > minimizer->MinValue() && minimizer->X()[0] > 1 && minimizer->X()[1] > 0 && minimizer->X()[2] >= 1)
        {   
            alpha = minimizer->X()[0];
            beta = minimizer->X()[1];
            r = minimizer->X()[2];
        }
    }
}

void ZeroAdjustedBetaNegativeBinomial::precalculate_log_pmf(const double* params, std::vector<double>& logpmf)
{
    BetaNegativeBinomial::precalculate_log_pmf(params, logpmf);

    double link = lp::log1x(-lp::ext_exp(logpmf[0]));
    std::for_each(logpmf.begin(), logpmf.end(), [link](double& p){p = p - link;});
    logpmf[0] = 0.0;
}

Discrete::Discrete(std::vector<double> p): prob(p)
{
    if (fabs(std::accumulate(prob.begin(), prob.end(), 0.0) - 1.0) > std::numeric_limits<double>::epsilon())
        throw std::invalid_argument("Probability vector must sum up to 1.");
}

double Discrete::pmf(int x) const
{
    return prob[x];
}

double Discrete::log_pmf(int x) const
{
    return lp::ext_log(prob[x]);
}

std::vector<double> Discrete::get_parameters() const
{
    return prob;
}

void Discrete::update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{
    double denom = 0.0;
    for (auto it = gammaBegin; it != gammaEnd; ++it)
    {
        denom += *it;
    }
    
    if (denom <= 0.0)
        throw std::logic_error("Membership coefficients are all zero for at least one state.");

    for (size_t i = 0; i < prob.size(); ++i)
    {
        double num = 0.0;
        auto it2 = obsBegin;
        for (auto it = gammaBegin; it != gammaEnd; ++it)
        {
            if (*it2 == i)
            {
                num += *it;
            }
            ++it2;
        }
        prob[i] = num / denom;
    }
}

Gaussian::Gaussian(double mu, double sigma) : mu(mu), sigma(sigma)
{
    if (sigma <= 0)
        throw std::invalid_argument("Variance of a normal distibution must be > 0.");
    lognorm = -lp::ext_log(sigma * sqrt(2 * M_PI));
}

double Gaussian::pmf(int x) const
{
    return (1.0 / (sigma * sqrt(2 * M_PI))) * lp::ext_exp(-0.5 * pow((x - mu) / sigma, 2));
}

double Gaussian::log_pmf(int x) const
{
    return lognorm - 0.5 * pow((x - mu) / sigma, 2);
}

std::vector<double> Gaussian::get_parameters() const
{
    return std::vector<double>{mu, sigma};
}

void Gaussian::update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{
    mu = stat::sample_mean(gammaBegin, gammaEnd, obsBegin, obsEnd);

    double var = stat::sample_variance(mu, gammaBegin, gammaEnd, obsBegin, obsEnd);
    double newSigma = sqrt(var);
    if (newSigma > 0)
    {
        sigma = newSigma;
        lognorm = -lp::ext_log(sigma * sqrt(2 * M_PI));
    }
}

Sichel::Sichel(double mu, double sigma, double v, int max) : mu(mu), sigma(sigma), v(v), max(max)
{
    if (sigma <= 0.0 || mu <= 0.0)
        throw std::invalid_argument("Mu and Sigma must be > 0.");

    pmf_ = std::vector<double>(max+1);
    precalculate_pmf(mu, sigma, v, pmf_);

    minimizer = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    minimizer->SetMaxFunctionCalls(300);
    minimizer->SetTolerance(1e-16);
}

double Sichel::pmf(int x) const
{
    if (x < pmf_.size())
    {
        return pmf_[x];
    }
    
    double alpha = sqrt(1.0 / pow(sigma, 2) + 2.0 * mu / sigma);
    double w = sqrt(pow(mu, 2) + pow(alpha, 2)) - mu;

    if (x == 0)
    {
        return pow(w / alpha, v) * bmath::cyl_bessel_k(v, alpha) / bmath::cyl_bessel_k(v, w);
    }
    
    if (x == 1)
    {
        double p1 = pow(w / alpha, v) * bmath::cyl_bessel_k(v, alpha) / bmath::cyl_bessel_k(v, w);
        return p1 * (mu * w / alpha) * bmath::cyl_bessel_k(v+1, alpha) / bmath::cyl_bessel_k(v, alpha);
    }

    int i = pmf_.size();
    double p1 = pmf_[i-2];
    double p2 = pmf_[i-1];
    double p = 0.0;
    for (; i <= x; ++i)
    {
        p = (2.0 * mu * w / pow(alpha, 2)) * ((i + v - 1.0) / i) * p2 + pow(mu * w / alpha, 2) / (i * (i - 1.0)) * p1;
        p1 = p2;
        p2 = p;
    }
    return p;
}

double Sichel::log_pmf(int x) const
{
    double p = pmf(x);
    return p > 0.0 ? lp::ext_log(p) : lp::ext_log(std::numeric_limits<double>::min());
}

std::vector<double> Sichel::get_parameters() const
{
    return std::vector<double>{mu, sigma, v};
}

void Sichel::update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{
    using namespace ROOT::Math;

    max = *std::max_element(obsBegin, obsEnd);
    if (pmf_.size() <= max)
    {
        pmf_.resize(max+1);
    }

    LogLikelihoodFunction L (gammaBegin, gammaEnd, obsBegin, obsEnd, pmf_, precalculate_log_pmf);
    Functor f(L, 3);

    try
    {
        minimizer->SetFunction(f);
        minimizer->SetLowerLimitedVariable(0, "mu", mu, 0.1, std::numeric_limits<double>::min());
        minimizer->SetLowerLimitedVariable(1, "sigma", sigma, 0.1, std::numeric_limits<double>::min());
        minimizer->SetVariable(2, "v", v, 0.1);
        minimizer->Minimize();

        double p [] = {mu, sigma, v};
        double logL = L(&p[0]);

        if (logL > minimizer->MinValue() && minimizer->X()[0] > 0 && minimizer->X()[1] > 0)
        {  
            precalculate_pmf(minimizer->X()[0], minimizer->X()[1], minimizer->X()[2], pmf_);

            if (std::find_if(pmf_.begin(), pmf_.end(), [] (double p) {return p > 1.0 || p < 0.0;}) != pmf_.end())
            {
                throw std::logic_error("Illegal probabilities encountered.");
            }
            mu = minimizer->X()[0];
            sigma = minimizer->X()[1];
            v = minimizer->X()[2];
        }
        else
        {
            precalculate_pmf(mu, sigma, v, pmf_);
        }
    }
    catch(const std::exception& e)
    {
        std::clog << "Updating Sichel parameters not successful.\n";
        precalculate_pmf(mu, sigma, v, pmf_);
    }
}

void Sichel::precalculate_pmf(double mu, double sigma, double v, std::vector<double>& pmf)
{
    double alpha = sqrt(1.0 / pow(sigma, 2) + 2.0 * mu / sigma);
    double w = sqrt(pow(mu, 2) + pow(alpha, 2)) - mu;

    for (size_t i = 0; i < pmf.size(); ++i)
    {
        if (i == 0)
        {
            pmf[0] = pow(w / alpha, v) * bmath::cyl_bessel_k(v, alpha) / bmath::cyl_bessel_k(v, w);
        }
        else if (i == 1)
        {
            pmf[1] = pmf[0] * (mu * w / alpha) * bmath::cyl_bessel_k(v+1, alpha) / bmath::cyl_bessel_k(v, alpha);
        }
        else
        {
            pmf[i] = (2.0 * mu * w / pow(alpha, 2)) * ((i + v - 1.0) / i) * pmf[i-1] + pow(mu * w / alpha, 2) / (i * (i - 1.0)) * pmf[i-2];
        }
    }
}

void Sichel::precalculate_log_pmf(const double* params, std::vector<double>& logpmf)
{
    precalculate_pmf(params[0], params[1], params[2], logpmf);
    std::for_each(logpmf.begin(), logpmf.end(), 
                    [](double& p){p = p <= 0 ? lp::ext_log(std::numeric_limits<double>::min()) : lp::ext_log(p);});
}

ZeroAdjustedSichel::ZeroAdjustedSichel(double mu, double sigma, double v, double pi, int max) : mu(mu), sigma(sigma), v(v), pi(pi), max(max)
{
    if (sigma <= 0.0 || mu <= 0.0)
        throw std::invalid_argument("Mu and Sigma must be > 0.");

    if (pi < 0.0 || pi > 1.0)
        throw std::invalid_argument("Zero probability must be in the interval [0, 1].");

    pmf_ = std::vector<double>(max+1);
    Sichel::precalculate_pmf(mu, sigma, v, pmf_);

    minimizer = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    minimizer->SetMaxFunctionCalls(300);
    minimizer->SetTolerance(1e-16);
}

double ZeroAdjustedSichel::pmf(int x) const
{
    if (x == 0)
        return pi;

    if (pi > 1.0 - std::numeric_limits<double>::epsilon())
        return 0.0;
    
    if (x < pmf_.size())
    {
        return (1.0 - pi) * pmf_[x] / (1.0 - pmf_[0]);
    }
    
    double alpha = sqrt(1.0 / pow(sigma, 2) + 2.0 * mu / sigma);
    double w = sqrt(pow(mu, 2) + pow(alpha, 2)) - mu;
    double p0 = pow(w / alpha, v) * bmath::cyl_bessel_k(v, alpha) / bmath::cyl_bessel_k(v, w);
    
    if (x == 1)
    {
        double p1 = pow(w / alpha, v) * bmath::cyl_bessel_k(v, alpha) / bmath::cyl_bessel_k(v, w);
        double res = p1 * (mu * w / alpha) * bmath::cyl_bessel_k(v+1, alpha) / bmath::cyl_bessel_k(v, alpha);
        return (1.0 - pi) * res / (1.0 - p0);
    }

    int i = pmf_.size();
    double p1 = pmf_[i-2];
    double p2 = pmf_[i-1];
    double p = 0.0;
    for (; i <= x; ++i)
    {
        p = (2.0 * mu * w / pow(alpha, 2)) * ((i + v - 1.0) / i) * p2 + pow(mu * w / alpha, 2) / (i * (i - 1.0)) * p1;
        p1 = p2;
        p2 = p;
    }
    return (1.0 - pi) * p / (1.0 - p0);
}

double ZeroAdjustedSichel::log_pmf(int x) const
{
    if (x == 0)
        return lp::ext_log(pi);

    if (pi > 1.0 - std::numeric_limits<double>::epsilon())
        return NAN;

    double p = pmf(x);
    return p > 0.0 ? lp::ext_log(p) : lp::ext_log(std::numeric_limits<double>::min());
}

std::vector<double> ZeroAdjustedSichel::get_parameters() const
{
    return std::vector<double>{mu, sigma, v, pi};
}

void ZeroAdjustedSichel::update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{
    pi = stat::zero_frequency(gammaBegin, gammaEnd, obsBegin, obsEnd);
    if (pi > 1.0 - std::numeric_limits<double>::epsilon())
    {
        pi = 1.0;
        v = -0.5;
    }
    else
    {
        using namespace ROOT::Math;

        max = *std::max_element(obsBegin, obsEnd);
        if (pmf_.size() <= max)
        {
            pmf_.resize(max+1);
        }

        LogLikelihoodFunction L (gammaBegin, gammaEnd, obsBegin, obsEnd, pmf_, precalculate_log_pmf);
        Functor f(L, 3);

        try
        {
            minimizer->SetFunction(f);
            minimizer->SetLowerLimitedVariable(0, "mu", mu, 0.1, std::numeric_limits<double>::min());
            minimizer->SetLowerLimitedVariable(1, "sigma", sigma, 0.1, std::numeric_limits<double>::min());
            minimizer->SetVariable(2, "v", v, 0.1);
            minimizer->Minimize();

            double p [] = {mu, sigma, v};
            double logL = L(&p[0]);

            if (logL > minimizer->MinValue() && minimizer->X()[0] > 0 && minimizer->X()[1] > 0)
            {   
                Sichel::precalculate_pmf(minimizer->X()[0], minimizer->X()[1], minimizer->X()[2], pmf_);

                if (std::find_if(pmf_.begin(), pmf_.end(), [] (double p) {return p > 1.0 || p < 0.0;}) != pmf_.end())
                {
                    throw std::logic_error("Illegal probabilities encountered.");
                }

                mu = minimizer->X()[0];
                sigma = minimizer->X()[1];
                v = minimizer->X()[2];

            }
            else
            {
                Sichel::precalculate_pmf(mu, sigma, v, pmf_);
            }

        }
        catch(const std::exception& e)
        {
            std::clog << "Updating Sichel parameters not successful.\n";
            Sichel::precalculate_pmf(mu, sigma, v, pmf_);
        }
    }
}

void ZeroAdjustedSichel::precalculate_log_pmf(const double* params, std::vector<double>& logpmf)
{
    Sichel::precalculate_log_pmf(params, logpmf);

    double link = lp::log1x(-lp::ext_exp(logpmf[0]));
    std::for_each(logpmf.begin(), logpmf.end(), [link](double& p){p = p - link;});
    logpmf[0] = 0.0;
}

AdjustedBeta::AdjustedBeta(double a, double b, double p): alpha(a), beta(b), pi(p)
{   
    if (alpha <= 0.0 || beta <= 0.0)
        throw std::invalid_argument("Both alpha and beta of the beta distribution must be > 0.");

    if (pi < 0.0 || pi > 1.0)
        throw std::invalid_argument("Zero probability must be between 0 and 1");
}

double AdjustedBeta::pmf(int n, int x) const
{
    if (n < 0)
        throw std::invalid_argument("Number of trials must be > 0 for the beta distribution.");

    if (x > n)
        return 0.0;

    if (x == 0 && n == 0)
        return pi;

    double lower = std::floor((double)x / n * 1000.0) / 1000.0;
    if (lower > 0.999)
    {
        return (1.0 - pi) * (1.0 - bmath::ibeta(alpha, beta, 0.999));
    }

    double upper = lower + 0.001;
    return (1.0 - pi) * (bmath::ibeta(alpha, beta, upper) - bmath::ibeta(alpha, beta, lower));
}

double AdjustedBeta::log_pmf(int n, int x) const
{   
    return lp::ext_log(pmf(n, x));
}

std::vector<double> AdjustedBeta::get_parameters() const
{
    return std::vector<double>{alpha, beta, pi};
}

void AdjustedBeta::update_methylation(dna_g_iterator gammaBegin, dna_g_iterator gammaEnd, dna_o_iterator covBegin, dna_o_iterator methBegin)
{
    // method of moment estimates
    pi = 0.0;
    double denom = 0.0;
    double mean = 0.0;
    auto covIt = covBegin;
    auto methIt = methBegin;
    for (auto it = gammaBegin; it != gammaEnd; ++it)
    {
        if (*covIt == 0 && *methIt == 0)
        {
            pi += *it;
        }
        else
        {
            mean += (*it) * ((double) *methIt / *covIt);
        }
        denom += *it;
        ++covIt;
        ++methIt;
    }
    pi /= denom;
    if (pi >= 1.0 - std::pow(10, -12))
    {
        pi = 1.0;
        alpha = 1.0; beta = 1.0;
        return; 
    }

    mean /= denom;
    double var = 0.0;
    covIt = covBegin;
    methIt = methBegin;
    for (auto it = gammaBegin; it != gammaEnd; ++it)
    {
        int cov = *covIt;
        int meth = *methIt;
        if (!(cov == 0 && meth == 0))
        {
            double p = (double) meth / cov;
            var += (*it) * std::pow(p - mean, 2);
        }
        ++covIt;
        ++methIt;
    }
    var /= denom;

    double psi = (mean * (1-mean) / var) - 1.0;
    if (psi <= 0.0)
    {
        std::clog << "Adjusted beta distribution: alpha or beta became negative, no update performed! \n";
        return;
    }
    alpha = mean * psi;
    beta = (1.0 - mean) * psi;
}

Bernoulli::Bernoulli(double p)
{
    if (p < 0.0 || p > 1.0)
    {
        std::cout << p << std::endl;
        throw std::invalid_argument("Parameter p of a bernoulli distribution must be between 0 and 1.");
    }

    logp = lp::ext_log(p);
    log1_p = lp::log1x(-p);
}

double Bernoulli::pmf(int x) const
{
    return lp::ext_exp(log_pmf(x));
}

double Bernoulli::log_pmf(int x) const
{
    if (x == 1)
        return logp;
    if (x == 0)
        return log1_p;
    else
        throw std::invalid_argument("Bernoulli distribution is only valid for binarized data.");
}

std::vector<double> Bernoulli::get_parameters() const
{
    return std::vector<double>{lp::ext_exp(logp)};
}

void Bernoulli::update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd)
{    
    double sampleMean = stat::sample_mean(gammaBegin, gammaEnd, obsBegin, obsEnd);
    if (sampleMean >= 0.0 && sampleMean <= 1.0)
    {
        logp = lp::ext_log(sampleMean);
        log1_p = lp::log1x(-sampleMean);
    }
}