#ifndef _DISTRIBUTION_H
#define _DISTRIBUTION_H

#include <vector>
#include <numeric>
#include <tuple>
#include "matrix.h"
#include "Math/Minimizer.h"

/**
 * @brief Abstract base distribution all available distributions have to inherit from and override the required functionalities.
 */
class DiscreteDistribution
{
    public:
    using g_iterator = std::vector<double>::iterator;
    using o_iterator = Matrix<int>::const_col_iterator;

    DiscreteDistribution() {}; 
    /**
     * @brief Returns the probability mass function for x using the current parameters of the model.
     * 
     * @param x 
     * @return P(X = x) 
     */
    virtual double pmf(int x) const = 0;

    /**
     * @brief Returns the logarithm of the probability mass function for x using the current parameters of the model.
     * 
     * @param x 
     * @return log(P(X = x))
     */
    virtual double log_pmf(int x) const = 0;

    /**
     * @brief updates the parameters of the distribution based on the membership coefficients gamma and the current observation.
     * 
     * @param gamma_begin 
     * @param gamma_end 
     * @param obs_begin 
     * @param obs_end 
     */
    virtual void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd) = 0;

    /**
     * @brief Returns all parameters of the distribution.
     * 
     * @return std::vector<double> 
     */
    virtual std::vector<double> get_parameters() const = 0;

    /**
     * @brief Returns the name of the distribution
     * 
     * @return std::string 
     */
    virtual std::string get_name() const = 0;
};

/**
 * @brief Abstract base distribution for two value observations (trials and successes), e.g. useful for DNA methylation data
 */
class TwoValueDiscreteDistribution
{
    public:
    using dna_g_iterator = std::vector<double>::iterator;
    using dna_o_iterator = Matrix<int>::const_col_iterator;

    TwoValueDiscreteDistribution() {}; 
    /**
     * @brief Returns the probability mass function for n trials and x successes using the current parameters of the model.
     * 
     * @param x 
     * @return P(X = x) 
     */
    virtual double pmf(int n, int x) const = 0;

    /**
     * @brief Returns the logarithm of the probability mass function for n trials and x successes using the current parameters of the model.
     * 
     * @param n
     * @param x 
     * @return log(P(X = x))
     */
    virtual double log_pmf(int n, int x) const = 0;

    /**
     * @brief updates the parameters of the distribution based on the membership coefficients gamma and the current observation.
     * 
     * @param gamma_begin 
     * @param gamma_end 
     * @param obs_begin 
     * @param obs_end 
     */
    virtual void update_methylation(dna_g_iterator gammaBegin, dna_g_iterator gammaEnd, dna_o_iterator covBegin, dna_o_iterator methBegin) = 0;
};

class Poisson : public DiscreteDistribution
{
    public:
    /**
     * @brief Initializes a Poisson distribution with parameter lambda.
     * 
     * @param l 
     */
    Poisson(double l);

    double pmf(int x) const;

    double log_pmf(int x) const;

    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "PO";};

    private:
    double lambda;
};

class ZeroAdjustedPoisson : public DiscreteDistribution
{
    public:
    /**
     * @brief Initializes a zero adjusted Poisson distribution with parameter lambda and zero probability pi.
     * 
     * @param l 
     * @param pi
     */
    ZeroAdjustedPoisson(double l, double pi);

    double pmf(int x) const;

    double log_pmf(int x) const;

    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "ZAP";};

    static void precalculate_log_pmf(const double* params, std::vector<double>& logpmf);

    private:
    double lambda;
    double pi;

    std::unique_ptr<ROOT::Math::Minimizer> minimizer;  
};

class Binomial : public DiscreteDistribution, public TwoValueDiscreteDistribution
{
    public:
    /**
     * @brief Initializes a Binomial distribution with success parameter p and trials n.
     * The parameter n is assumed to be known thus not updated (e.g. by setting n to the maximum value)
     * 
     * @param p 
     * @param n
     */
    Binomial(double p, int n);

    /**
     * @brief Initializes a Binomial distribution with success parameter p assuming trials and successes as input for pmf and logpmf.
     * 
     * @param p 
     */
    Binomial(double p);

    double pmf(int x) const;

    double log_pmf(int x) const;

    double pmf(int n, int x) const;

    double log_pmf(int n, int x) const;

    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    void update_methylation(dna_g_iterator gammaBegin, dna_g_iterator gammaEnd, dna_o_iterator covBegin, dna_o_iterator methBeging);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "BI";};

    private:
    double p;
    int n;
    double lgammaN;
};

class NegativeBinomial : public DiscreteDistribution
{
    public:
    /**
     * @brief Initializes a negative binomial distribution with success probability p and number of successes r.
     * 
     * @param p 
     * @param r 
     */
    NegativeBinomial(double p, double r);

    double pmf(int x) const;

    double log_pmf(int x) const;

    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "NBI";};

    static void precalculate_log_pmf(const double* params, std::vector<double>& logpmf);

    private:

    void find_mle(double initP, double initR, g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    double p;
    double r;
    double lgammaR;

    std::unique_ptr<ROOT::Math::Minimizer> minimizer;    
};

class ZeroAdjustedNegativeBinomial : public DiscreteDistribution
{
    public:
    /**
     * @brief Initializes a zero adjusted negative binomial distribution with success probability p, number of successes r and zero probability pi.
     * 
     * @param p 
     * @param r 
     * @param pi
     */
    ZeroAdjustedNegativeBinomial(double p, double r, double pi);

    double pmf(int x) const;

    double log_pmf(int x) const;

    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "ZANBI";};

    static void precalculate_log_pmf(const double* params, std::vector<double>& logpmf);

    private:

    double p;
    double r;
    double pi;
    double lgammaR;

    std::unique_ptr<ROOT::Math::Minimizer> minimizer;
};


class BetaBinomial : public DiscreteDistribution, public TwoValueDiscreteDistribution
{
    public:
    /**
     * @brief Initializes a beta binomial distribution with parameters alpha, beta, n
     * The parameter n is assumed to be known and not updated (e.g., use sample maximum as an approximation for the real n)
     * 
     * @param alpha
     * @param beta
     * @param n
     */
    BetaBinomial(double alpha, double beta, int n);

    /**
     * @brief Initializes a beta binomial distribution with parameters alpha and beta assuming trials and successes as input for pmf and logpmf.
     * 
     * @param alpha
     * @param beta
     */
    BetaBinomial(double alpha, double beta);

    double pmf(int x) const;

    double log_pmf(int x) const;

    double pmf(int n, int x) const;

    double log_pmf(int n, int x) const;

    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    void update_methylation(dna_g_iterator gammaBegin, dna_g_iterator gammaEnd, dna_o_iterator covBegin, dna_o_iterator methBegin);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "BB";};

    static void precalculate_log_pmf(const double* params, int n, std::vector<double>& pmf);

    static double calculate_log_pmf(const double* params, int n, int x);

    private:
    void find_mle(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd, double initA, double initB);

    double alpha;
    double beta;
    int n;

    double lgammaN1;
    double lgammaNAB;
    double lgammaAB;
    double lgammaA;
    double lgammaB;

    std::unique_ptr<ROOT::Math::Minimizer> minimizer;
};

class BetaNegativeBinomial : public DiscreteDistribution
{
    public:
    /**
     * @brief Initializes a beta negative binomial distribution with parameters alpha, beta, r
     * 
     * @param alpha
     * @param beta
     * @param r
     */
    BetaNegativeBinomial(double alpha, double beta, double r);

    double pmf(int x) const;

    double log_pmf(int x) const;

    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "BNB";};

    static void precalculate_log_pmf(const double* params, std::vector<double>& pmf);

    private:
    void find_mle(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    double alpha;
    double beta;
    double r;

    double lgammaR;
    double lgammaAR;
    double lgammaAB;
    double lgammaA;
    double lgammaB;

    std::unique_ptr<ROOT::Math::Minimizer> minimizer;
};

class ZeroAdjustedBetaNegativeBinomial : public DiscreteDistribution
{
    public:
    /**
     * @brief Initializes a zero adjusted beta negative binomial distribution with parameters alpha, beta, r and zero probability pi
     * 
     * @param alpha
     * @param beta
     * @param r
     * @param pi
     */
    ZeroAdjustedBetaNegativeBinomial(double alpha, double beta, double r, double pi);

    double pmf(int x) const;

    double log_pmf(int x) const;

    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "ZABNB";};

    static void precalculate_log_pmf(const double* params, std::vector<double>& logpmf);

    private:
    void find_mle(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    double alpha;
    double beta;
    double r;
    double pi;

    double lgammaR;
    double lgammaAR;
    double lgammaAB;
    double lgammaA;
    double lgammaB;

    std::unique_ptr<ROOT::Math::Minimizer> minimizer;
};


class Sichel : public DiscreteDistribution
{
    public:
    /**
     * @brief Initializes a Sichel distribution with parameters mu > 0, sigma > 0, -infty < v < infty and precalculates probabilities up to max
     * 
     * @param mu
     * @param sigma
     * @param v
     * @param max
     */
    Sichel(double mu, double sigma, double v, int max = 1000);

    double pmf(int x) const;

    double log_pmf(int x) const;

    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "SI";};

    static void precalculate_pmf(double mu, double sigma, double v, std::vector<double>& pmf);

    static void precalculate_log_pmf(const double * params, std::vector<double>& logpmf);

    private:

    double mu;
    double sigma;
    double v;
    int max;

    std::vector<double> pmf_;
    std::unique_ptr<ROOT::Math::Minimizer> minimizer;
};

class ZeroAdjustedSichel : public DiscreteDistribution
{
    public:
    /**
     * @brief Initializes a zero adjusted Sichel distribution with parameters mu > 0, sigma > 0, -infty < v < infty and zero probability pi (probabilities are precalculted up to max)
     * 
     * @param mu
     * @param sigma
     * @param v 
     * @param pi
     * @param max
     */
    ZeroAdjustedSichel(double mu, double sigma, double v, double pi, int max = 1000);

    double pmf(int x) const;

    double log_pmf(int x) const;

    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "ZASI";};

    static void precalculate_log_pmf(const double* params, std::vector<double>& logpmf);

    private:

    double mu;
    double sigma;
    double v;
    double pi;
    int max;

    std::vector<double> pmf_;

    std::unique_ptr<ROOT::Math::Minimizer> minimizer;
};

class Discrete : public DiscreteDistribution
{
    public:
    /**
     * @brief Initializes a discrete distribution with probability vector p (with values 0, ..., len(p))
     * 
     * @param p 
     */
    Discrete(std::vector<double> p);

    double pmf(int x) const;

    double log_pmf(int x) const;

    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "DIS";};

    private:
    std::vector<double> prob;
};


class Gaussian : public DiscreteDistribution
{
    public:
    /**
     * @brief Initializes a gaussian distribution with mean mu and variance sigma^2
     * 
     * @param mu 
     * @param sigma 
     */
    Gaussian(double mu, double sigma);

    double pmf(int x) const;

    double log_pmf(int x) const;

    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "GA";};

    private:
    double mu;
    double sigma;
    double lognorm;
};

class Bernoulli : public DiscreteDistribution
{
    public:
    /**
     * @brief Initializes a bernoulli distribution with parameter p
     * 
     * @param p  
     */
    Bernoulli(double p);

    double pmf(int x) const;

    double log_pmf(int x) const;

    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "B";};

    private:
    double logp;
    double log1_p;
};

class AdjustedBeta : public DiscreteDistribution, public TwoValueDiscreteDistribution
{
    public:
    /**
     * @brief Initializes a Beta distribution with parameters alpha and beta and probability pi for zero coverage.
     * 
     * @param alpha
     * @param beta
     */
    AdjustedBeta(double alpha, double beta, double pi);

    double pmf(int n, int x) const;

    double log_pmf(int n, int x) const;

    void update_methylation(dna_g_iterator gammaBegin, dna_g_iterator gammaEnd, dna_o_iterator covBegin, dna_o_iterator methBegin);

    std::vector<double> get_parameters() const;

    std::string get_name() const {return "AB";};

    double pmf(int x) const {throw std::invalid_argument("Only two column input valid for beta distribution."); return 0.0;};
    double log_pmf(int x) const {throw std::invalid_argument("Only two column input valid for beta distribution."); return 0.0;};
    void update(g_iterator gammaBegin, g_iterator gammaEnd, o_iterator obsBegin, o_iterator obsEnd) 
        {throw std::invalid_argument("Only two column input valid for beta distribution.");};

    private:
    double alpha;
    double beta;
    double pi;
};

#endif