#ifndef _LOG_PROB_H
#define _LOG_PROB_H

#include <vector>
#include <numeric>
#include <math.h>
#include <cmath>

namespace lp
{
    /**
     * @brief Standard exponential function but can handle ln(0)==nan.
     * 
     * @param x 
     * @return exp(x)
     */
    double ext_exp(double x);

    /**
     * @brief Standard logarithm but returns nan for ln(0).
     * 
     * @param x 
     * @return log(x)
     */
    double ext_log(double x);

    /**
     * @brief Returns the logarithm of the sum of x and y, expecting ln(x) and ln(y) as input.
    Using the following formula for numerical stability: ln(x) + ln(1+e^(ln(y) - ln(x))), where ln(x) > ln(y)
     * 
     * @param logx 
     * @param logy 
     * @return log(x+y)
     */
    double log_add(double logx, double logy);

    /**
     * @brief Returns the logarithm of the product of x and y, expecting ln(x) and ln(y) as input.
     * 
     * @param logx 
     * @param logy 
     * @return log(x * y)
     */
    double log_mul(double logx, double logy);

    /**
     * @brief Comparison operator for logarithmic probabilities (e.g. NAN = log(0) < log(0.5) = -0.69)
     * 
     * @param logx 
     * @param logy 
     * @return log(x) > log(y)
     */
    bool log_greater(double logx, double logy);

    /**
     * @brief calculate x * log(y), with handling the case of x = 0 (return 0)
     * 
     * @param x 
     * @param y 
     * @return double 
     */
    double xlogy(double x, double y);

    /**
     * @brief calculate log(1+x), with handling the case of x = -1 (return 0)
     * 
     * @param x 
     * @return double 
     */
    double log1x(double x);

    /**
     * @brief Computes log(x1 + ... + xn) assuming an iterator to the first and last position of [log(x1), ..., log(xn)] as input.
     * 
     * @param Iterator begin, Iterator end 
     * @return log(x1 + ... + xn)
     */
    template<typename Iterator>
    double log_sum(Iterator begin, Iterator end);

    /**
     * @brief Computes log(x1 * ... * xn) assuming an iterator to the first and last position of [log(x1), ..., log(xn)] as input.
     * 
     * @param Iterator begin, Iterator end
     * @return log(x1 * ... * xn)
     */
    template<typename Iterator>
    double log_prod(Iterator begin, Iterator end);
};

template<typename Iterator>
double lp::log_sum(Iterator begin, Iterator end)
{
    double max = std::numeric_limits<int>::lowest();
    int maxIndex = -1;
    int i = 0;
    for (Iterator it = begin; it != end; ++it)
    {
        if (!isnan(*it) && *it > max)
        {
            max = *it;
            maxIndex = i;
        }
        ++i;
    }

    if (maxIndex == -1)
    {
        return NAN;
    }

    double res = 0.0;
    i = 0;
    for (Iterator it = begin; it != end; ++it)
    {
        if (i != maxIndex)
        {
            res += ext_exp(*it - max);
        }
        ++i;
    }

    return max + log1p(res);
}

template<typename Iterator>
inline double lp::log_prod(Iterator begin, Iterator end)
{
    return std::accumulate(begin, end, 0.0);
}

inline double lp::ext_exp(double x)
{
    if (isnan(x)) 
        return 0;
    return exp(x);
}

inline double lp::ext_log(double x)
{
    if (x <= 0)
        return NAN;
    return log(x);
}

inline double lp::log_add(double logx, double logy)
{
    if (isnan(logx))
        return logy;

    if (isnan(logy))
        return logx;
    
    if (logx > logy)
        return logx + log1p(ext_exp(logy-logx));
    
    return logy + log1p(ext_exp(logx-logy));
}

inline double lp::log_mul(double logx, double logy)
{
    return logx + logy;
}

inline bool lp::log_greater(double logx, double logy)
{
    if (isnan(logx) && isnan(logy))
        return false;

    if (isnan(logy))
        return true;
    
    if (isnan(logx))
        return false;

    return logx > logy;
}

inline double lp::xlogy(double x, double y)
{
    if (x <= 0)
        return 0;
    return x * ext_log(y);
}

inline double lp::log1x(double x)
{
    if (1+x <= 0)
        return NAN;
    return log1p(x);
}

#endif