#ifndef _SAMPLE_STATISTICS_H
#define _SAMPLE_STATISTICS_H

namespace stat 
{
    /**
     * @brief calculates the membership weighted sample mean
     * 
     * @param gammaBegin 
     * @param gammaEnd 
     * @param obsBegin 
     * @param obsEnd 
     * @return mean 
     */
    template<typename It1, typename It2>
    double sample_mean(It1 gammaBegin, It1 gammaEnd, It2 obsBegin, It2 obsEnd);

    /**
     * @brief calculates the membership weighted sample variance
     * 
     * @param mean 
     * @param gammaBegin 
     * @param gammaEnd 
     * @param obsBegin 
     * @param obsEnd 
     * @return variance 
     */
    template<typename It1, typename It2>
    double sample_variance(double mean, It1 gammaBegin, It1 gammaEnd, It2 obsBegin, It2 obsEnd);

    /**
     * @brief calculates the membership weighted zero frequency
     * 
     * @param gammaBegin 
     * @param gammaEnd 
     * @param obsBegin 
     * @param obsEnd 
     * @return zero_frequency
     */
    template<typename It1, typename It2>
    double zero_frequency(It1 gammaBegin, It1 gammaEnd, It2 obsBegin, It2 obsEnd);

    /**
     * @brief calculates the membership weighted mean ignoring zero values
     * 
     * @param gammaBegin 
     * @param gammaEnd 
     * @param obsBegin 
     * @param obsEnd 
     * @return mean
     */
    template<typename It1, typename It2>
    double sample_mean_without_zeros(It1 gammaBegin, It1 gammaEnd, It2 obsBegin, It2 obsEnd);

    /**
     * @brief calculates the membership weighted variance ignoring zero values
     * 
     * @param gammaBegin 
     * @param gammaEnd 
     * @param obsBegin 
     * @param obsEnd 
     * @return mean
     */
    template<typename It1, typename It2>
    double sample_variance_without_zeros(double mean, It1 gammaBegin, It1 gammaEnd, It2 obsBegin, It2 obsEnd);
}

template<typename It1, typename It2>
double stat::sample_mean(It1 gammaBegin, It1 gammaEnd, It2 obsBegin, It2 obsEnd)
{
    double denom = 0.0;
    double num = 0.0;

    auto obsIt = obsBegin;
    for (auto it = gammaBegin; it != gammaEnd; ++it)
    {
        denom += *it;
        num += (*it) * (*obsIt);
        ++obsIt;
    }

    if (denom <= 0.0)
        throw std::logic_error("Membership coefficients are all zero for at least one state.");

    return num / denom;
}

template<typename It1, typename It2>
double stat::sample_variance(double mean, It1 gammaBegin, It1 gammaEnd, It2 obsBegin, It2 obsEnd)
{
    double num = 0.0;
    double denom = 0.0;

    auto obsIt = obsBegin;
    for (auto it = gammaBegin; it != gammaEnd; ++it)
    {
        num += (*it) * std::pow(*obsIt - mean, 2);
        denom += *it;
        ++obsIt;
    }

    return num / denom;
}

template<typename It1, typename It2>
double stat::zero_frequency(It1 gammaBegin, It1 gammaEnd, It2 obsBegin, It2 obsEnd)
{
    double num = 0.0;
    double denom = 0.0;

    auto obsIt = obsBegin;
    for (auto it = gammaBegin; it != gammaEnd; ++it)
    {
        if (*obsIt == 0)
        {
            num += *it;
        }
        denom += *it;
        ++obsIt;
    }

    if (denom <= 0.0)
        throw std::logic_error("Membership coefficients are all zero for at least one state.");
    return num / denom;
}

template<typename It1, typename It2>
double stat::sample_mean_without_zeros(It1 gammaBegin, It1 gammaEnd, It2 obsBegin, It2 obsEnd)
{
    double denom = 0.0;
    double num = 0.0;

    auto obsIt = obsBegin;
    for (auto it = gammaBegin; it != gammaEnd; ++it)
    {
        if (*obsIt != 0)
        {
            denom += *it;
            num += (*it) * (*obsIt);   
        }
        ++obsIt;
    }

    if (denom <= 0.0)
        throw std::logic_error("Membership coefficients are all zero for at least one state.");

    return num / denom;
}

template<typename It1, typename It2>
double stat::sample_variance_without_zeros(double mean, It1 gammaBegin, It1 gammaEnd, It2 obsBegin, It2 obsEnd)
{
    double num = 0.0;
    double denom = 0.0;

    auto obsIt = obsBegin;
    for (auto it = gammaBegin; it != gammaEnd; ++it)
    {
        if (*obsIt != 0)
        {
            num += (*it) * std::pow(*obsIt - mean, 2);
            denom += *it;
        }
        ++obsIt;
    }
    return num / denom;
}

#endif