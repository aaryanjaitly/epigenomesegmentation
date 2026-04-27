#include <gtest/gtest.h>
#include <cmath>
#include <math.h>

#include "../distribution.h"

double EPS = std::pow(10, -8);

TEST(Poisson, pmf)
{
    Poisson dis = Poisson(3.0);
    
    EXPECT_NEAR(0.049787068, dis.pmf(0), EPS);
    EXPECT_NEAR(0.021604031, dis.pmf(7), EPS);
    EXPECT_NEAR(0.2240418076, dis.pmf(3), EPS);

    EXPECT_NEAR(-3.0, dis.log_pmf(0), EPS);
    EXPECT_NEAR(-3.83487534, dis.log_pmf(7), EPS);
    EXPECT_NEAR(-1.495922603, dis.log_pmf(3), EPS);

    dis = Poisson(7.5);
    EXPECT_NEAR(0.0005530843, dis.pmf(0), EPS);
    EXPECT_NEAR(0.1464838321, dis.pmf(7), EPS);
    EXPECT_NEAR(0.03888874, dis.pmf(3), EPS);

    EXPECT_NEAR(-7.5, dis.log_pmf(0), EPS);
    EXPECT_NEAR(-1.920840217, dis.log_pmf(7), EPS);
    EXPECT_NEAR(-3.247050407, dis.log_pmf(3), EPS);
}

TEST(Binomial, pmf)
{
    Binomial dis = Binomial(0.5, 20);
    
    EXPECT_NEAR(0.176197052, dis.pmf(10), EPS);
    EXPECT_NEAR(0.0, dis.pmf(21), EPS);
    EXPECT_NEAR(0.001087188, dis.pmf(17), EPS);
    EXPECT_NEAR(9.5 * pow(10, -7), dis.pmf(0), EPS);

    EXPECT_NEAR(-1.736152296, dis.log_pmf(10), EPS);
    EXPECT_TRUE(isnan(dis.log_pmf(21)));
    EXPECT_NEAR(-6.8241600698, dis.log_pmf(17), EPS);
    EXPECT_NEAR(-13.862943611, dis.log_pmf(0), EPS);

    dis = Binomial(0.2, 20);
    EXPECT_NEAR(-4.4628710262, dis.log_pmf(0), EPS);
    EXPECT_NEAR(-2.9086403075, dis.log_pmf(7), EPS);

    dis = Binomial(1.0, 10);
    EXPECT_NEAR(0.0, dis.pmf(7), EPS);
    EXPECT_NEAR(1.0, dis.pmf(10), EPS);
    EXPECT_TRUE(isnan(dis.log_pmf(7)));
    EXPECT_NEAR(0.0, dis.log_pmf(10), EPS);

    dis = Binomial(0.0, 10);
    EXPECT_NEAR(0.0, dis.pmf(7), EPS);
    EXPECT_NEAR(1.0, dis.pmf(0), EPS);
    EXPECT_TRUE(isnan(dis.log_pmf(7)));
    EXPECT_NEAR(0.0, dis.log_pmf(0), EPS);

    dis = Binomial(0.5);
    EXPECT_NEAR(-6.8241600698, dis.log_pmf(20, 17), EPS);
    EXPECT_NEAR(-13.862943611, dis.log_pmf(20, 0), EPS);

    dis = Binomial(1.0);
    EXPECT_NEAR(0.0, dis.pmf(10, 7), EPS);
    EXPECT_NEAR(1.0, dis.pmf(10, 10), EPS);

    dis = Binomial(0.0);
    EXPECT_NEAR(0.0, dis.pmf(10, 7), EPS);
    EXPECT_NEAR(1.0, dis.pmf(10, 0), EPS);
}

TEST(Gaussian, pmf)
{
    Gaussian dis = Gaussian(10, 5);
    
    EXPECT_NEAR(0.079788456, dis.pmf(10), EPS);
    EXPECT_NEAR(0.0070949185, dis.pmf(21), EPS);
    EXPECT_NEAR(0.0299454931, dis.pmf(17), EPS);
    EXPECT_NEAR(0.0107981933, dis.pmf(0), EPS);

    EXPECT_NEAR(-2.5283764456, dis.log_pmf(10), EPS);
    EXPECT_NEAR(-4.948376445, dis.log_pmf(21), EPS);
    EXPECT_NEAR(-3.50837644, dis.log_pmf(17), EPS);
    EXPECT_NEAR(-4.52837644, dis.log_pmf(0), EPS);
}

TEST(NegativeBinomial, pmf)
{
    NegativeBinomial dis = NegativeBinomial(0.5, 10);
    
    EXPECT_NEAR(0.06109619140, dis.pmf(5), EPS);
    EXPECT_NEAR(0.0009765625, dis.pmf(0), EPS);
    EXPECT_NEAR(0.070078372, dis.pmf(12), EPS);

    EXPECT_NEAR(-2.795305748, dis.log_pmf(5), EPS);
    EXPECT_NEAR(-6.93147180, dis.log_pmf(0), EPS);
    EXPECT_NEAR(-2.658141049, dis.log_pmf(12), EPS);

    dis = NegativeBinomial(1.0, 2);
    EXPECT_NEAR(0.0, dis.pmf(1), EPS);
    EXPECT_NEAR(1.0, dis.pmf(0), EPS);

    EXPECT_TRUE(isnan(dis.log_pmf(1)));
    EXPECT_NEAR(0.0, dis.log_pmf(0), EPS);

    dis = NegativeBinomial(0.0, 2);
    EXPECT_NEAR(0.0, dis.pmf(0), EPS);
    EXPECT_TRUE(isnan(dis.log_pmf(1)));
}

TEST(BetaBinomial, pmf)
{
    BetaBinomial dis = BetaBinomial(0.2, 0.25, 10);
    
    EXPECT_NEAR(0.3435192871, dis.pmf(0), EPS);
    EXPECT_NEAR(0.035248498, dis.pmf(4), EPS);
    EXPECT_NEAR(0.0451332432, dis.pmf(8), EPS);
    EXPECT_NEAR(0.2421227996, dis.pmf(10), EPS);

    EXPECT_NEAR(-1.0685120199, dis.log_pmf(0), EPS);
    EXPECT_NEAR(-3.345332359, dis.log_pmf(4), EPS);
    EXPECT_NEAR(-3.098136203, dis.log_pmf(8), EPS);
    EXPECT_NEAR(-1.418310244, dis.log_pmf(10), EPS);

    EXPECT_TRUE(isnan(dis.log_pmf(11)));
    EXPECT_NEAR(0.0, dis.pmf(11), EPS);

    dis = BetaBinomial(0.7, 2, 10);
    EXPECT_NEAR(0.2736474225, dis.pmf(0), EPS);
    EXPECT_NEAR(0.086257710, dis.pmf(4), EPS);
    EXPECT_NEAR(0.0304123583, dis.pmf(8), EPS);
    EXPECT_NEAR(0.0095055515, dis.pmf(10), EPS);

    dis = BetaBinomial(600, 400, 10);

    EXPECT_NEAR(-9.096145715, dis.log_pmf(0), EPS);
    EXPECT_NEAR(-2.19148762, dis.log_pmf(4), EPS);
    EXPECT_NEAR(-2.108411831, dis.log_pmf(8), EPS);
    EXPECT_NEAR(-5.0785071464, dis.log_pmf(10), EPS);

    dis = BetaBinomial(0.2, 0.25);
    EXPECT_NEAR(0.3435192871, dis.pmf(10, 0), EPS);
    EXPECT_NEAR(0.035248498, dis.pmf(10, 4), EPS);
    EXPECT_NEAR(0.2421227996, dis.pmf(10, 10), EPS);

    dis = BetaBinomial(600, 400);
    EXPECT_NEAR(-2.19148762, dis.log_pmf(10, 4), EPS);
}

TEST(Sichel, pmf)
{
    Sichel dis = Sichel(80, 2, 0.5, 10);
    
    EXPECT_NEAR(0.00001184083, dis.pmf(0), EPS);
    EXPECT_NEAR(0.00082845141, dis.pmf(5), EPS);
    EXPECT_NEAR(0.002318059, dis.pmf(10), EPS);
    EXPECT_NEAR(0.002575507, dis.pmf(11), EPS);

    dis = Sichel(50, 10, -0.8, 5);

    EXPECT_NEAR(0.08020148, dis.pmf(0), EPS);
    EXPECT_NEAR(0.05927062, dis.pmf(5), EPS);
    EXPECT_NEAR(0.11180423, dis.pmf(2), EPS);
    EXPECT_NEAR(0.01206296, dis.pmf(15), EPS);
    EXPECT_NEAR(0.00000269, dis.pmf(1000), EPS);   
}

TEST(Beta, pmf)
{
    AdjustedBeta dis = AdjustedBeta(2, 2, 0.2);
    EXPECT_NEAR(0.8 * 0.000542398, dis.pmf(10, 1), EPS);
    EXPECT_NEAR(0.2, dis.pmf(0, 0), EPS);
    EXPECT_NEAR(2.3984 * std::pow(10, -6), dis.pmf(3, 3), EPS);
    EXPECT_NEAR(0.8 * 2.997999 * std::pow(10, -6), dis.pmf(3, 0), EPS); 

    dis = AdjustedBeta(0.5, 0.5, 0.0);
    EXPECT_NEAR(0.001058687, dis.pmf(10, 1), EPS);
    EXPECT_NEAR(0.0, dis.pmf(0, 0), EPS);
    EXPECT_NEAR(0.020135041, dis.pmf(3, 3), EPS);
    EXPECT_NEAR(0.020135041, dis.pmf(3, 0), EPS);  
    EXPECT_NEAR(0.000640578, dis.pmf(9, 4), EPS); 

    dis = AdjustedBeta(5, 1, 0.0);
    EXPECT_NEAR(5.10100501 * std::pow(10, -7), dis.pmf(10, 1), EPS);
    EXPECT_NEAR(0.004990009, dis.pmf(3, 3), EPS);
    EXPECT_NEAR(0.0001951902, dis.pmf(9, 4), EPS); 
}