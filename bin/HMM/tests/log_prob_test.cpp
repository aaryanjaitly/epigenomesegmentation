#include <gtest/gtest.h>
#include <cmath>
#include <math.h>
#include <vector>
#include "../log_prob.h"

double EPS = std::pow(10, -8);

TEST(LogProb, ext_exp)
{
    EXPECT_NEAR(0.0, lp::ext_exp(NAN), EPS);
    EXPECT_NEAR(0.00551656442, lp::ext_exp(-5.2), EPS);
    EXPECT_NEAR(15.958634, lp::ext_exp(2.77), EPS);
}

TEST(LogProb, ext_log)
{
    EXPECT_NEAR(lp::ext_exp(lp::ext_log(0)), 0.0, EPS);
    EXPECT_NEAR(-6.074846156, lp::ext_log(0.0023), EPS);
    EXPECT_NEAR(1.5475625, lp::ext_log(4.7), EPS);
}

TEST(LogProb, log_add)
{
    EXPECT_NEAR(-3.45, lp::log_add(NAN, -3.45), EPS);
    EXPECT_NEAR(-23.1, lp::log_add(-23.1, NAN), EPS);
    EXPECT_NEAR(-7.8, lp::log_add(-43.2, -7.8), EPS);
    EXPECT_NEAR(-2.3600466668, lp::log_add(-2.4, -5.6), EPS);
}

TEST(LogProb, log_mul)
{
    EXPECT_TRUE(isnan(lp::log_mul(NAN, -3.45)));
    EXPECT_TRUE(isnan(lp::log_mul(-23.1, NAN)));
    EXPECT_NEAR(-51.0, lp::log_mul(-43.2, -7.8), EPS);
}

TEST(LogProb, log_sum)
{
    std::vector<double> v1 {-43.2, -7.8};
    EXPECT_NEAR(-7.8, lp::log_sum(v1.begin(), v1.end()), EPS);
    std::vector<double> v2 {-2.4, -5.6, NAN, -23.2};
    EXPECT_NEAR(-2.36004666, lp::log_sum(v2.begin(), v2.end()), EPS);
    std::vector<double> v3 {-2.4, -5.6, NAN, -23.2, -45.2};
    EXPECT_NEAR(-2.36004666, lp::log_sum(v3.begin(), v3.end()), EPS);
    std::vector<double> v4 {-2.4, -12.4, -5.6, NAN, -23.2, -45.2, -3.9};
    EXPECT_NEAR(-2.165767942, lp::log_sum(v4.begin(), v4.end()), EPS);
}

TEST(LogProb, log_prob)
{
    std::vector<double> v1 {-2, -12, -3, -1, -2};
    EXPECT_NEAR(-20.0, lp::log_prod(v1.begin(), v1.end()), EPS);
}

TEST(LogProb, log_greater)
{
    EXPECT_TRUE(lp::log_greater(-10.2, -53.7));
    EXPECT_TRUE(lp::log_greater(-10.2, NAN));
    EXPECT_FALSE(lp::log_greater(NAN, -42.0));
    EXPECT_FALSE(lp::log_greater(NAN, NAN));
}