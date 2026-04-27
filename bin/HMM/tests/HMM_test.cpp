#include <gtest/gtest.h>

#include <memory>
#include <iostream>
#include <math.h>

#include "../HMM.h"
#include "TROOT.h"
#include <omp.h>

double EPS = std::pow(10, -8);
std::vector<size_t> start = {0};

TEST(HMM, trainPoisson)
{
    Matrix<int> obs {{0, 2, 2}, {1, 10, 15}, {20, 5, 2}, {30, 27, 1}, {2, 0, 5}, {0, 0, 0}, {3, 4, 3}};

    Matrix<std::shared_ptr<DiscreteDistribution>> em (2, 3);
    em(0, 0) = std::make_shared<Poisson>(Poisson(2));
    em(0, 1) = std::make_shared<Poisson>(Poisson(20));
    em(0, 2) = std::make_shared<Poisson>(Poisson(0.5));
    em(1, 0) = std::make_shared<Poisson>(Poisson(12));
    em(1, 1) = std::make_shared<Poisson>(Poisson(1));
    em(1, 2) = std::make_shared<Poisson>(Poisson(2));

    HMM model (2, em);

    HMM::const_matrix_ptr<int> p = std::make_shared<Matrix<int>>(obs);
    double old_likelihood = model.log_likelihood(p);
    HMM::const_matrix_ptr<int> empty = std::make_shared<Matrix<int>>(Matrix<int>());
    model.train(p, empty, start, EPS, 5);
    EXPECT_GE(model.log_likelihood(p), old_likelihood);

    EXPECT_FALSE(model.has_methylation());
}

TEST(HMM, Discrete_one_iteration)
{
    Matrix<int> obs {{0}, {0}, {1}, {2}, {1}, {0}, {0}, {0}, {1}, {1}};

    Matrix<std::shared_ptr<DiscreteDistribution>> em (2, 1);
    em(0, 0) = std::make_shared<Discrete>(Discrete(std::vector<double>{1.0/9, 1.0/3, 5.0/9}));
    em(1, 0) = std::make_shared<Discrete>(Discrete(std::vector<double>{1.0/6, 1.0/3, 1.0/2}));

    HMM model (2, em);
    HMM::const_matrix_ptr<int> p = std::make_shared<Matrix<int>>(obs);
    HMM::const_matrix_ptr<int> empty = std::make_shared<Matrix<int>>(Matrix<int>());
    model.train(p, empty, start, EPS, 1);
    Matrix<double> A = model.get_transition_matrix();
    double precision = std::pow(10, -7);
    EXPECT_NEAR(A(0, 0), 0.4610458, precision);
    EXPECT_NEAR(A(0, 1), 0.5389542, precision);
    EXPECT_NEAR(A(1, 0), 0.4564021, precision);
    EXPECT_NEAR(A(1, 1), 0.5435979, precision);

    EXPECT_NEAR(em(0, 0)->get_parameters()[0], 0.441860465116279, EPS);
    EXPECT_NEAR(em(0, 0)->get_parameters()[1], 0.4418604651162791, EPS);
    EXPECT_NEAR(em(0, 0)->get_parameters()[2], 0.11627906976744189, EPS);

    EXPECT_NEAR(em(1, 0)->get_parameters()[0], 0.5480769230769231, EPS);
    EXPECT_NEAR(em(1, 0)->get_parameters()[1], 0.3653846153846154, EPS);
    EXPECT_NEAR(em(1, 0)->get_parameters()[2], 0.08653846153846155, EPS);
}

TEST(HMM, Discrete_multiple_iterations)
{
    Matrix<int> obs {{0}, {0}, {1}, {2}, {1}, {0}, {0}, {0}, {1}, {1}};

    Matrix<std::shared_ptr<DiscreteDistribution>> em (2, 1);
    em(0, 0) = std::make_shared<Discrete>(Discrete(std::vector<double>{1.0/9, 1.0/3, 5.0/9}));
    em(1, 0) = std::make_shared<Discrete>(Discrete(std::vector<double>{1.0/6, 1.0/3, 1.0/2}));

    HMM model (2, em);
    HMM::const_matrix_ptr<int> p = std::make_shared<Matrix<int>>(obs);
    HMM::const_matrix_ptr<int> empty = std::make_shared<Matrix<int>>(Matrix<int>());
    model.train(p, empty, start, EPS, 50, false);

    Matrix<double> A = model.get_transition_matrix();
    EXPECT_NEAR(A(0, 0), 0.82848337, EPS);
    EXPECT_NEAR(A(0, 1), 0.17151663, EPS);
    EXPECT_NEAR(A(1, 0), 0.53705022, EPS);
    EXPECT_NEAR(A(1, 1), 0.46294978, EPS);

    EXPECT_NEAR(em(0, 0)->get_parameters()[0], 0.24546197010540405, EPS);
    EXPECT_NEAR(em(0, 0)->get_parameters()[1], 0.6036247459992231, EPS);
    EXPECT_NEAR(em(0, 0)->get_parameters()[2], 0.15091328389537284, EPS);

    EXPECT_NEAR(em(1, 0)->get_parameters()[0], 0.9999442393411482, EPS);
    EXPECT_NEAR(em(1, 0)->get_parameters()[1], 5.5760658743875645 * std::pow(10, -5), EPS);
    EXPECT_NEAR(em(1, 0)->get_parameters()[2], 1.0800420126317268 * std::pow(10, -13), EPS);
}

TEST(HMM, viterbi_decoding)
{
    Matrix<int> obs {{2},{2},{1},{0},{1},{3},{2},{0},{0}};
    Matrix<std::shared_ptr<DiscreteDistribution>> em (2, 1);
    em(0, 0) = std::make_shared<Discrete>(Discrete(std::vector<double>{0.2, 0.3, 0.3, 0.2}));
    em(1, 0) = std::make_shared<Discrete>(Discrete(std::vector<double>{0.3, 0.2, 0.2, 0.3}));

    HMM model (2, em);

    Matrix<double> transitions {{0.5, 0.5}, {0.4, 0.6}};
    model.set_transitions(transitions);

    HMM::const_matrix_ptr<int> p = std::make_shared<Matrix<int>>(obs);
    std::pair<double, std::vector<int>> res = model.viterbi_decoding(std::make_pair<size_t, size_t>(0, 9), p);
    std::vector<int> correct {1, 1, 1, 2, 2, 2, 2, 2, 2};

    for (size_t i = 0; i < correct.size(); ++i)
    {
        EXPECT_EQ(res.second[i], correct[i]);
    }
    EXPECT_NEAR(res.first, -16.9734, std::pow(10, -4));
}

TEST(HMM, posterior_decoding)
{
    Matrix<int> obs {{2},{2},{1},{0},{1},{3},{2},{0},{0}};

    Matrix<std::shared_ptr<DiscreteDistribution>> em (2, 1);
    em(0, 0) = std::make_shared<Discrete>(Discrete(std::vector<double>{0.2, 0.3, 0.3, 0.2}));
    em(1, 0) = std::make_shared<Discrete>(Discrete(std::vector<double>{0.3, 0.2, 0.2, 0.3}));

    HMM model (2, em);
    Matrix<double> transitions {{0.7, 0.3}, {0.1, 0.9}};
    model.set_transitions(transitions);

    HMM::const_matrix_ptr<int> p = std::make_shared<Matrix<int>>(obs);
    std::vector<std::pair<int, double>> res = model.posterior_decoding(std::make_pair<size_t, size_t>(0, 9), p);
    std::vector<int> correct {1, 1, 2, 2, 2, 2, 2, 2, 2};

    for (size_t i = 0; i < correct.size(); ++i)
    {
         EXPECT_EQ(res[i].first, correct[i]);
    }
}

TEST(HMM, Gaussian)
{
    Matrix<int> obs {{0}, {2}, {3}, {2}, {10}, {15}, {20}};

    Matrix<std::shared_ptr<DiscreteDistribution>> em (2, 1);
    em(0, 0) = std::make_shared<Gaussian>(Gaussian(1, 10));
    em(1, 0) = std::make_shared<Gaussian>(Gaussian(20, 5));

    HMM model (2, em);
    HMM::const_matrix_ptr<int> p = std::make_shared<Matrix<int>>(obs);
    double old_likelihood = model.log_likelihood(p);
    for (size_t i = 0; i < 5; ++i)
    {
        HMM::const_matrix_ptr<int> empty = std::make_shared<Matrix<int>>(Matrix<int>());
        model.train(p, empty, start, EPS, 1);
        double new_likelihood = model.log_likelihood(p);
        EXPECT_GE(new_likelihood, old_likelihood);
        old_likelihood = new_likelihood;
    }
}

TEST(HMM, Methylation)
{
    ROOT::EnableThreadSafety();
    Matrix<int> obs {{0, 2, 2}, {1, 10, 15}, {20, 5, 2}, {30, 27, 1}, {2, 0, 5}, {0, 0, 0}, {3, 4, 3}};
    Matrix<int> nObs {{10, 0}, {10, 8}, {0, 0}, {8, 0}, {5, 3}, {6, 0}, {4, 0}};

    HMM::const_matrix_ptr<int> o = std::make_shared<Matrix<int>>(obs);
    HMM::const_matrix_ptr<int> nO = std::make_shared<Matrix<int>>(nObs);

    Matrix<std::shared_ptr<DiscreteDistribution>> em (2, 4);
    em(0, 0) = std::make_shared<Poisson>(Poisson(2));
    em(0, 1) = std::make_shared<Poisson>(Poisson(20));
    em(0, 2) = std::make_shared<Poisson>(Poisson(0.5));
    em(0, 3) = std::make_shared<Binomial>(Binomial(0.1));
    em(1, 0) = std::make_shared<Poisson>(Poisson(12));
    em(1, 1) = std::make_shared<Poisson>(Poisson(1));
    em(1, 2) = std::make_shared<Poisson>(Poisson(2));
    em(1, 3) = std::make_shared<Binomial>(Binomial(0.5));

    HMM model (2, em, true);

    EXPECT_TRUE(model.has_methylation());

    double old_likelihood = model.log_likelihood(o, nO);
    model.train(o, nO, start, EPS, 5);
    EXPECT_GE(model.log_likelihood(o, nO), old_likelihood);

    EXPECT_NO_THROW(model.viterbi_decoding(std::make_pair<size_t, size_t>(0,7), o, nO));
    EXPECT_NO_THROW(model.posterior_decoding(std::make_pair<size_t, size_t>(0,7), o, nO));
    EXPECT_THROW(model.viterbi_decoding(std::make_pair<size_t, size_t>(0,7), o), std::invalid_argument);
    EXPECT_THROW(model.posterior_decoding(std::make_pair<size_t, size_t>(0,7), o), std::invalid_argument);

    em(0, 3) = std::make_shared<BetaBinomial>(BetaBinomial(0.7, 2));
    em(1, 3) = std::make_shared<BetaBinomial>(BetaBinomial(0.2, 0.2));

    HMM model2 (2, em, true);
    old_likelihood = model2.log_likelihood(o, nO);
    model2.train(o, nO, start, EPS, 5);
    EXPECT_GE(model2.log_likelihood(o, nO), old_likelihood);
}