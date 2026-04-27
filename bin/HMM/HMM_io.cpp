#include <stdexcept>
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

#include "HMM.h"
#include "adjustableDurationHMM.h"
#include "log_prob.h"

using namespace boost::property_tree;

/**
 * @brief parse parameters for the given distribution
 * 
 * @param states 
 * @param m 
 * @param distribution 
 * @param parameters 
 * @param emission 
 */
void parse_distribution(int states, int m, const std::string& distribution, bool methylation, ptree& parameters, Matrix<std::shared_ptr<DiscreteDistribution>>& emission);

/**
 * @brief read initial state distribution
 * 
 * @param tree 
 * @return std::vector<double> 
 */
std::vector<double> parse_initial_state_distribution(boost::property_tree::ptree& tree);

/**
 * @brief read transition matrix
 * 
 * @return matrix<double> 
 */
Matrix<double> parse_transition_matrix(int, boost::property_tree::ptree&);

/**
 * @brief create property tree for distribution parameters
 * 
 * @return ptree 
 */
ptree get_distribution_parameters(const Matrix<std::shared_ptr<DiscreteDistribution>>&, int, int);

std::ostream& operator<<(std::ostream& out, const HMM& model)
{
    ptree root;

    root.put("states", model.N);

    ptree marker;
    for (auto& name: model.marker)
    {
        ptree marker_entry;
        marker_entry.put("", name);
        marker.push_back(std::make_pair("", marker_entry));
    }
    root.add_child("marker", marker);

    root.put("methylation", model.methylation);

    ptree emission;
    for (size_t i = 0; i < model.emission.nrows(); ++i)
    {
        ptree row;
        for (size_t j = 0; j < model.emission.ncols(); ++j)
        {
            ptree cell;
            cell.put("distribution", model.emission(i, j)->get_name());
            ptree param = get_distribution_parameters(model.emission, i, j);
            cell.add_child("parameters", param);
            row.push_back(std::make_pair("", cell));
        }
        emission.push_back(std::make_pair("", row));
    }
    root.add_child("emission", emission);

    ptree initial;
    for (size_t i = 0; i < model.N; ++i)
    {
        ptree initial_entry;
        double prob = lp::ext_exp(model.logPi[i]);
        if (prob > 1e-300)
        {
            initial_entry.put_value(prob);
        }
        else
        {
            initial_entry.put_value(0.0);
        }
        initial.push_back(std::make_pair("", initial_entry));
    }
    root.add_child("initial_state_distribution", initial);

    ptree transition;
    for (size_t i = 0; i < model.N; i++)
    {
        ptree row;
        for (size_t j = 0; j < model.N; j++)
        {
            ptree cell;
            double prob = lp::ext_exp(model.logA(i, j));
            if (prob > 1e-300)
            {
                cell.put_value(prob);
            }
            else
            {
                cell.put_value(0.0);
            }
            row.push_back(std::make_pair("", cell));
        }
        transition.push_back(std::make_pair("", row));
    }
    root.add_child("transition_matrix", transition);

    write_json(out, root);
    return out;
}

std::istream& operator>>(std::istream& in, HMM& model)
{
    ptree root;
    try
    {
        read_json(in, root);
    }
    catch(const std::exception& e)
    {
        std::string ex = e.what();
        throw std::ios_base::failure("Wrong json format: " + ex);
    }

    if (root.find("states") == root.not_found())
        throw std::ios_base::failure("Missing <states> parameter.");
    model.N = root.get<int>("states");
    
    if (root.find("marker") == root.not_found())
        throw std::ios_base::failure("Missing <marker> parameter.");
    for (ptree::value_type& p : root.get_child("marker"))
    {
        model.marker.push_back(p.second.get_value<std::string>());
    }
    model.m = model.marker.size();
    
    model.methylation = false;
    if (root.find("methylation") != root.not_found())
    {
        model.methylation = root.get<bool>("methylation");
    }

    if (root.find("emission") == root.not_found())
        throw std::ios_base::failure("Missing <emission> parameter.");

    size_t dim = model.methylation ? model.m + 1 : model.m;
    model.emission = Matrix<std::shared_ptr<DiscreteDistribution>> (model.N, dim);
    int s = 0;
    for (ptree::value_type& row : root.get_child("emission"))
    {
        if (s >= model.N)
            throw std::ios_base::failure("Too many states.");
        
        int m = 0;
        for (ptree::value_type& cell : row.second)
        {
            if (m >= dim)
                throw std::ios_base::failure("Distributions per state and number of markers must be the same (plus one if methylation is true).");
            
            std::string distribution = cell.second.get<std::string>("distribution");
            bool meth = false;
            if (model.methylation && m == dim-1)
                meth = true;
            parse_distribution(s, m, distribution, meth, cell.second.get_child("parameters"), model.emission);
            ++m;
        }
        if (m != dim)
            throw std::ios_base::failure("Distributions per state and number of markers must be the same.");
        ++s;
    }
    if (s != model.N)
        throw std::ios_base::failure("Not enough states in emission matrix.");


    if (root.find("initial_state_distribution") != root.not_found())
    {
        std::vector<double> pi = parse_initial_state_distribution(root);
        model.set_initial(pi);
    }
    else
    {
        model.init_initial();
    }

    if (root.find("transition_matrix") != root.not_found())
    {
        Matrix<double> transition = parse_transition_matrix(model.N, root);
        model.set_transitions(transition);
    }
    else
    {
        model.init_transitions();
    }

    return in;
}

void parse_distribution(int s, int m, const std::string& distribution, bool methylation, ptree& parameters, Matrix<std::shared_ptr<DiscreteDistribution>>& emission)
{
    if (!methylation)
    {
        if (distribution == "PO")
        {
            if (parameters.find("lambda") == parameters.not_found())
                throw std::ios_base::failure("Poisson distribution requires parameter <lambda>.");
            emission(s, m) = std::make_shared<Poisson>(Poisson(parameters.get<double>("lambda")));
        }
        else if (distribution == "NBI")
        {
            if (parameters.find("p") == parameters.not_found() || parameters.find("r") == parameters.not_found())
                throw std::ios_base::failure("Negative binomial distribution requires parameters <p> and <r>.");
            emission(s, m) = std::make_shared<NegativeBinomial>(NegativeBinomial(parameters.get<double>("p"), parameters.get<double>("r")));
        }
        else if (distribution == "BI")
        {
            if (parameters.find("p") == parameters.not_found() || parameters.find("n") == parameters.not_found())
                throw std::ios_base::failure("Binomial distribution requires parameters <p> and <n>.");
            emission(s, m) = std::make_shared<Binomial>(Binomial(parameters.get<double>("p"), parameters.get<int>("n")));
        }
        else if (distribution == "BB")
        {
            if (parameters.find("alpha") == parameters.not_found() || parameters.find("beta") == parameters.not_found() || parameters.find("n") == parameters.not_found())
                throw std::ios_base::failure("Beta binomial distribution requires parameters <alpha>, <beta> and <n>.");
            emission(s, m) = std::make_shared<BetaBinomial>(BetaBinomial(parameters.get<double>("alpha"), parameters.get<double>("beta"), parameters.get<int>("n")));
        }
        else if (distribution == "BNB")
        {
            if (parameters.find("alpha") == parameters.not_found() || parameters.find("beta") == parameters.not_found() || parameters.find("r") == parameters.not_found())
                throw std::ios_base::failure("Beta negative binomial distribution requires parameters <alpha>, <beta> and <r>.");
            
            emission(s, m) = std::make_shared<BetaNegativeBinomial>(BetaNegativeBinomial(parameters.get<double>("alpha"), parameters.get<double>("beta"), parameters.get<double>("r")));
        }
        else if (distribution == "ZAP")
        {
            if (parameters.find("lambda") == parameters.not_found() || parameters.find("pi") == parameters.not_found())
                throw std::ios_base::failure("Zero Adjusted Poisson requires parameters <lambda> and <pi>.");
            emission(s, m) = std::make_shared<ZeroAdjustedPoisson>(ZeroAdjustedPoisson(parameters.get<double>("lambda"), parameters.get<double>("pi")));
        }
        else if (distribution == "ZANBI")
        {
            if (parameters.find("p") == parameters.not_found() || parameters.find("r") == parameters.not_found() || parameters.find("pi") == parameters.not_found())
                throw std::ios_base::failure("Zero adjusted negative binomial distribution requires parameters <p> and <r> and <pi>.");
            
            emission(s, m) = std::make_shared<ZeroAdjustedNegativeBinomial>(ZeroAdjustedNegativeBinomial(parameters.get<double>("p"), parameters.get<double>("r"), parameters.get<double>("pi")));
        }
        else if (distribution == "ZABNB")
        {
            if (parameters.find("alpha") == parameters.not_found() || parameters.find("beta") == parameters.not_found() || parameters.find("r") == parameters.not_found() || parameters.find("pi") == parameters.not_found())
                throw std::ios_base::failure("Zero adjusted beta negative binomial distribution requires parameters <alpha>, <beta>, <r> and <pi>.");
            
            emission(s, m) = std::make_shared<ZeroAdjustedBetaNegativeBinomial>(ZeroAdjustedBetaNegativeBinomial(parameters.get<double>("alpha"), parameters.get<double>("beta"), parameters.get<double>("r"), parameters.get<double>("pi")));
        }
        else if (distribution == "GA")
        {        
            if (parameters.find("mean") == parameters.not_found() || parameters.find("std") == parameters.not_found())
                throw std::ios_base::failure("Gaussian distribution requires parameters <mean> and <std>.");
            emission(s, m) = std::make_shared<Gaussian>(Gaussian(parameters.get<double>("mean"), parameters.get<double>("std")));
        }
        else if (distribution == "SI")
        {        
            if (parameters.find("mu") == parameters.not_found() || parameters.find("sigma") == parameters.not_found() || parameters.find("v") == parameters.not_found())
                throw std::ios_base::failure("Sichel distribution requires parameters <mu> and <sigma> and <v>.");
            emission(s, m) = std::make_shared<Sichel>(Sichel(parameters.get<double>("mu"), parameters.get<double>("sigma"), parameters.get<double>("v")));
        }
        else if (distribution == "ZASI")
        {        
            if (parameters.find("mu") == parameters.not_found() || parameters.find("sigma") == parameters.not_found() || parameters.find("v") == parameters.not_found() || parameters.find("pi") == parameters.not_found())
                throw std::ios_base::failure("Zero adjusted Sichel distribution requires parameters <mu>, <sigma>, <v> and <pi>.");
            emission(s, m) = std::make_shared<ZeroAdjustedSichel>(ZeroAdjustedSichel(parameters.get<double>("mu"), parameters.get<double>("sigma"), parameters.get<double>("v"), parameters.get<double>("pi")));
        }
        else if (distribution == "B")
        {
            if (parameters.find("p") == parameters.not_found())
                throw std::ios_base::failure("Bernoulli distribution requires parameter <p>.");
            emission(s, m) = std::make_shared<Bernoulli>(Bernoulli(parameters.get<double>("p")));
        }
        else 
        {
            throw std::ios_base::failure("Unknown distribution.");
        }
    }
    else
    {
        if (distribution == "BI")
        {
            if (parameters.find("p") == parameters.not_found())
                throw std::ios_base::failure("Binomial distribution requires parameter <p>.");
            emission(s, m) = std::make_shared<Binomial>(Binomial(parameters.get<double>("p")));
        }
        else if (distribution == "BB")
        {
            if (parameters.find("alpha") == parameters.not_found() || parameters.find("beta") == parameters.not_found())
                throw std::ios_base::failure("Beta binomial distribution requires parameters <alpha> and <beta>.");
            emission(s, m) = std::make_shared<BetaBinomial>(BetaBinomial(parameters.get<double>("alpha"), parameters.get<double>("beta")));
        }
        else if (distribution == "AB")
        {        
            if (parameters.find("alpha") == parameters.not_found() || parameters.find("beta") == parameters.not_found() || parameters.find("pi") == parameters.not_found())
                throw std::ios_base::failure("Adjusted beta distribution requires parameters <alpha>, <beta> and <pi>.");
            emission(s, m) = std::make_shared<AdjustedBeta>(AdjustedBeta(parameters.get<double>("alpha"), parameters.get<double>("beta"), parameters.get<double>("pi")));
        }
        else 
        {
            throw std::ios_base::failure("Unsupported distribution for DNA methylation.");
        }
    }
}

ptree get_distribution_parameters(const Matrix<std::shared_ptr<DiscreteDistribution>>& emission, int i, int j)
{
    ptree parameters;
    size_t N = emission.nrows();

    if (emission(i, j)->get_name() == "PO")
    {
        parameters.put("lambda", emission(i, j)->get_parameters()[0]);
    }
    else if (emission(i, j)->get_name() == "NBI")
    {
        parameters.put("p", emission(i, j)->get_parameters()[0]);
        parameters.put("r", emission(i, j)->get_parameters()[1]);
    }
    else if (emission(i, j)->get_name() == "BI")
    {
        parameters.put("p", emission(i, j)->get_parameters()[0]);
        int n = emission(i, j)->get_parameters()[1];
        if (n > 0)
        {
            parameters.put("n", n);
        }
    }
    else if (emission(i, j)->get_name() == "BB")
    {
        parameters.put("alpha", emission(i, j)->get_parameters()[0]);
        parameters.put("beta", emission(i, j)->get_parameters()[1]);

        int n = emission(i, j)->get_parameters()[2];
        if (n > 0)
        {
            parameters.put("n", n);
        }
    }
    else if (emission(i, j)->get_name() == "BNB")
    {
        parameters.put("alpha", emission(i, j)->get_parameters()[0]);
        parameters.put("beta", emission(i, j)->get_parameters()[1]);
        parameters.put("r", emission(i, j)->get_parameters()[2]);
    }
    else if (emission(i, j)->get_name() == "ZAP")
    {
        parameters.put("lambda", emission(i, j)->get_parameters()[0]);
        parameters.put("pi", emission(i, j)->get_parameters()[1]);
    }
    else if (emission(i, j)->get_name() == "ZANBI")
    {
        parameters.put("p", emission(i, j)->get_parameters()[0]);
        parameters.put("r", emission(i, j)->get_parameters()[1]);
        parameters.put("pi", emission(i, j)->get_parameters()[2]);
    }
    else if (emission(i, j)->get_name() == "ZABNB")
    {
        parameters.put("alpha", emission(i, j)->get_parameters()[0]);
        parameters.put("beta", emission(i, j)->get_parameters()[1]);
        parameters.put("r", emission(i, j)->get_parameters()[2]);
        parameters.put("pi", emission(i, j)->get_parameters()[3]);
    }
    else if (emission(i, j)->get_name() == "GA")
    {        
        parameters.put("mean", emission(i, j)->get_parameters()[0]);
        parameters.put("std", emission(i, j)->get_parameters()[1]);
    }
    else if (emission(i, j)->get_name() == "SI")
    {        
        parameters.put("mu", emission(i, j)->get_parameters()[0]);
        parameters.put("sigma", emission(i, j)->get_parameters()[1]);
        parameters.put("v", emission(i, j)->get_parameters()[2]);
    }
    else if (emission(i, j)->get_name() == "ZASI")
    {        
        parameters.put("mu", emission(i, j)->get_parameters()[0]);
        parameters.put("sigma", emission(i, j)->get_parameters()[1]);
        parameters.put("v", emission(i, j)->get_parameters()[2]);
        parameters.put("pi", emission(i, j)->get_parameters()[3]);
    }
    else if (emission(i, j)->get_name() == "AB")
    {        
        parameters.put("alpha", emission(i, j)->get_parameters()[0]);
        parameters.put("beta", emission(i, j)->get_parameters()[1]);
        parameters.put("pi", emission(i, j)->get_parameters()[2]);
    }
    else if (emission(i, j)->get_name() == "B")
    {        
        parameters.put("p", emission(i, j)->get_parameters()[0]);
    }

    return parameters;
}

Matrix<double> parse_transition_matrix(int states, ptree& tree)
{
    Matrix<double> transitionMatrix = Matrix<double>(states, states);
    int x = 0;
    for (ptree::value_type& row : tree.get_child("transition_matrix"))
    {
        if (x >= states)
                throw std::ios_base::failure("Transition matrix must have shape states x states");

        int y = 0;
        for (ptree::value_type& cell : row.second)
        {
            if (y >= states)
                throw std::ios_base::failure("Transition matrix must have shape states x states");
            transitionMatrix(x, y) = cell.second.get_value<double>();
            ++y;
        }
        if (y != states)
            throw std::ios_base::failure("Transition matrix must have shape states x states");

        ++x;
    }
    if (x != states)
        throw std::ios_base::failure("Transition matrix must have shape states x states");
    return transitionMatrix;
}

std::vector<double> parse_initial_state_distribution(ptree& tree)
{
    std::vector<double> pi;
    for (ptree::value_type& p : tree.get_child("initial_state_distribution"))
    {
        pi.push_back(p.second.get_value<double>());
    }
    return pi;
}

// Adjustable duration HMM
std::ostream& operator<<(std::ostream& out, const AdjustableDurationHMM& model)
{
    ptree root;

    root.put("states", model.stateIndices.size());

    ptree topology;
    for (size_t s = 0; s < model.stateIndices.size(); ++s)
    {
        ptree topology_entry;
        for (auto subState : model.stateIndices[s])
        {
            ptree subHMM_entry;
            subHMM_entry.put("", subState);
            topology_entry.push_back(std::make_pair("", subHMM_entry));
        }
        topology.add_child(std::to_string(s+1), topology_entry);
    }
    root.add_child("topology", topology);

    ptree marker;
    for (auto& name: model.marker)
    {
        ptree marker_entry;
        marker_entry.put("", name);
        marker.push_back(std::make_pair("", marker_entry));
    }
    root.add_child("marker", marker);

    root.put("methylation", model.methylation);

    ptree emission;
    for (size_t i = 0; i < model.emission.nrows(); ++i)
    {
        ptree row;
        for (size_t j = 0; j < model.emission.ncols(); ++j)
        {
            ptree cell;
            cell.put("distribution", model.emission(i, j)->get_name());
            ptree param = get_distribution_parameters(model.emission, i, j);
            cell.add_child("parameters", param);
            row.push_back(std::make_pair("", cell));
        }
        emission.push_back(std::make_pair("", row));
    }
    root.add_child("emission", emission);

    ptree initial;
    for (size_t i = 0; i < model.N; ++i)
    {
        ptree initial_entry;
        double prob = lp::ext_exp(model.logPi[i]);
        if (prob > 1e-300)
        {
            initial_entry.put_value(prob);
        }
        else
        {
            initial_entry.put_value(0.0);
        }
        initial.push_back(std::make_pair("", initial_entry));
    }
    root.add_child("initial_state_distribution", initial);

    ptree transition;
    for (size_t i = 0; i < model.N; i++)
    {
        ptree row;
        for (size_t j = 0; j < model.N; j++)
        {
            ptree cell;
            double prob = lp::ext_exp(model.logA(i, j));
            if (prob > 1e-300)
            {
                cell.put_value(prob);
            }
            else
            {
                cell.put_value(0.0);
            }
            row.push_back(std::make_pair("", cell));
        }
        transition.push_back(std::make_pair("", row));
    }
    root.add_child("transition_matrix", transition);

    write_json(out, root);
    return out;
}

std::istream& operator>>(std::istream& in, AdjustableDurationHMM& model)
{
    ptree root;
    try
    {
        read_json(in, root);
    }
    catch(const std::exception& e)
    {
        std::string ex = e.what();
        throw std::ios_base::failure("Wrong json format: " + ex);
    }

    if (root.find("marker") == root.not_found())
        throw std::ios_base::failure("Missing <marker> parameter.");
    for (ptree::value_type& p : root.get_child("marker"))
    {
        model.marker.push_back(p.second.get_value<std::string>());
    }
    model.m = model.marker.size();

    if (root.find("states") == root.not_found())
        throw std::ios_base::failure("Missing <states> parameter.");
    size_t states = root.get<int>("states");

    if (root.find("topology") == root.not_found())
        throw std::ios_base::failure("Missing <topology> parameter.");

    model.N = 0;
    for (ptree::value_type& state : root.get_child("topology"))
    { 
        std::vector<size_t> subHMM;
        for (ptree::value_type& s : state.second)
        {
            subHMM.push_back(s.second.get_value<size_t>());
            ++model.N;
        }
        model.stateIndices.push_back(subHMM);
    }
    if (states != model.stateIndices.size())
        throw std::ios_base::failure("Missing state in topology description.");
    model.stateAssignment = std::vector<size_t>(model.N);
    for (size_t i = 0; i < states; ++i)
    {
        for (auto s : model.stateIndices[i])
        {
            model.stateAssignment[s] = i;
        }
    }

    model.methylation = false;
    if (root.find("methylation") != root.not_found())
    {
        model.methylation = root.get<bool>("methylation");
    }

    if (root.find("emission") == root.not_found())
        throw std::ios_base::failure("Missing <emission> parameter.");

    size_t dim = model.methylation ? model.m + 1 : model.m;

    model.emission = Matrix<std::shared_ptr<DiscreteDistribution>> (model.stateIndices.size(), dim);
    int s = 0;
    for (ptree::value_type& row : root.get_child("emission"))
    {
        if (s >= model.stateIndices.size())
            throw std::ios_base::failure("Too many states.");
        
        int m = 0;
        for (ptree::value_type& cell : row.second)
        {
            if (m >= dim)
                throw std::ios_base::failure("Distributions per state and number of markers must be the same (plus one if methylation is true).");
            
            std::string distribution = cell.second.get<std::string>("distribution");
            bool meth = false;
            if (model.methylation && m == dim-1)
                meth = true;
            parse_distribution(s, m, distribution, meth, cell.second.get_child("parameters"), model.emission);
            ++m;
        }
        if (m != dim)
            throw std::ios_base::failure("Distributions per state and number of markers must be the same.");
        ++s;
    }
    if (s != model.stateIndices.size())
        throw std::ios_base::failure("Not enough states in emission matrix.");


    if (root.find("initial_state_distribution") != root.not_found())
    {
        std::vector<double> pi = parse_initial_state_distribution(root);
        model.set_initial(pi);
    }
    else
    {
        model.init_initial();
    }

    if (root.find("transition_matrix") != root.not_found())
    {
        Matrix<double> transition = parse_transition_matrix(model.N, root);
        model.set_transitions(transition);
    }
    else
    {
        model.init_transitions();
    }

    return in;
}