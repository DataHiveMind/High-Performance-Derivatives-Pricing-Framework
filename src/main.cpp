#include "option_pricer.hpp"
#include "black_scholes.hpp"
#include "binomial_tree.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>

using namespace hpdf;

namespace {

void print_result(const OptionSpec& spec, const PriceResult& res, const std::string& model, bool csv = false) {
    if (csv) {
        std::cout << spec.S << "," << spec.K << "," << spec.r << "," << spec.q << "," << spec.sigma << "," << spec.T
                  << "," << (spec.type == OptionType::Call ? "Call" : "Put")
                  << "," << model << "," << res.price << "," << res.greeks.delta << "," << res.greeks.gamma
                  << "," << res.greeks.vega << "," << res.greeks.theta << "," << res.greeks.rho << "\n";
    } else {
        std::cout << "Model: " << model << "\n"
                  << "Type: " << (spec.type == OptionType::Call ? "Call" : "Put") << "\n"
                  << "S=" << spec.S << " K=" << spec.K << " r=" << spec.r << " q=" << spec.q
                  << " sigma=" << spec.sigma << " T=" << spec.T << "\n"
                  << "Price: " << res.price << "\n"
                  << "Delta: " << res.greeks.delta << "  Gamma: " << res.greeks.gamma
                  << "  Vega: " << res.greeks.vega << "  Theta: " << res.greeks.theta
                  << "  Rho: " << res.greeks.rho << "\n";
    }
}

OptionSpec parse_option_from_args(int argc, char** argv) {
    OptionSpec spec;
    spec.type = OptionType::Call;
    spec.S = 100.0;
    spec.K = 100.0;
    spec.r = constants::DEFAULT_R;
    spec.q = constants::DEFAULT_Q;
    spec.sigma = 0.2;
    spec.T = 1.0;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--put") spec.type = OptionType::Put;
        else if (arg == "--call") spec.type = OptionType::Call;
        else if (arg == "--S" && i + 1 < argc) spec.S = std::stod(argv[++i]);
        else if (arg == "--K" && i + 1 < argc) spec.K = std::stod(argv[++i]);
        else if (arg == "--r" && i + 1 < argc) spec.r = std::stod(argv[++i]);
        else if (arg == "--q" && i + 1 < argc) spec.q = std::stod(argv[++i]);
        else if (arg == "--sigma" && i + 1 < argc) spec.sigma = std::stod(argv[++i]);
        else if (arg == "--T" && i + 1 < argc) spec.T = std::stod(argv[++i]);
        // ignore dividends for CLI
    }
    return spec;
}

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [--call|--put] [--S val] [--K val] [--r val] [--q val] [--sigma val] [--T val] "
              << "[--steps N] [--model bs|binomial] [--csv] [--batch file.csv] [--bench]\n";
}

int main(int argc, char** argv) {
    OptionSpec spec = parse_option_from_args(argc, argv);
    int steps = 200;
    std::string model = "bs";
    bool csv = false, batch = false, bench = false;
    std::string batch_file = "data/sample_options.csv";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--steps" && i + 1 < argc) steps = std::stoi(argv[++i]);
        else if (arg == "--model" && i + 1 < argc) model = argv[++i];
        else if (arg == "--csv") csv = true;
        else if (arg == "--batch" && i + 1 < argc) { batch = true; batch_file = argv[++i]; }
        else if (arg == "--bench") bench = true;
        else if (arg == "--help" || arg == "-h") { print_usage(argv[0]); return 0; }
    }

    if (batch) {
        auto rows = utils::load_csv(batch_file);
        if (!csv) {
            std::cout << "Batch mode: " << rows.size() << " options from " << batch_file << "\n";
            std::cout << "S,K,r,q,sigma,T,Type,Model,Price,Delta,Gamma,Vega,Theta,Rho\n";
        }
        std::vector<OptionSpec> specs;
        for (const auto& row : rows) {
            if (row.size() < 7) continue;
            OptionSpec s;
            s.S = utils::parse_param(row[0], 100.0);
            s.K = utils::parse_param(row[1], 100.0);
            s.r = utils::parse_param(row[2], constants::DEFAULT_R);
            s.q = utils::parse_param(row[3], constants::DEFAULT_Q);
            s.sigma = utils::parse_param(row[4], 0.2);
            s.T = utils::parse_param(row[5], 1.0);
            s.type = (row[6] == "Put" || row[6] == "put") ? OptionType::Put : OptionType::Call;
            specs.push_back(s);
        }
        std::vector<PriceResult> results(specs.size());
        auto t0 = std::chrono::high_resolution_clock::now();
        if (model == "bs") {
            for (size_t i = 0; i < specs.size(); ++i)
                results[i] = black_scholes::BlackScholesPricer().price(specs[i]);
        } else {
            for (size_t i = 0; i < specs.size(); ++i)
                results[i] = price_binomial(specs[i], steps);
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < specs.size(); ++i)
            print_result(specs[i], results[i], model, csv);
        if (bench) {
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            std::cout << "Batch latency: " << ms << " ms (" << specs.size() << " options)\n";
        }
        return 0;
    }

    // Single price mode
    auto t0 = std::chrono::high_resolution_clock::now();
    PriceResult res;
    if (model == "bs") {
        res = black_scholes::BlackScholesPricer().price(spec);
    } else {
        res = price_binomial(spec, steps);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    print_result(spec, res, model, csv);
    if (bench) {
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << "Latency: " << ms << " ms\n";
    }
    return 0;
}