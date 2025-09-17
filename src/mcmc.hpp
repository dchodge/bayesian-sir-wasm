#pragma once
#include <vector>
#include <random>
#include <cmath>
#include <array>
#include <memory>

struct SIRData {
    std::vector<double> time;
    std::vector<double> S;
    std::vector<double> I;
    std::vector<double> R;
};

struct MCMCParams {
    double beta;
    double gamma; 
    double sigma;
};

struct MCMCPriors {
    double beta_mean, beta_sd;
    double gamma_mean, gamma_sd;
    double sigma_a, sigma_b; // inverse gamma parameters
};

struct ProposalSDs {
    double beta_sd;
    double gamma_sd;
    double sigma_sd;
};

struct ObservedData {
    std::vector<double> time;
    std::vector<double> I_observed;
};

// 3x3 Matrix for covariance operations
using Matrix3x3 = std::array<std::array<double, 3>, 3>;
using Vector3 = std::array<double, 3>;

// Adaptive MCMC structures
struct AdaptiveData {
    // Empirical statistics
    std::vector<Vector3> samples;
    Vector3 empirical_mean{};
    Matrix3x3 empirical_cov{};
    Matrix3x3 proposal_cov{};
    Matrix3x3 cholesky_L{};
    
    // Robbins-Monro scaling
    std::vector<double> acceptance_history;
    double scaling_factor = 1.0;
    double log_scaling = 0.0;
    
    // Control flags
    bool adaptive = false;
    bool robbins_monro_active = false;
    bool has_valid_cholesky = false;
    
    // Parameters
    static constexpr double target_acceptance = 0.234;
    static constexpr double adaptation_scale = 2.38 * 2.38 / 3.0; // Optimal scaling
    static constexpr double epsilon = 1e-6;
    static constexpr int min_samples = 10;
    static constexpr int max_samples = 1000;
};

// Enhanced proposal structure
struct EnhancedProposalSDs {
    double beta_sd;
    double gamma_sd;
    double sigma_sd;
    
    // Enhanced features
    bool use_adaptive = false;
    AdaptiveData* adaptive_data = nullptr;
};

class MCMCSampler {
private:
    std::mt19937 rng;
    std::normal_distribution<double> normal_dist;
    std::uniform_real_distribution<double> uniform_dist;
    
    // Matrix operations
    Matrix3x3 choleskyDecomposition(const Matrix3x3& matrix);
    Vector3 multiplyMatrixVector(const Matrix3x3& matrix, const Vector3& vector);
    Matrix3x3 scaleMatrix(const Matrix3x3& matrix, double scale);
    Matrix3x3 addRegularization(const Matrix3x3& matrix, double epsilon);
    
    // Statistical operations
    Vector3 calculateMean(const std::vector<Vector3>& samples);
    Matrix3x3 calculateCovariance(const std::vector<Vector3>& samples, const Vector3& mean);
    
public:
    MCMCSampler(unsigned seed = 12345);
    
    // SIR model simulation
    SIRData simulateSIR(double beta, double gamma, double initial_infected_pct, 
                        int time_period, int population = 100000);
    
    // Likelihood calculations
    double logLikelihood(const MCMCParams& params, const ObservedData& observed, 
                        const SIRData& predicted);
    
    // Prior calculations  
    double logPrior(const MCMCParams& params, const MCMCPriors& priors);
    
    // Parameter proposals (legacy)
    MCMCParams proposeParameters(const MCMCParams& current, const ProposalSDs& proposal_sds);
    
    // Enhanced parameter proposals with adaptive capability
    MCMCParams proposeParametersAdaptive(const MCMCParams& current, AdaptiveData& adaptive_data);
    
    // Adaptive MCMC management
    void updateEmpiricalStats(AdaptiveData& adaptive_data);
    void updateRobbinsMonroScaling(AdaptiveData& adaptive_data, bool accepted);
    bool enableAdaptiveSampling(AdaptiveData& adaptive_data);
    
    // MCMC step results
    struct MCMCResult {
        MCMCParams params;
        double log_posterior;
        bool accepted;
        SIRData predicted_data;
    };
    
    // Enhanced MCMC result with adaptive info
    struct AdaptiveMCMCResult {
        MCMCParams params;
        double log_posterior;
        bool accepted;
        SIRData predicted_data;
        bool used_adaptive;
        double acceptance_rate;
        double scaling_factor;
    };
    
    // MCMC step methods
    MCMCResult mcmcStep(const MCMCParams& current_params, double current_log_post,
                       const ObservedData& observed, const MCMCPriors& priors,
                       const ProposalSDs& proposal_sds, double initial_infected_pct, 
                       int time_period);
                       
    AdaptiveMCMCResult mcmcStepAdaptive(const MCMCParams& current_params, double current_log_post,
                                       const ObservedData& observed, const MCMCPriors& priors,
                                       AdaptiveData& adaptive_data, double initial_infected_pct,
                                       int time_period);
    
    // Full MCMC chain
    struct MCMCChain {
        std::vector<MCMCParams> params;
        std::vector<double> log_likelihood;
        std::vector<bool> accepted;
        int total_steps;
        int accepted_steps;
    };
    
    MCMCChain runMCMC(const MCMCParams& initial_params, const ObservedData& observed,
                     const MCMCPriors& priors, const ProposalSDs& proposal_sds,
                     int n_steps, double initial_infected_pct, int time_period);
                     
    // Utility functions
    Vector3 paramsToVector(const MCMCParams& params);
    MCMCParams vectorToParams(const Vector3& vec);
};

// C interface for Emscripten
extern "C" {
    // Sampler management
    MCMCSampler* create_sampler(unsigned seed);
    void destroy_sampler(MCMCSampler* sampler);
    
    // Adaptive data management
    AdaptiveData* create_adaptive_data();
    void destroy_adaptive_data(AdaptiveData* adaptive_data);
    void reset_adaptive_data(AdaptiveData* adaptive_data);
    
    // SIR simulation
    int run_sir_simulation(MCMCSampler* sampler,
                          double beta, double gamma, double initial_infected_pct,
                          int time_period, int population,
                          double* out_times, double* out_S, double* out_I, double* out_R,
                          int* out_length);
    
    // Basic MCMC step (legacy)
    int mcmc_step(MCMCSampler* sampler, 
                  double current_beta, double current_gamma, double current_sigma,
                  double current_log_post,
                  double* obs_times, double* obs_infected, int n_obs,
                  double beta_prior_mean, double beta_prior_sd,
                  double gamma_prior_mean, double gamma_prior_sd,
                  double sigma_prior_a, double sigma_prior_b,
                  double beta_prop_sd, double gamma_prop_sd, double sigma_prop_sd,
                  double initial_infected_pct, int time_period,
                  // Output parameters
                  double* new_beta, double* new_gamma, double* new_sigma,
                  double* new_log_post, int* accepted,
                  double* pred_times, double* pred_infected, int* pred_length);
    
    // Adaptive MCMC step
    int mcmc_step_adaptive(MCMCSampler* sampler, AdaptiveData* adaptive_data,
                          double current_beta, double current_gamma, double current_sigma,
                          double current_log_post,
                          double* obs_times, double* obs_infected, int n_obs,
                          double beta_prior_mean, double beta_prior_sd,
                          double gamma_prior_mean, double gamma_prior_sd,
                          double sigma_prior_a, double sigma_prior_b,
                          double initial_infected_pct, int time_period,
                          // Output parameters
                          double* new_beta, double* new_gamma, double* new_sigma,
                          double* new_log_post, int* accepted, int* used_adaptive,
                          double* acceptance_rate, double* scaling_factor,
                          double* pred_times, double* pred_infected, int* pred_length);
    
    // Adaptive MCMC management
    int enable_adaptive_sampling(AdaptiveData* adaptive_data);
    int update_empirical_stats(MCMCSampler* sampler, AdaptiveData* adaptive_data);
    void add_sample_to_adaptive(AdaptiveData* adaptive_data, 
                               double beta, double gamma, double sigma);
    
    // Utility functions
    void get_adaptive_info(AdaptiveData* adaptive_data,
                          int* is_adaptive, int* is_robbins_monro_active,
                          double* scaling_factor, double* acceptance_rate,
                          int* num_samples);
}
