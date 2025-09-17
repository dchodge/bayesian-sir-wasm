#include "mcmc.hpp"
#include <algorithm>
#include <cstring>

MCMCSampler::MCMCSampler(unsigned seed) 
    : rng(seed), normal_dist(0.0, 1.0), uniform_dist(0.0, 1.0) {}

SIRData MCMCSampler::simulateSIR(double beta, double gamma, double initial_infected_pct, 
                                int time_period, int population) {
    SIRData data;
    double dt = 1.0; // time step in days
    
    // Initial conditions
    double S = population - (initial_infected_pct * population / 100.0);
    double I = initial_infected_pct * population / 100.0;
    double R = 0.0;
    
    for (int t = 0; t <= time_period; ++t) {
        data.time.push_back(t);
        data.S.push_back((S / population) * 100.0);
        data.I.push_back((I / population) * 100.0);
        data.R.push_back((R / population) * 100.0);
        
        if (t < time_period) {
            // SIR differential equations (Euler's method)
            double dS = -beta * S * I / population;
            double dI = beta * S * I / population - gamma * I;
            double dR = gamma * I;
            
            S += dS * dt;
            I += dI * dt;
            R += dR * dt;
        }
    }
    
    return data;
}

double MCMCSampler::logLikelihood(const MCMCParams& params, const ObservedData& observed, 
                                 const SIRData& predicted) {
    if (params.beta <= 0 || params.gamma <= 0 || params.sigma <= 0) {
        return -std::numeric_limits<double>::infinity();
    }
    
    double log_lik = 0.0;
    
    // Match observed timepoints with predicted data
    for (size_t i = 0; i < observed.time.size(); ++i) {
        double t = observed.time[i];
        double obs_I = observed.I_observed[i];
        
        // Find closest predicted value
        auto it = std::min_element(predicted.time.begin(), predicted.time.end(),
                                  [t](double a, double b) { 
                                      return std::abs(a - t) < std::abs(b - t); 
                                  });
        
        if (it != predicted.time.end()) {
            size_t pred_idx = std::distance(predicted.time.begin(), it);
            if (pred_idx < predicted.I.size()) {
                double pred_I = predicted.I[pred_idx];
                
                // Normal likelihood
                double residual = (obs_I - pred_I) / params.sigma;
                log_lik += -0.5 * std::log(2 * M_PI * params.sigma * params.sigma) - 
                          0.5 * residual * residual;
            }
        }
    }
    
    return log_lik;
}

double MCMCSampler::logPrior(const MCMCParams& params, const MCMCPriors& priors) {
    if (params.beta <= 0 || params.gamma <= 0 || params.sigma <= 0) {
        return -std::numeric_limits<double>::infinity();
    }
    
    double log_prior = 0.0;
    
    // Beta ~ Normal(mean, sd)
    double beta_residual = (params.beta - priors.beta_mean) / priors.beta_sd;
    log_prior += -0.5 * std::log(2 * M_PI * priors.beta_sd * priors.beta_sd) - 
                0.5 * beta_residual * beta_residual;
    
    // Gamma ~ Normal(mean, sd)
    double gamma_residual = (params.gamma - priors.gamma_mean) / priors.gamma_sd;
    log_prior += -0.5 * std::log(2 * M_PI * priors.gamma_sd * priors.gamma_sd) - 
                0.5 * gamma_residual * gamma_residual;
    
    // Sigma ~ InverseGamma(a, b)
    log_prior += priors.sigma_a * std::log(priors.sigma_b) - 
                std::lgamma(priors.sigma_a) - 
                (priors.sigma_a + 1) * std::log(params.sigma) - 
                priors.sigma_b / params.sigma;
    
    return log_prior;
}

MCMCParams MCMCSampler::proposeParameters(const MCMCParams& current, const ProposalSDs& proposal_sds) {
    MCMCParams proposed = current;
    
    proposed.beta += normal_dist(rng) * proposal_sds.beta_sd;
    proposed.gamma += normal_dist(rng) * proposal_sds.gamma_sd;
    proposed.sigma += normal_dist(rng) * proposal_sds.sigma_sd;
    
    return proposed;
}

MCMCSampler::MCMCResult MCMCSampler::mcmcStep(const MCMCParams& current_params, 
                                             double current_log_post,
                                             const ObservedData& observed, 
                                             const MCMCPriors& priors,
                                             const ProposalSDs& proposal_sds, 
                                             double initial_infected_pct, 
                                             int time_period) {
    
    // Propose new parameters
    MCMCParams proposed_params = proposeParameters(current_params, proposal_sds);
    
    // Run SIR simulation with proposed parameters
    SIRData proposed_data = simulateSIR(proposed_params.beta, proposed_params.gamma, 
                                       initial_infected_pct, time_period);
    
    // Calculate proposed log-posterior
    double proposed_log_lik = logLikelihood(proposed_params, observed, proposed_data);
    double proposed_log_prior = logPrior(proposed_params, priors);
    double proposed_log_post = proposed_log_lik + proposed_log_prior;
    
    // Accept/reject
    double log_accept_ratio = proposed_log_post - current_log_post;
    bool accept = std::log(uniform_dist(rng)) < log_accept_ratio;
    
    MCMCResult result;
    if (accept) {
        result.params = proposed_params;
        result.log_posterior = proposed_log_post;
        result.accepted = true;
        result.predicted_data = proposed_data;
    } else {
        result.params = current_params;
        result.log_posterior = current_log_post;
        result.accepted = false;
        // Don't need predicted data for rejected step
    }
    
    return result;
}

// ============= MATRIX OPERATIONS =============

Matrix3x3 MCMCSampler::choleskyDecomposition(const Matrix3x3& matrix) {
    Matrix3x3 L{};
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            if (j == i) {
                for (int k = 0; k < j; k++) {
                    sum += L[j][k] * L[j][k];
                }
                double diagonal = matrix[j][j] - sum;
                if (diagonal <= 0) {
                    throw std::runtime_error("Non-positive diagonal in Cholesky decomposition");
                }
                L[j][j] = std::sqrt(diagonal);
            } else {
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }
                if (std::abs(L[j][j]) < 1e-12) {
                    throw std::runtime_error("Near-zero diagonal in Cholesky decomposition");
                }
                L[i][j] = (matrix[i][j] - sum) / L[j][j];
                
                if (!std::isfinite(L[i][j])) {
                    throw std::runtime_error("Invalid Cholesky element");
                }
            }
        }
    }
    return L;
}

Vector3 MCMCSampler::multiplyMatrixVector(const Matrix3x3& matrix, const Vector3& vector) {
    Vector3 result{};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    return result;
}

Matrix3x3 MCMCSampler::scaleMatrix(const Matrix3x3& matrix, double scale) {
    Matrix3x3 result{};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i][j] = matrix[i][j] * scale;
        }
    }
    return result;
}

Matrix3x3 MCMCSampler::addRegularization(const Matrix3x3& matrix, double epsilon) {
    Matrix3x3 result = matrix;
    for (int i = 0; i < 3; i++) {
        result[i][i] += epsilon;
    }
    return result;
}

Vector3 MCMCSampler::calculateMean(const std::vector<Vector3>& samples) {
    Vector3 mean{};
    if (samples.empty()) return mean;
    
    for (const auto& sample : samples) {
        for (int i = 0; i < 3; i++) {
            mean[i] += sample[i];
        }
    }
    
    double n = static_cast<double>(samples.size());
    for (int i = 0; i < 3; i++) {
        mean[i] /= n;
    }
    
    return mean;
}

Matrix3x3 MCMCSampler::calculateCovariance(const std::vector<Vector3>& samples, const Vector3& mean) {
    Matrix3x3 cov{};
    if (samples.size() <= 1) return cov;
    
    for (const auto& sample : samples) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                cov[i][j] += (sample[i] - mean[i]) * (sample[j] - mean[j]);
            }
        }
    }
    
    double n = static_cast<double>(samples.size() - 1);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cov[i][j] /= n;
        }
    }
    
    return cov;
}

// ============= ADAPTIVE MCMC METHODS =============

MCMCParams MCMCSampler::proposeParametersAdaptive(const MCMCParams& current, AdaptiveData& adaptive_data) {
    Vector3 current_vec = paramsToVector(current);
    MCMCParams proposed;
    
    try {
        if (adaptive_data.adaptive && adaptive_data.has_valid_cholesky) {
            // Generate standard normal samples
            Vector3 z = {normal_dist(rng), normal_dist(rng), normal_dist(rng)};
            
            // Apply Robbins-Monro scaling
            Matrix3x3 scaled_cov = scaleMatrix(adaptive_data.proposal_cov, adaptive_data.scaling_factor);
            
            // Decompose and sample
            Matrix3x3 L = choleskyDecomposition(scaled_cov);
            Vector3 delta = multiplyMatrixVector(L, z);
            
            // Create proposed parameters
            Vector3 proposed_vec = {
                current_vec[0] + delta[0],
                current_vec[1] + delta[1], 
                current_vec[2] + delta[2]
            };
            
            proposed = vectorToParams(proposed_vec);
            
        } else {
            // Fall back to scaled fixed proposals
            double base_scale = 0.01;
            double effective_scale = base_scale * adaptive_data.scaling_factor;
            
            proposed.beta = current.beta + normal_dist(rng) * effective_scale;
            proposed.gamma = current.gamma + normal_dist(rng) * effective_scale;
            proposed.sigma = current.sigma + normal_dist(rng) * effective_scale;
        }
        
        // Apply bounds
        proposed.beta = std::max(0.001, std::min(5.0, proposed.beta));
        proposed.gamma = std::max(0.001, std::min(2.0, proposed.gamma));
        proposed.sigma = std::max(0.001, std::min(1.0, proposed.sigma));
        
    } catch (const std::exception& e) {
        // Fallback to non-adaptive if there's an error
        adaptive_data.adaptive = false;
        adaptive_data.has_valid_cholesky = false;
        
        double base_scale = 0.01;
        proposed.beta = current.beta + normal_dist(rng) * base_scale;
        proposed.gamma = current.gamma + normal_dist(rng) * base_scale;
        proposed.sigma = current.sigma + normal_dist(rng) * base_scale;
        
        // Apply bounds
        proposed.beta = std::max(0.001, std::min(5.0, proposed.beta));
        proposed.gamma = std::max(0.001, std::min(2.0, proposed.gamma));
        proposed.sigma = std::max(0.001, std::min(1.0, proposed.sigma));
    }
    
    return proposed;
}

void MCMCSampler::updateEmpiricalStats(AdaptiveData& adaptive_data) {
    if (adaptive_data.samples.size() < AdaptiveData::min_samples) {
        return;
    }
    
    try {
        // Calculate empirical mean and covariance
        adaptive_data.empirical_mean = calculateMean(adaptive_data.samples);
        adaptive_data.empirical_cov = calculateCovariance(adaptive_data.samples, adaptive_data.empirical_mean);
        
        // Add regularization
        adaptive_data.empirical_cov = addRegularization(adaptive_data.empirical_cov, AdaptiveData::epsilon);
        
        // Scale for proposal
        adaptive_data.proposal_cov = scaleMatrix(adaptive_data.empirical_cov, AdaptiveData::adaptation_scale);
        
        // Test Cholesky decomposition
        choleskyDecomposition(adaptive_data.proposal_cov);
        adaptive_data.has_valid_cholesky = true;
        
    } catch (const std::exception& e) {
        adaptive_data.has_valid_cholesky = false;
        adaptive_data.adaptive = false;
    }
}

void MCMCSampler::updateRobbinsMonroScaling(AdaptiveData& adaptive_data, bool accepted) {
    adaptive_data.acceptance_history.push_back(accepted ? 1.0 : 0.0);
    
    // Keep only recent history
    if (adaptive_data.acceptance_history.size() > 100) {
        adaptive_data.acceptance_history.erase(adaptive_data.acceptance_history.begin());
    }
    
    if (adaptive_data.acceptance_history.size() >= 10) {
        // Calculate recent acceptance rate
        double recent_acceptance = 0.0;
        int recent_count = std::min(50, static_cast<int>(adaptive_data.acceptance_history.size()));
        for (int i = adaptive_data.acceptance_history.size() - recent_count; 
             i < adaptive_data.acceptance_history.size(); i++) {
            recent_acceptance += adaptive_data.acceptance_history[i];
        }
        recent_acceptance /= recent_count;
        
        // Robbins-Monro update
        double step_size = 1.0 / std::pow(adaptive_data.acceptance_history.size(), 0.6);
        double error = recent_acceptance - AdaptiveData::target_acceptance;
        adaptive_data.log_scaling += step_size * error;
        
        // Apply bounds to log scaling
        adaptive_data.log_scaling = std::max(-10.0, std::min(2.0, adaptive_data.log_scaling));
        adaptive_data.scaling_factor = std::exp(adaptive_data.log_scaling);
        
        adaptive_data.robbins_monro_active = true;
    }
}

bool MCMCSampler::enableAdaptiveSampling(AdaptiveData& adaptive_data) {
    if (adaptive_data.samples.size() >= AdaptiveData::min_samples) {
        updateEmpiricalStats(adaptive_data);
        adaptive_data.adaptive = adaptive_data.has_valid_cholesky;
        return adaptive_data.adaptive;
    }
    return false;
}

MCMCSampler::AdaptiveMCMCResult MCMCSampler::mcmcStepAdaptive(
    const MCMCParams& current_params, double current_log_post,
    const ObservedData& observed, const MCMCPriors& priors,
    AdaptiveData& adaptive_data, double initial_infected_pct, int time_period) {
    
    // Add current sample to adaptive data
    Vector3 current_vec = paramsToVector(current_params);
    adaptive_data.samples.push_back(current_vec);
    
    // Keep sample history manageable
    if (adaptive_data.samples.size() > AdaptiveData::max_samples) {
        adaptive_data.samples.erase(adaptive_data.samples.begin());
    }
    
    // Propose new parameters
    MCMCParams proposed_params = proposeParametersAdaptive(current_params, adaptive_data);
    
    // Run SIR simulation with proposed parameters
    SIRData proposed_data = simulateSIR(proposed_params.beta, proposed_params.gamma, 
                                       initial_infected_pct, time_period);
    
    // Calculate proposed log-posterior
    double proposed_log_lik = logLikelihood(proposed_params, observed, proposed_data);
    double proposed_log_prior = logPrior(proposed_params, priors);
    double proposed_log_post = proposed_log_lik + proposed_log_prior;
    
    // Accept/reject
    double log_accept_ratio = proposed_log_post - current_log_post;
    bool accept = std::log(uniform_dist(rng)) < log_accept_ratio;
    
    // Update Robbins-Monro scaling
    updateRobbinsMonroScaling(adaptive_data, accept);
    
    // Calculate acceptance rate
    double acceptance_rate = 0.0;
    if (!adaptive_data.acceptance_history.empty()) {
        for (double acc : adaptive_data.acceptance_history) {
            acceptance_rate += acc;
        }
        acceptance_rate /= adaptive_data.acceptance_history.size();
    }
    
    // Prepare result
    AdaptiveMCMCResult result;
    if (accept) {
        result.params = proposed_params;
        result.log_posterior = proposed_log_post;
        result.accepted = true;
        result.predicted_data = proposed_data;
    } else {
        result.params = current_params;
        result.log_posterior = current_log_post;
        result.accepted = false;
    }
    
    result.used_adaptive = adaptive_data.adaptive;
    result.acceptance_rate = acceptance_rate;
    result.scaling_factor = adaptive_data.scaling_factor;
    
    return result;
}

// ============= UTILITY FUNCTIONS =============

Vector3 MCMCSampler::paramsToVector(const MCMCParams& params) {
    return {params.beta, params.gamma, params.sigma};
}

MCMCParams MCMCSampler::vectorToParams(const Vector3& vec) {
    MCMCParams params;
    params.beta = vec[0];
    params.gamma = vec[1];
    params.sigma = vec[2];
    return params;
}

// C interface implementations
extern "C" {

// ============= SAMPLER MANAGEMENT =============
    
MCMCSampler* create_sampler(unsigned seed) {
    return new MCMCSampler(seed);
}

void destroy_sampler(MCMCSampler* sampler) {
    delete sampler;
}

// ============= ADAPTIVE DATA MANAGEMENT =============

AdaptiveData* create_adaptive_data() {
    return new AdaptiveData();
}

void destroy_adaptive_data(AdaptiveData* adaptive_data) {
    delete adaptive_data;
}

void reset_adaptive_data(AdaptiveData* adaptive_data) {
    if (!adaptive_data) return;
    
    adaptive_data->samples.clear();
    adaptive_data->empirical_mean = {};
    adaptive_data->empirical_cov = {};
    adaptive_data->proposal_cov = {};
    adaptive_data->cholesky_L = {};
    adaptive_data->acceptance_history.clear();
    adaptive_data->scaling_factor = 1.0;
    adaptive_data->log_scaling = 0.0;
    adaptive_data->adaptive = false;
    adaptive_data->robbins_monro_active = false;
    adaptive_data->has_valid_cholesky = false;
}

// ============= SIR SIMULATION =============

int run_sir_simulation(MCMCSampler* sampler,
                      double beta, double gamma, double initial_infected_pct,
                      int time_period, int population,
                      double* out_times, double* out_S, double* out_I, double* out_R,
                      int* out_length) {
    
    if (!sampler || !out_times || !out_S || !out_I || !out_R || !out_length) {
        return 0;
    }
    
    try {
        SIRData data = sampler->simulateSIR(beta, gamma, initial_infected_pct, time_period, population);
        
        *out_length = std::min(static_cast<int>(data.time.size()), time_period + 1);
        
        for (int i = 0; i < *out_length; ++i) {
            out_times[i] = data.time[i];
            out_S[i] = data.S[i];
            out_I[i] = data.I[i];
            out_R[i] = data.R[i];
        }
        
        return 1; // success
    } catch (...) {
        return 0; // error
    }
}

// ============= BASIC MCMC STEP (LEGACY) =============

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
              double* pred_times, double* pred_infected, int* pred_length) {
    
    try {
        // Setup current parameters
        MCMCParams current_params = {current_beta, current_gamma, current_sigma};
        
        // Setup observed data
        ObservedData observed;
        for (int i = 0; i < n_obs; ++i) {
            observed.time.push_back(obs_times[i]);
            observed.I_observed.push_back(obs_infected[i]);
        }
        
        // Setup priors
        MCMCPriors priors = {beta_prior_mean, beta_prior_sd, 
                            gamma_prior_mean, gamma_prior_sd,
                            sigma_prior_a, sigma_prior_b};
        
        // Setup proposal standard deviations
        ProposalSDs proposal_sds = {beta_prop_sd, gamma_prop_sd, sigma_prop_sd};
        
        // Perform MCMC step
        auto result = sampler->mcmcStep(current_params, current_log_post, observed, 
                                       priors, proposal_sds, initial_infected_pct, time_period);
        
        // Copy results
        *new_beta = result.params.beta;
        *new_gamma = result.params.gamma;
        *new_sigma = result.params.sigma;
        *new_log_post = result.log_posterior;
        *accepted = result.accepted ? 1 : 0;
        
        // Copy predicted data if accepted
        if (result.accepted) {
            *pred_length = std::min((int)result.predicted_data.time.size(), time_period + 1);
            for (int i = 0; i < *pred_length; ++i) {
                pred_times[i] = result.predicted_data.time[i];
                pred_infected[i] = result.predicted_data.I[i];
            }
        } else {
            *pred_length = 0;
        }
        
        return 1; // success
        
    } catch (...) {
        return 0; // error
    }
}

// ============= ADAPTIVE MCMC STEP =============

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
                      double* pred_times, double* pred_infected, int* pred_length) {
    
    if (!sampler || !adaptive_data) return 0;
    
    try {
        // Setup current parameters
        MCMCParams current_params = {current_beta, current_gamma, current_sigma};
        
        // Setup observed data
        ObservedData observed;
        for (int i = 0; i < n_obs; ++i) {
            observed.time.push_back(obs_times[i]);
            observed.I_observed.push_back(obs_infected[i]);
        }
        
        // Setup priors
        MCMCPriors priors = {beta_prior_mean, beta_prior_sd, 
                            gamma_prior_mean, gamma_prior_sd,
                            sigma_prior_a, sigma_prior_b};
        
        // Perform adaptive MCMC step
        auto result = sampler->mcmcStepAdaptive(current_params, current_log_post, observed, 
                                               priors, *adaptive_data, initial_infected_pct, time_period);
        
        // Copy results
        *new_beta = result.params.beta;
        *new_gamma = result.params.gamma;
        *new_sigma = result.params.sigma;
        *new_log_post = result.log_posterior;
        *accepted = result.accepted ? 1 : 0;
        *used_adaptive = result.used_adaptive ? 1 : 0;
        *acceptance_rate = result.acceptance_rate;
        *scaling_factor = result.scaling_factor;
        
        // Copy predicted data if accepted
        if (result.accepted) {
            *pred_length = std::min((int)result.predicted_data.time.size(), time_period + 1);
            for (int i = 0; i < *pred_length; ++i) {
                pred_times[i] = result.predicted_data.time[i];
                pred_infected[i] = result.predicted_data.I[i];
            }
        } else {
            *pred_length = 0;
        }
        
        return 1; // success
        
    } catch (...) {
        return 0; // error
    }
}

// ============= ADAPTIVE MCMC MANAGEMENT =============

int enable_adaptive_sampling(AdaptiveData* adaptive_data) {
    if (!adaptive_data) return 0;
    
    if (adaptive_data->samples.size() >= AdaptiveData::min_samples) {
        adaptive_data->adaptive = true;
        return 1;
    }
    return 0;
}

int update_empirical_stats(MCMCSampler* sampler, AdaptiveData* adaptive_data) {
    if (!sampler || !adaptive_data) return 0;
    
    try {
        sampler->updateEmpiricalStats(*adaptive_data);
        return adaptive_data->has_valid_cholesky ? 1 : 0;
    } catch (...) {
        return 0;
    }
}

void add_sample_to_adaptive(AdaptiveData* adaptive_data, 
                           double beta, double gamma, double sigma) {
    if (!adaptive_data) return;
    
    Vector3 sample = {beta, gamma, sigma};
    adaptive_data->samples.push_back(sample);
    
    // Keep sample history manageable
    if (adaptive_data->samples.size() > AdaptiveData::max_samples) {
        adaptive_data->samples.erase(adaptive_data->samples.begin());
    }
}

// ============= UTILITY FUNCTIONS =============

void get_adaptive_info(AdaptiveData* adaptive_data,
                      int* is_adaptive, int* is_robbins_monro_active,
                      double* scaling_factor, double* acceptance_rate,
                      int* num_samples) {
    
    if (!adaptive_data) {
        if (is_adaptive) *is_adaptive = 0;
        if (is_robbins_monro_active) *is_robbins_monro_active = 0;
        if (scaling_factor) *scaling_factor = 1.0;
        if (acceptance_rate) *acceptance_rate = 0.0;
        if (num_samples) *num_samples = 0;
        return;
    }
    
    if (is_adaptive) *is_adaptive = adaptive_data->adaptive ? 1 : 0;
    if (is_robbins_monro_active) *is_robbins_monro_active = adaptive_data->robbins_monro_active ? 1 : 0;
    if (scaling_factor) *scaling_factor = adaptive_data->scaling_factor;
    if (num_samples) *num_samples = static_cast<int>(adaptive_data->samples.size());
    
    if (acceptance_rate) {
        if (adaptive_data->acceptance_history.empty()) {
            *acceptance_rate = 0.0;
        } else {
            double sum = 0.0;
            for (double acc : adaptive_data->acceptance_history) {
                sum += acc;
            }
            *acceptance_rate = sum / adaptive_data->acceptance_history.size();
        }
    }
}

} // extern "C"
