#ifndef GUARD_funs_h
#define GUARD_funs_h

#include <vector>
#include <map>
#include <set>
#include <string>
#include <functional>

#include "tree.h"


void tree_traversal(suff_stat &ss, tree &t, data_info &di);

// void fit_single_tree(double* ftemp, tree &t, data_info &di);
void fit_ensemble(std::vector<double> &fit, std::vector<tree> &t_vec, data_info &di);

void compute_suff_stat_grow(suff_stat &orig_suff_stat, suff_stat &new_suff_stat,
                            int &nx_nid, rule_t &rule, tree &t, data_info &di);

void compute_suff_stat_prune(suff_stat &orig_suff_stat, suff_stat &new_suff_stat,
                             int &nl_nid, int &nr_nid, int &np_nid, tree &t, data_info &di);

void compute_p_theta_ind(int &j, arma::mat &P, arma::vec &Theta,
                         std::map<int,int> &leaf_map, suff_stat &ss,
                         double &sigma, data_info &di, tree_prior_info &tree_pi);

void compute_p_theta_cs(int &j, arma::mat &P, arma::vec &Theta,
                        std::map<int,int> &leaf_map, suff_stat &ss,
                        double &rho, double &sigma, data_info &di, tree_prior_info &tree_pi);

double compute_lil(arma::mat &P, arma::vec &Theta, tree_prior_info &tree_pi);

void draw_mu(tree &t, arma::mat &P, arma::vec &Theta, std::map<int,int> &leaf_map, RNG &gen);

void draw_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);

std::string write_tree(tree &t, tree_prior_info &tree_pi, set_str_conversion &set_str);
void read_tree(tree &t, std::string &tree_string, set_str_conversion &set_str);

void build_symmetric_edge_map(edge_map &emap, std::vector<edge> &edges, std::set<int> &vertices);
std::vector<edge> get_induced_edges(std::vector<edge> &edges, std::set<int> &vertex_subset);
void dfs(int v, std::map<int, bool> &visited, std::vector<int> &comp, edge_map &emap);
void find_components(std::vector<std::vector<int> > &components, std::vector<edge> &edges, std::set<int> &vertices);
void get_unique_edges(std::vector<edge> &edges);
void boruvka(std::vector<edge> &mst_edges, std::vector<edge> &edges, std::set<int> &vertices);
void wilson(std::vector<edge> &mst_edges, std::vector<edge> &edges, std::set<int> &vertices, RNG &gen);
void graph_partition(std::set<int> &avail_levels, std::set<int> &l_vals, std::set<int> &r_vals,
                     std::vector<edge> &orig_edges, int &K, int &cut_type, RNG &gen);

void update_theta_u(std::vector<double> &theta, double &u, std::vector<int> &var_count, int &R,
                    double &a_u, double &b_u, RNG &gen);

void update_theta_rc(double& theta_rc, int &rc_var_count, int &rc_rule_count,
                     double &a_rc, double &b_rc, int &R_cont, RNG &gen);

void update_rho(double &rho, double &sigma, data_info &di, RNG &gen);

void update_sigma_ind(double &sigma, double &nu, double &lambda, data_info &di, RNG &gen);
void update_sigma_cs(double &sigma, double &rho, double &nu, double &lambda, data_info &di, RNG &gen);


//  helpers for global–local (regularized horseshoe) shrinkage
// ------------------------------------------------------------------
//
// Keep tree updates Gaussian by using the effective RHS prior variance
//   s2_j = (tau^2 * lambda_j^2 * c2 / (c2 + tau^2 * lambda_j^2)) / M
// shared by all leaves that belong to predictor j.
//
// Overloads of P,Theta builders that take a scalar prior variance per leaf.
//    Use these from update_trees when GL shrinkage is enabled.
//
void compute_p_theta_ind_gl(int &j, arma::mat &P, arma::vec &Theta,
                            std::map<int,int> &leaf_map, suff_stat &ss,
                            double &sigma, data_info &di, tree_prior_info &tree_pi,
                            double prior_var_per_leaf);

void compute_p_theta_cs_gl(int &j, arma::mat &P, arma::vec &Theta,
                           std::map<int,int> &leaf_map, suff_stat &ss,
                           double &rho, double &sigma, data_info &di, tree_prior_info &tree_pi,
                           double prior_var_per_leaf);

// Small utilities for RHS variance/shrink factor.
inline double rhs_shrink_factor(double tau, double lambda, double c2) {
  // sqrt( c2 / (c2 + tau^2 * lambda^2) )
  double tl = tau * lambda;
  return std::sqrt( c2 / (c2 + tl*tl) );
}

inline double rhs_prior_var(double tau, double lambda, double c2, int M) {
  // (tau^2 * lambda^2 * c2 / (c2 + tau^2 * lambda^2)) / M
  double tl = tau * lambda;
  double num = tl*tl * c2;
  double den = c2 + tl*tl;
  return (num / den) / static_cast<double>(M);
}

// Log–posterior kernels on the log–scale for slice/MH updates.
//    • log λ_j | μ_j, τ, c2  (uses sum of squares across ALL leaves of predictor j)
//    • log τ   | {μ_j}, {λ_j}, c2  (global)
//    • log c2  | {μ_j}, {λ_j}, τ    (slab IG prior: a = df/2, b = df*scale^2/2)
//
// All return the log unnormalized posterior density at the provided log-parameter.
//
double logpost_loglambda_rhs(double log_lambda_j,
                             double tau, double c2,
                             double sumsq_mu_j, int n_mu_j,
                             int M);

double logpost_logtau_rhs(double log_tau,
                          const std::vector<double>& sumsq_by_j,
                          const std::vector<int>& n_by_j,
                          const std::vector<double>& lambda_by_j,
                          double c2, int M, double tau0);

double logpost_logc2_rhs(double log_c2,
                         const std::vector<double>& sumsq_by_j,
                         const std::vector<int>& n_by_j,
                         const std::vector<double>& lambda_by_j,
                         double tau, int M,
                         double a_slab, double b_slab);

// 1-D slice sampler on log–scale for the above kernels.
//    Inputs:
//      • logf:   functor returning log-density at x (already on log-scale of parameter)
//      • x0:     current value (log-parameter)
//      • w:      step-out width
//      • m:      max step-out steps
//      • lower/upper: hard bounds on x (use -inf/+inf if none)
//
double slice1d_on_log(const std::function<double(double)>& logf,
                      double x0, RNG& gen,
                      double w, int m,
                      double lower, double upper);

#endif /* GUARD_funs_h */
