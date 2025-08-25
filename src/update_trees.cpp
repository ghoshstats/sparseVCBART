#include "update_trees.h"
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

// ---------------------------------------------------------------------
// For numerical robustness
// ---------------------------------------------------------------------
static constexpr double EPS_DIAG = 1e-10;  // tiny ridge for P when assembling
static constexpr int    MIN_LEAF = 3;      // avoid tiny leaves on grow (this is optional)

// ---------------------------------------------------------------------
// Robust log-marginal likelihood and posterior draw of mu via Cholesky
// K = P + (1/tau^2) I + tiny ridge
// ---------------------------------------------------------------------
static inline double compute_lil_safe(const arma::mat& P,
                                      const arma::vec& Theta,
                                      const tree_prior_info& pi)
{
  arma::mat K = P;
  K.diag() += 1.0 / (pi.tau * pi.tau);  // prior precision
  K.diag() += 1e-8;                     // numerical ridge

  arma::mat L;
  bool ok = arma::chol(L, K, "lower");
  if(!ok){
    K.diag() += 1e-6;
    ok = arma::chol(L, K, "lower");
    if(!ok) Rcpp::stop("[compute_lil_safe]: chol failed");
  }

  const double logdet = 2.0 * arma::sum(arma::log(L.diag()));
  arma::vec y = arma::solve(arma::trimatl(L), Theta);
  const double quad = arma::dot(y, y);
  return -0.5 * (logdet - quad);
}

static inline void draw_mu_safe(tree& t,
                                const arma::mat& P, const arma::vec& Theta,
                                const std::map<int,int>& leaf_map,
                                RNG& gen, const tree_prior_info& pi)
{
  arma::mat K = P;
  K.diag() += 1.0 / (pi.tau * pi.tau);
  K.diag() += 1e-8;

  arma::mat L;
  bool ok = arma::chol(L, K, "lower");
  if(!ok){
    K.diag() += 1e-6;
    ok = arma::chol(L, K, "lower");
    if(!ok) Rcpp::stop("[draw_mu_safe]: chol failed");
  }

  // posterior mean m = K^{-1} Theta
  arma::vec y = arma::solve(arma::trimatl(L), Theta);
  arma::vec m = arma::solve(arma::trimatu(L.t()), y);

  // sample eps = L^{-T} z with z ~ N(0, I)
  arma::vec z = arma::randn<arma::vec>(Theta.n_rows);
  arma::vec eps = arma::solve(arma::trimatu(L.t()), z);

  int k = 0;
  for(auto it = leaf_map.begin(); it != leaf_map.end(); ++it, ++k){
    t.get_ptr(it->first)->set_mu( m(k) + eps(k) );
  }
}


// Build P (diagonal) and Theta for IND likelihood for ensemble j
// P_ll = sum_i (x_ij^2) / sigma^2 ;  Theta_l = sum_i (x_ij * r_i) / sigma^2
static inline void
compute_p_theta_ind_local(int j,
                          arma::mat& P, arma::vec& Theta, std::map<int,int>& leaf_map,
                          const suff_stat& ss, double sigma,
                          const data_info& di)
{
  const double invs2 = 1.0 / (sigma * sigma);
  const int L = static_cast<int>(ss.size());
  P.zeros(L, L);
  Theta.zeros(L);
  leaf_map.clear();

  int k = 0;
  for (const auto& kv : ss) {
    const int leaf_id = kv.first;
    const std::vector<int>& idx = kv.second;

    double S1 = 0.0, S2 = 0.0;
    for (int i : idx) {
      const double xji = di.x[j + i * di.p]; // di.x is t(X): row=j, col=i
      const double ri  = di.rp[i];
      S1 += xji * ri;
      S2 += xji * xji;
    }
    const double Qxx = invs2 * std::max(S2, EPS_DIAG);
    P(k,k)   = Qxx;
    Theta(k) = invs2 * S1;

    leaf_map[leaf_id] = k;
    ++k;
  }
}

// Build P (diagonal) and Theta for CS likelihood for ensemble j
// Work subject-by-subject to apply Σ^{-1} under CS structure.
static inline void
compute_p_theta_cs_local(int j,
                         arma::mat& P, arma::vec& Theta, std::map<int,int>& leaf_map,
                         const suff_stat& ss, double rho, double sigma,
                         const data_info& di)
{
  const int L = static_cast<int>(ss.size());
  P.zeros(L, L);
  Theta.zeros(L);
  leaf_map.clear();

  const double inv_scale = 1.0 / (sigma * sigma * (1.0 - rho));

  int k = 0;
  for (const auto& kv : ss) {
    const int leaf_id = kv.first;
    const std::vector<int>& idx = kv.second;

    // per-subject aggregates for this leaf
    std::vector<double> a(di.n, 0.0);  // sum xji * r_i
    std::vector<double> b(di.n, 0.0);  // sum xji
    std::vector<double> c(di.n, 0.0);  // sum xji^2
    std::vector<double> rs(di.n, 0.0); // sum r_i
    std::vector<int>    ns(di.n, 0);

    for (int i : idx) {
      const int s   = di.subj_id[i];          // subject index (0-based)
      const double xji = di.x[j + i * di.p];
      const double ri  = di.rp[i];

      a[s]  += xji * ri;
      b[s]  += xji;
      c[s]  += xji * xji;
      rs[s] += ri;
      ns[s] += 1;
    }

    double Qxr_raw = 0.0, Qxx_raw = 0.0;
    for (int s = 0; s < di.n; ++s) if (ns[s] > 0) {
      const double alpha = rho / (1.0 - rho + ns[s] * rho);
      const double adj = c[s] - alpha * b[s] * b[s];
      Qxx_raw += (adj > 0.0 ? adj : 0.0);
      Qxr_raw += (a[s]  - alpha * b[s] * rs[s]);
    }

    const double Qxx = inv_scale * std::max(Qxx_raw, EPS_DIAG);
    const double Qxr = inv_scale * Qxr_raw;

    P(k,k)   = Qxx;
    Theta(k) = Qxr;

    leaf_map[leaf_id] = k;
    ++k;
  }
}

// ---------------------------------------------------------------------
//  IND moves
// ---------------------------------------------------------------------

void grow_tree_ind(tree &t, int &accept, int &j, double &sigma, suff_stat &ss_train, suff_stat &ss_test,
                   data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  std::vector<int> bn_nid_vec; // bottom-node ids
  bn_nid_vec.reserve(ss_train.size());
  for(const auto& kv : ss_train) bn_nid_vec.push_back(kv.first);

  int ni = static_cast<int>(std::floor(gen.uniform() * bn_nid_vec.size()));
  int nx_nid = bn_nid_vec[ni];
  int nxl_nid = 2*nx_nid;
  int nxr_nid = 2*nx_nid+1;
  tree::tree_p nx = t.get_ptr(nx_nid);
  tree::tree_cp nxp = nx->get_p();

  double q_grow_old = tree_pi.prob_b;
  double q_prune_new = 1.0 - tree_pi.prob_b;
  int nleaf_old = t.get_nbots();
  int nnog_old  = t.get_nnogs();
  int nnog_new  = nnog_old;

  if(nxp == 0) q_grow_old = 1.0;
  else if(!nxp->is_nog()) nnog_new = 1 + nnog_old;

  double log_trans_ratio = (std::log(q_prune_new) - std::log((double)nnog_new))
                         - (std::log(q_grow_old)  - std::log((double)nleaf_old));

  double p_grow_nx  = tree_pi.alpha/std::pow(1.0 + (double) nx->get_depth(), tree_pi.beta);
  double p_grow_nxl = tree_pi.alpha/std::pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  double p_grow_nxr = tree_pi.alpha/std::pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  double log_prior_ratio = std::log(p_grow_nx) + std::log(1.0 - p_grow_nxl) + std::log(1.0 - p_grow_nxr)
                         - std::log(1.0 - p_grow_nx);

  // draw a rule and build proposed sufficient-stat map
  rule_t rule;
  draw_rule(rule, t, nx_nid, di_train, tree_pi, gen);

  suff_stat prop_ss_train;
  compute_suff_stat_grow(ss_train, prop_ss_train, nx_nid, rule, t, di_train);

  suff_stat prop_ss_test;
  if(di_test.N > 0) compute_suff_stat_grow(ss_test, prop_ss_test, nx_nid, rule, t, di_test);

  // Likelihood parts via x-weighted stats (for current/original tree)
  std::map<int,int> orig_leaf_map;
  arma::mat orig_P (ss_train.size(), ss_train.size(), arma::fill::zeros);
  arma::vec orig_Theta(ss_train.size(), arma::fill::zeros);
  compute_p_theta_ind_local(j, orig_P, orig_Theta, orig_leaf_map, ss_train, sigma, di_train);

  // Early reject if children too small (optional safety)
  auto nxl_chk = prop_ss_train.find(nxl_nid);
  auto nxr_chk = prop_ss_train.find(nxr_nid);
  if(nxl_chk == prop_ss_train.end() || nxr_chk == prop_ss_train.end() ||
     (int)nxl_chk->second.size() < MIN_LEAF || (int)nxr_chk->second.size() < MIN_LEAF){
    draw_mu_safe(t, orig_P, orig_Theta, orig_leaf_map, gen, tree_pi);
    accept = 0;
    return;
  }

  // Proposed likelihood
  std::map<int,int> prop_leaf_map;
  arma::mat prop_P (prop_ss_train.size(), prop_ss_train.size(), arma::fill::zeros);
  arma::vec prop_Theta(prop_ss_train.size(), arma::fill::zeros);
  compute_p_theta_ind_local(j, prop_P, prop_Theta, prop_leaf_map, prop_ss_train, sigma, di_train);

  const double orig_lil = compute_lil_safe(orig_P, orig_Theta, tree_pi);
  const double prop_lil = compute_lil_safe(prop_P, prop_Theta, tree_pi);

  const double log_like_ratio = prop_lil - orig_lil;
  const double log_alpha = log_like_ratio + log_prior_ratio + log_trans_ratio;

  if(gen.log_uniform() <= log_alpha){
    // bookkeeping for var/rule counts
    ++(*tree_pi.rule_count);
    if(rule.is_aa && !rule.is_cat){
      ++(tree_pi.var_count->at(rule.v_aa));
    } else if(!rule.is_aa && rule.is_cat){
      int v_raw = rule.v_cat + di_train.R_cont;
      ++(tree_pi.var_count->at(v_raw));
    } else if(!rule.is_aa && !rule.is_cat){
      ++(*tree_pi.rc_rule_count);
      for(auto it = rule.rc_weight.begin(); it != rule.rc_weight.end(); ++it){
        ++(*tree_pi.rc_var_count);
      }
    } else{
      Rcpp::stop("[grow_tree_ind]: unable to resolve rule type!");
    }

    // update suff stats maps
    auto nxl_it = prop_ss_train.find(nxl_nid);
    auto nxr_it = prop_ss_train.find(nxr_nid);
    if(nxl_it == prop_ss_train.end() || nxr_it == prop_ss_train.end()){
      Rcpp::stop("[grow_tree_ind]: missing child ids in proposed training suff-stat map");
    } else{
      ss_train.insert(std::make_pair(nxl_nid, nxl_it->second));
      ss_train.insert(std::make_pair(nxr_nid, nxr_it->second));
      ss_train.erase(nx_nid);
    }
    if(di_test.N > 0){
      nxl_it = prop_ss_test.find(nxl_nid);
      nxr_it = prop_ss_test.find(nxr_nid);
      if(nxl_it == prop_ss_test.end() || nxr_it == prop_ss_test.end()){
        Rcpp::stop("[grow_tree_ind]: missing child ids in proposed testing suff-stat map");
      } else{
        ss_test.insert(std::make_pair(nxl_nid, nxl_it->second));
        ss_test.insert(std::make_pair(nxr_nid, nxr_it->second));
        ss_test.erase(nx_nid);
      }
    }

    t.birth(nx_nid, rule);
    draw_mu_safe(t, prop_P, prop_Theta, prop_leaf_map, gen, tree_pi);
    accept = 1;
  } else{
    draw_mu_safe(t, orig_P, orig_Theta, orig_leaf_map, gen, tree_pi);
    accept = 0;
  }
}

void prune_tree_ind(tree &t, int &accept, int &j, double &sigma, suff_stat &ss_train, suff_stat &ss_test,
                    data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  tree::npv nogs_vec;
  t.get_nogs(nogs_vec);

  int ni = static_cast<int>(std::floor(gen.uniform() * nogs_vec.size()));
  tree::tree_p nx = nogs_vec[ni];
  tree::tree_p nxl = nx->get_l();
  tree::tree_p nxr = nx->get_r();

  const double q_prune_old = 1.0 - tree_pi.prob_b;
  const double q_grow_new  = (nx->get_p()==0 ? 1.0 : tree_pi.prob_b);

  const int nleaf_new = t.get_nbots() - 1;
  const int nnog_old  = t.get_nnogs();

  const double log_trans_ratio = (std::log(q_grow_new) - std::log((double)nleaf_new))
                               - (std::log(q_prune_old) - std::log((double)nnog_old));

  const double p_grow_nx  = tree_pi.alpha/std::pow(1.0 + (double) nx->get_depth(), tree_pi.beta);
  const double p_grow_nxl = tree_pi.alpha/std::pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  const double p_grow_nxr = tree_pi.alpha/std::pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  const double log_prior_ratio = std::log(1.0 - p_grow_nx)
                               - (std::log(1.0 - p_grow_nxl) + std::log(1.0 - p_grow_nxr) + std::log(p_grow_nx));

  // proposed suff stats after prune (children -> parent)
  suff_stat prop_ss_train;
  int nx_nid  = nx->get_nid();
  int nxl_nid = nxl->get_nid();
  int nxr_nid = nxr->get_nid();

  compute_suff_stat_prune(ss_train, prop_ss_train, nxl_nid, nxr_nid, nx_nid, t, di_train);

  suff_stat prop_ss_test;
  if(di_test.N > 0) compute_suff_stat_prune(ss_test, prop_ss_test, nxl_nid, nxr_nid, nx_nid, t, di_test);

  // Likelihood parts via x-weighted stats
  std::map<int,int> orig_leaf_map;
  arma::mat orig_P (ss_train.size(), ss_train.size(), arma::fill::zeros);
  arma::vec orig_Theta(ss_train.size(), arma::fill::zeros);
  compute_p_theta_ind_local(j, orig_P, orig_Theta, orig_leaf_map, ss_train, sigma, di_train);

  std::map<int,int> prop_leaf_map;
  arma::mat prop_P (prop_ss_train.size(), prop_ss_train.size(), arma::fill::zeros);
  arma::vec prop_Theta(prop_ss_train.size(), arma::fill::zeros);
  compute_p_theta_ind_local(j, prop_P, prop_Theta, prop_leaf_map, prop_ss_train, sigma, di_train);

  const double orig_lil = compute_lil_safe(orig_P, orig_Theta, tree_pi);
  const double prop_lil = compute_lil_safe(prop_P, prop_Theta, tree_pi);

  const double log_like_ratio = prop_lil - orig_lil;
  const double log_alpha = log_like_ratio + log_prior_ratio + log_trans_ratio;

  if(gen.log_uniform() <= log_alpha){
    --(*tree_pi.rule_count);
    if(nx->get_is_aa() && !nx->get_is_cat()){
      --(tree_pi.var_count->at(nx->get_v_aa()));
    } else if(!nx->get_is_aa() && nx->get_is_cat()){
      int v_raw = di_train.R_cont + nx->get_v_cat();
      --(tree_pi.var_count->at(v_raw));
    } else if(!nx->get_is_aa() && !nx->get_is_cat()){
      std::map<int,double> rc_weight = nx->get_rc_weight();
      --(*tree_pi.rc_rule_count);
      for(auto it = rc_weight.begin(); it != rc_weight.end(); ++it) --(*tree_pi.rc_var_count);
    } else{
      Rcpp::Rcout << "[prune tree]: cannot resolve rule type at nog node " << nx_nid << std::endl;
      t.print();
      Rcpp::stop("Cannot resolve rule type!");
    }

    // update suff stats maps
    auto nx_it = prop_ss_train.find(nx_nid);
    if(nx_it == prop_ss_train.end()){
      Rcpp::stop("[prune_tree_ind]: new leaf id not found in proposed training suff-stat map");
    } else{
      ss_train.erase(nxl_nid);
      ss_train.erase(nxr_nid);
      ss_train.insert(std::make_pair(nx_nid, nx_it->second));
    }

    if(di_test.N > 0){
      nx_it = prop_ss_test.find(nx_nid);
      if(nx_it == prop_ss_test.end()){
        Rcpp::stop("[prune_tree_ind]: new leaf id not found in proposed testing suff-stat map");
      } else{
        ss_test.erase(nxl_nid);
        ss_test.erase(nxr_nid);
        ss_test.insert(std::make_pair(nx_nid, nx_it->second));
      }
    }
    t.death(nx_nid);
    accept = 1;
    draw_mu_safe(t, prop_P, prop_Theta, prop_leaf_map, gen, tree_pi);
  } else{
    draw_mu_safe(t, orig_P, orig_Theta, orig_leaf_map, gen, tree_pi);
    accept = 0;
  }
}

void update_tree_ind(tree &t, int &accept, int &j, double &sigma, suff_stat &ss_train, suff_stat &ss_test,
                     data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  accept = 0;
  double PBx = tree_pi.prob_b;
  if(t.get_treesize() == 1) PBx = 1.0;
  if(gen.uniform() < PBx) grow_tree_ind(t, accept, j, sigma, ss_train, ss_test, di_train, di_test, tree_pi, gen);
  else prune_tree_ind(t, accept, j, sigma, ss_train, ss_test, di_train, di_test, tree_pi, gen);
}

// ---------------------------------------------------------------------
//  CS moves
// ---------------------------------------------------------------------

void grow_tree_cs(tree &t, int &accept, int &j, double &rho, double &sigma, suff_stat &ss_train, suff_stat &ss_test,
                  data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  std::vector<int> bn_nid_vec;
  bn_nid_vec.reserve(ss_train.size());
  for(const auto& kv : ss_train) bn_nid_vec.push_back(kv.first);

  int ni = static_cast<int>(std::floor(gen.uniform() * bn_nid_vec.size()));
  int nx_nid = bn_nid_vec[ni];
  int nxl_nid = 2*nx_nid;
  int nxr_nid = 2*nx_nid+1;

  tree::tree_p nx  = t.get_ptr(nx_nid);
  tree::tree_cp nxp = nx->get_p();

  double q_grow_old  = tree_pi.prob_b;
  double q_prune_new = 1.0 - tree_pi.prob_b;
  int nleaf_old = t.get_nbots();
  int nnog_old  = t.get_nnogs();
  int nnog_new  = nnog_old;

  if(nxp == 0) q_grow_old = 1.0;
  else if(!nxp->is_nog()) nnog_new = 1 + nnog_old;

  double log_trans_ratio = (std::log(q_prune_new) - std::log((double)nnog_new))
                         - (std::log(q_grow_old)  - std::log((double)nleaf_old));

  double p_grow_nx  = tree_pi.alpha/std::pow(1.0 + (double) nx->get_depth(), tree_pi.beta);
  double p_grow_nxl = tree_pi.alpha/std::pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  double p_grow_nxr = tree_pi.alpha/std::pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  double log_prior_ratio = std::log(p_grow_nx) + std::log(1.0 - p_grow_nxl) + std::log(1.0 - p_grow_nxr)
                         - std::log(1.0 - p_grow_nx);

  rule_t rule;
  draw_rule(rule, t, nx_nid, di_train, tree_pi, gen);

  suff_stat prop_ss_train;
  compute_suff_stat_grow(ss_train, prop_ss_train, nx_nid, rule, t, di_train);

  suff_stat prop_ss_test;
  if(di_test.N > 0) compute_suff_stat_grow(ss_test, prop_ss_test, nx_nid, rule, t, di_test);

  // x-weighted CS likelihood pieces (original)
  std::map<int,int> orig_leaf_map;
  arma::mat orig_P (ss_train.size(), ss_train.size(), arma::fill::zeros);
  arma::vec orig_Theta(ss_train.size(), arma::fill::zeros);
  compute_p_theta_cs_local(j, orig_P, orig_Theta, orig_leaf_map, ss_train, rho, sigma, di_train);

  // Early reject if children too small (optional safety)
  auto nxl_chk = prop_ss_train.find(nxl_nid);
  auto nxr_chk = prop_ss_train.find(nxr_nid);
  if(nxl_chk == prop_ss_train.end() || nxr_chk == prop_ss_train.end() ||
     (int)nxl_chk->second.size() < MIN_LEAF || (int)nxr_chk->second.size() < MIN_LEAF){
    draw_mu_safe(t, orig_P, orig_Theta, orig_leaf_map, gen, tree_pi);
    accept = 0;
    return;
  }

  // proposed
  std::map<int,int> prop_leaf_map;
  arma::mat prop_P (prop_ss_train.size(), prop_ss_train.size(), arma::fill::zeros);
  arma::vec prop_Theta(prop_ss_train.size(), arma::fill::zeros);
  compute_p_theta_cs_local(j, prop_P, prop_Theta, prop_leaf_map, prop_ss_train, rho, sigma, di_train);

  const double orig_lil = compute_lil_safe(orig_P, orig_Theta, tree_pi);
  const double prop_lil = compute_lil_safe(prop_P, prop_Theta, tree_pi);

  const double log_like_ratio = prop_lil - orig_lil;
  const double log_alpha = log_like_ratio + log_prior_ratio + log_trans_ratio;

  if(gen.log_uniform() <= log_alpha){
    ++(*tree_pi.rule_count);
    if(rule.is_aa && !rule.is_cat){
      ++(tree_pi.var_count->at(rule.v_aa));
    } else if(!rule.is_aa && rule.is_cat){
      int v_raw = rule.v_cat + di_train.R_cont;
      ++(tree_pi.var_count->at(v_raw));
    } else if(!rule.is_aa && !rule.is_cat){
      ++(*tree_pi.rc_rule_count);
      for(auto it = rule.rc_weight.begin(); it != rule.rc_weight.end(); ++it){
        ++(*tree_pi.rc_var_count);
      }
    } else{
      Rcpp::stop("[grow_tree_cs]: unable to resolve rule type!");
    }

    auto nxl_it = prop_ss_train.find(nxl_nid);
    auto nxr_it = prop_ss_train.find(nxr_nid);
    if(nxl_it == prop_ss_train.end() || nxr_it == prop_ss_train.end()){
      Rcpp::stop("[grow_tree_cs]: missing child ids in proposed training suff-stat map");
    } else{
      ss_train.insert(std::make_pair(nxl_nid, nxl_it->second));
      ss_train.insert(std::make_pair(nxr_nid, nxr_it->second));
      ss_train.erase(nx_nid);
    }
    if(di_test.N > 0){
      nxl_it = prop_ss_test.find(nxl_nid);
      nxr_it = prop_ss_test.find(nxr_nid);
      if(nxl_it == prop_ss_test.end() || nxr_it == prop_ss_test.end()){
        Rcpp::stop("[grow_tree_cs]: missing child ids in proposed testing suff-stat map");
      } else{
        ss_test.insert(std::make_pair(nxl_nid, nxl_it->second));
        ss_test.insert(std::make_pair(nxr_nid, nxr_it->second));
        ss_test.erase(nx_nid);
      }
    }

    t.birth(nx_nid, rule);
    draw_mu_safe(t, prop_P, prop_Theta, prop_leaf_map, gen, tree_pi);
    accept = 1;
  } else{
    draw_mu_safe(t, orig_P, orig_Theta, orig_leaf_map, gen, tree_pi);
    accept = 0;
  }
}

void prune_tree_cs(tree &t, int &accept, int &j, double &rho, double &sigma, suff_stat &ss_train, suff_stat &ss_test,
                   data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  tree::npv nogs_vec;
  t.get_nogs(nogs_vec);

  int ni = static_cast<int>(std::floor(gen.uniform() * nogs_vec.size()));
  tree::tree_p nx  = nogs_vec[ni];
  tree::tree_p nxl = nx->get_l();
  tree::tree_p nxr = nx->get_r();

  const double q_prune_old = 1.0 - tree_pi.prob_b;
  const double q_grow_new  = (nx->get_p()==0 ? 1.0 : tree_pi.prob_b);

  const int nleaf_new = t.get_nbots() - 1;
  const int nnog_old  = t.get_nnogs();

  const double log_trans_ratio = (std::log(q_grow_new) - std::log((double)nleaf_new))
                               - (std::log(q_prune_old) - std::log((double)nnog_old));

  const double p_grow_nx  = tree_pi.alpha/std::pow(1.0 + (double) nx->get_depth(), tree_pi.beta);
  const double p_grow_nxl = tree_pi.alpha/std::pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  const double p_grow_nxr = tree_pi.alpha/std::pow(2.0 + (double) nx->get_depth(), tree_pi.beta);
  const double log_prior_ratio = std::log(1.0 - p_grow_nx)
                               - (std::log(1.0 - p_grow_nxl) + std::log(1.0 - p_grow_nxr) + std::log(p_grow_nx));

  suff_stat prop_ss_train;
  int nx_nid  = nx->get_nid();
  int nxl_nid = nxl->get_nid();
  int nxr_nid = nxr->get_nid();
  compute_suff_stat_prune(ss_train, prop_ss_train, nxl_nid, nxr_nid, nx_nid, t, di_train);

  suff_stat prop_ss_test;
  if(di_test.N > 0) compute_suff_stat_prune(ss_test, prop_ss_test, nxl_nid, nxr_nid, nx_nid, t, di_test);

  // x-weighted CS likelihood pieces
  std::map<int,int> orig_leaf_map;
  arma::mat orig_P (ss_train.size(), ss_train.size(), arma::fill::zeros);
  arma::vec orig_Theta(ss_train.size(), arma::fill::zeros);
  compute_p_theta_cs_local(j, orig_P, orig_Theta, orig_leaf_map, ss_train, rho, sigma, di_train);

  std::map<int,int> prop_leaf_map;
  arma::mat prop_P (prop_ss_train.size(), prop_ss_train.size(), arma::fill::zeros);
  arma::vec prop_Theta(prop_ss_train.size(), arma::fill::zeros);
  compute_p_theta_cs_local(j, prop_P, prop_Theta, prop_leaf_map, prop_ss_train, rho, sigma, di_train);

  const double orig_lil = compute_lil_safe(orig_P, orig_Theta, tree_pi);
  const double prop_lil = compute_lil_safe(prop_P, prop_Theta, tree_pi);

  const double log_like_ratio = prop_lil - orig_lil;
  const double log_alpha = log_like_ratio + log_prior_ratio + log_trans_ratio;

  if(gen.log_uniform() <= log_alpha){
    --(*tree_pi.rule_count);
    if(nx->get_is_aa() && !nx->get_is_cat()){
      --(tree_pi.var_count->at(nx->get_v_aa()));
    } else if(!nx->get_is_aa() && nx->get_is_cat()){
      int v_raw = di_train.R_cont + nx->get_v_cat();
      --(tree_pi.var_count->at(v_raw));
    } else if(!nx->get_is_aa() && !nx->get_is_cat()){
      std::map<int,double> rc_weight = nx->get_rc_weight();
      --(*tree_pi.rc_rule_count);
      for(auto it = rc_weight.begin(); it != rc_weight.end(); ++it) --(*tree_pi.rc_var_count);
    } else{
      Rcpp::Rcout << "[prune_tree_cs]: cannot resolve rule type at nog node " << nx_nid << std::endl;
      t.print();
      Rcpp::stop("Cannot resolve rule type!");
    }

    auto nx_it = prop_ss_train.find(nx_nid);
    if(nx_it == prop_ss_train.end()){
      Rcpp::stop("[prune_tree_cs]: new leaf id not found in proposed training suff-stat map");
    } else{
      ss_train.erase(nxl_nid);
      ss_train.erase(nxr_nid);
      ss_train.insert(std::make_pair(nx_nid, nx_it->second));
    }

    if(di_test.N > 0){
      nx_it = prop_ss_test.find(nx_nid);
      if(nx_it == prop_ss_test.end()){
        Rcpp::stop("[prune_tree_cs]: new leaf id not found in proposed testing suff-stat map");
      } else{
        ss_test.erase(nxl_nid);
        ss_test.erase(nxr_nid);
        ss_test.insert(std::make_pair(nx_nid, nx_it->second));
      }
    }
    t.death(nx_nid);
    accept = 1;
    draw_mu_safe(t, prop_P, prop_Theta, prop_leaf_map, gen, tree_pi);
  } else{
    draw_mu_safe(t, orig_P, orig_Theta, orig_leaf_map, gen, tree_pi);
    accept = 0;
  }
}

void update_tree_cs(tree &t, int &accept, int &j, double &rho, double &sigma, suff_stat &ss_train, suff_stat &ss_test,
                    data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  accept = 0;
  double PBx = tree_pi.prob_b;
  if(t.get_treesize() == 1) PBx = 1.0;

  if(gen.uniform() < PBx) grow_tree_cs(t, accept, j, rho, sigma, ss_train, ss_test, di_train, di_test, tree_pi, gen);
  else prune_tree_cs(t, accept, j, rho, sigma, ss_train, ss_test, di_train, di_test, tree_pi, gen);
}
