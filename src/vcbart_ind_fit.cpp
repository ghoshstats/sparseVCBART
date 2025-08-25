#include "update_trees.h"
#include "funs.h"
#include <cmath>
#include <algorithm>
#include <vector>


// Compute sum of squared leaf means and number of leaves for each predictor j
inline void compute_leaf_summaries(std::vector<double>& S, std::vector<int>& L,
                                   const std::vector<std::vector<tree>>& tree_vec) {
  const int p = (int)tree_vec.size();
  for(int j = 0; j < p; ++j) {
    double Sj = 0.0;
    int    Lj = 0;
    for(size_t m = 0; m < tree_vec[j].size(); ++m) {
      tree::cnpv bnv;                      // const leaf pointers
      tree_vec[j][m].get_bots(bnv);
      for(tree::cnpv_it it = bnv.begin(); it != bnv.end(); ++it) {
        const double mu = (*it)->get_mu();
        Sj += mu * mu;
        Lj += 1;
      }
    }
    S[j] = Sj;
    L[j] = Lj;
  }
}

// Map (τ, λ_j, c2) to per-leaf prior sd injected into tree prior:
//   tau_leaf_j = sqrt( (τ^2 λ_j^2) * c2 / (c2 + τ^2 λ_j^2) ) / sqrt(M)
inline double leaf_tau_from_gl(double tau, double lambda_j, double c2, int M) {
  const double tl2 = tau * tau * lambda_j * lambda_j;
  const double s2  = (tl2 * c2) / (c2 + tl2);
  const double s_leaf = std::sqrt( (s2 > 1e-20 ? s2 : 1e-20) ) / std::sqrt((double)M);
  return s_leaf;
}

// ----------------------------


// [[Rcpp::export(".vcbart_ind_fit")]]
Rcpp::List vcbart_ind_fit(Rcpp::NumericVector Y_train,
                          Rcpp::IntegerVector subj_id_train,
                          Rcpp::IntegerVector ni_train,
                          Rcpp::NumericMatrix tX_train,
                          Rcpp::NumericMatrix tZ_cont_train,
                          Rcpp::IntegerMatrix tZ_cat_train,
                          Rcpp::NumericMatrix tX_test,
                          Rcpp::NumericMatrix tZ_cont_test,
                          Rcpp::IntegerMatrix tZ_cat_test,
                          Rcpp::LogicalVector unif_cuts,
                          Rcpp::Nullable<Rcpp::List> cutpoints_list,
                          Rcpp::Nullable<Rcpp::List> cat_levels_list,
                          Rcpp::Nullable<Rcpp::List> edge_mat_list,
                          Rcpp::LogicalVector graph_split, int graph_cut_type,
                          bool rc_split, double prob_rc, double a_rc, double b_rc,
                          bool sparse, double a_u, double b_u,
                          Rcpp::NumericVector mu0, Rcpp::NumericVector tau,
                          double lambda, double nu,
                          int M,
                          int nd, int burn, int thin,
                          bool save_samples,
                          bool save_trees,
                          bool verbose, int print_every)
{
  Rcpp::RNGScope scope;
  RNG gen;

  set_str_conversion set_str;

  // -------------------------
  // Preprocessing
  // -------------------------
  int N_train = 0;
  int n_train = 0;
  int p = 0;
  int R_cont = 0;
  int R_cat = 0;
  int R = 0;

  int N_test = 0;

  parse_training_data(N_train, n_train, p, R, R_cont, R_cat, ni_train, tX_train, tZ_cont_train, tZ_cat_train);
  if(Y_train.size() != N_train) Rcpp::stop("Number of observations in Y_train does not match number of rows in matrix of training X's");
  parse_testing_data(N_test, p, R_cont, R_cat, tX_test, tZ_cont_test, tZ_cat_test);

  if(verbose){
    Rcpp::Rcout << " p = " << p << " R_cont = " << R_cont << " R_cat = " << R_cat << std::endl;
    Rcpp::Rcout << "N_train = " << N_train << " n_train = " << n_train << std::endl;
    Rcpp::Rcout << "N_test = " << N_test << std::endl;
  }

  std::vector<std::set<double>> cutpoints;
  if(R_cont > 0){
    if(cutpoints_list.isNotNull()){
      Rcpp::List tmp_cutpoints = Rcpp::List(cutpoints_list);
      parse_cutpoints(cutpoints, R_cont, tmp_cutpoints, unif_cuts);
    }
  }

  std::vector<std::set<int>> cat_levels;
  std::vector<int> K;
  std::vector<std::vector<edge>> edges;

  if(R_cat > 0){
    if(cat_levels_list.isNotNull()){
      Rcpp::List tmp_cat_levels = Rcpp::List(cat_levels_list);
      parse_cat_levels(cat_levels,K, R_cat, tmp_cat_levels);
    } else{
      Rcpp::stop("Must provide list of categorical levels");
    }
    if(edge_mat_list.isNotNull()){
      Rcpp::List tmp_edge_mat = Rcpp::List(edge_mat_list);
      parse_graphs(edges, R_cat, K, tmp_edge_mat, graph_split);
    }
  }

  double* allfit_train = new double[N_train];
  double* beta_fit_train = new double[p*N_train];
  double* residual = new double[N_train];

  double* r_sum = new double[n_train];
  double* r2_sum = new double[n_train];

  for(int subj_ix = 0; subj_ix < n_train; subj_ix++){
    r_sum[subj_ix] = 0.0;
    r2_sum[subj_ix] = 0.0;
  }

  for(int i = 0; i < N_train; i++){
    allfit_train[i] = 0.0;
    for(int j = 0; j < p; j++) beta_fit_train[j + i*p] = 0.0;
  }

  data_info di_train;
  di_train.N = N_train;
  di_train.n = n_train;
  di_train.p = p;
  di_train.R = R;
  di_train.R_cont = R_cont;
  di_train.R_cat = R_cat;
  di_train.subj_id = subj_id_train.begin();
  di_train.ni = ni_train.begin();
  di_train.x = tX_train.begin();
  if(R_cont > 0) di_train.z_cont = tZ_cont_train.begin();
  if(R_cat > 0) di_train.z_cat = tZ_cat_train.begin();
  di_train.rp = residual;
  di_train.r_sum = r_sum;
  di_train.r2_sum = r2_sum;

  int tmp_N_test = (N_test > 0 ? N_test : 1);
  double* allfit_test = new double[tmp_N_test];
  double* beta_fit_test = new double[p*tmp_N_test];
  if(N_test > 0){
    for(int ii = 0; ii < N_test; ii++){
      allfit_test[ii] = 0.0;
      for(int j = 0; j < p; j++) beta_fit_test[j + ii*p] = 0.0;
    }
  }

  data_info di_test;
  if(N_test > 0){
    di_test.N = N_test;
    di_test.p = p;
    di_test.R = R;
    di_test.R_cont = R_cont;
    di_test.R_cat = R_cat;
    di_test.x = tX_test.begin();
    if(R_cont > 0) di_test.z_cont = tZ_cont_test.begin();
    if(R_cat > 0) di_test.z_cat = tZ_cat_test.begin();
  }

  // -------------------------
  // Modifier selection parts
  // -------------------------
  std::vector<std::vector<double>> theta(p, std::vector<double>(R, 1.0/((double)R)));
  double* u = new double[p];
  for(int j = 0; j < p; j++) u[j] = 1.0/(1.0 + (double) R);
  std::vector<std::vector<int>> var_count(p, std::vector<int>(R, 0));
  int* rule_count = new int[p];
  for(int j = 0; j < p; j++) rule_count[j] = 0;
  int* rc_rule_count = new int[p];
  int* rc_var_count  = new int[p];
  double* theta_rc   = new double[p];
  if(R_cont >= 2 && rc_split){
    for(int j = 0; j < p; j++) theta_rc[j] = 2.0/((double) R_cont);
  }

  // -------------------------
  // Tree priors (one per predictor ensemble)
  // -------------------------
  std::vector<tree_prior_info> tree_pi_vec(p);
  for(int j = 0; j < p; j++){
    tree_pi_vec[j].theta      = &(theta[j]);
    tree_pi_vec[j].var_count  = &(var_count[j]);
    tree_pi_vec[j].rule_count = &(rule_count[j]);
    tree_pi_vec[j].unif_cuts  = unif_cuts.begin();

    if(R_cont > 0){
      tree_pi_vec[j].unif_cuts = unif_cuts.begin();
      tree_pi_vec[j].cutpoints = &cutpoints;
      if(rc_split){
        tree_pi_vec[j].rc_split      = rc_split;
        tree_pi_vec[j].prob_rc       = &prob_rc;
        tree_pi_vec[j].theta_rc      = &(theta_rc[j]);
        tree_pi_vec[j].rc_var_count  = &(rc_var_count[j]);
        tree_pi_vec[j].rc_rule_count = &(rc_rule_count[j]);
      }
    }
    if(R_cat > 0){
      tree_pi_vec[j].cat_levels = &cat_levels;
      tree_pi_vec[j].edges      = &edges;
      tree_pi_vec[j].K          = &K;
      tree_pi_vec[j].graph_split= graph_split.begin();
    }
    tree_pi_vec[j].mu0 = 0.0;                         // GL center at 0
    tree_pi_vec[j].tau = 1.0 / std::sqrt((double)M);  // placeholder; overwritten below
  }

  // -------------------------
  // Residuals & GL state (independent errors)
  // -------------------------
  double sigma = 1.0;

  // Regularized HS hyperparameters
  const double slab_df    = 4.0;
  const double slab_scale = 2.0;

  // τ prior scale (Piironen & Vehtari): tau0 ≈ (p0/(p-p0)) * sd(y)/sqrt(N)
  int p0 = std::min(10, std::max(1, p/4));
  double y_mean = 0.0;
  for(int ii=0; ii<N_train; ++ii) y_mean += Y_train[ii];
  y_mean /= (double)N_train;
  double y_ss = 0.0;
  for(int ii=0; ii<N_train; ++ii){ double d = Y_train[ii]-y_mean; y_ss += d*d; }
  const double y_sd  = std::sqrt(y_ss / (double)(N_train-1));
  const double tau0  = ( (double)p0 / std::max(1.0, (double)(p - p0)) ) * (y_sd/std::sqrt((double)N_train));

  // Initialize τ, λ, c2
  double tau_gl = std::max(1e-2, tau0);
  std::vector<double> lambda_loc(p, 1.0);
  double c2 = slab_scale * slab_scale;

  // Inject initial leaf prior SDs
  for(int j = 0; j < p; ++j){
    tree_pi_vec[j].tau = leaf_tau_from_gl(tau_gl, lambda_loc[j], c2, M);
  }

  // -------------------------
  // Init trees & running sums
  // -------------------------
  int i = 0;
  double tmp_mu;
  std::vector<std::vector<tree>>      tree_vec(p, std::vector<tree>(M));
  std::vector<std::vector<suff_stat>> ss_train_vec(p, std::vector<suff_stat>(M));
  std::vector<std::vector<suff_stat>> ss_test_vec(p, std::vector<suff_stat>(M));

  for(int j = 0; j < p; j++){
    for(int m = 0; m < M; m++){
      tree_traversal(ss_train_vec[j][m], tree_vec[j][m], di_train);
      for(suff_stat_it l_it = ss_train_vec[j][m].begin(); l_it != ss_train_vec[j][m].end(); ++l_it){
        tmp_mu = tree_vec[j][m].get_ptr(l_it->first)->get_mu();
        for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
          i = *it;
          allfit_train[i] += tX_train(j,i) * tmp_mu;
          beta_fit_train[j + i*p] += tmp_mu;
        }
      }
      if(N_test > 0){
        tree_traversal(ss_test_vec[j][m], tree_vec[j][m], di_test);
        for(suff_stat_it l_it = ss_test_vec[j][m].begin(); l_it != ss_test_vec[j][m].end(); ++l_it){
          tmp_mu = tree_vec[j][m].get_ptr(l_it->first)->get_mu();
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
            i = *it;
            allfit_test[i] += tX_test(j,i) * tmp_mu;
            beta_fit_test[j+i*p] += tmp_mu;
          }
        }
      }
    }
  }

  for(int ii = 0; ii < N_train; ii++){
    residual[ii] = Y_train[ii] - allfit_train[ii];
    r_sum[subj_id_train[ii]]  += residual[ii];
    r2_sum[subj_id_train[ii]] += residual[ii]*residual[ii];
  }

  // -------------------------
  // Output containers
  // -------------------------
  int total_draws = 1 + burn + (nd - 1) * thin;
  int sample_index = 0;
  int accept = 0;
  int* total_accept = new int[p];
  for(int j = 0; j < p; j++) total_accept[j] = 0;

  arma::vec fit_train_mean = arma::zeros<arma::vec>(N_train);
  arma::vec fit_test_mean  = arma::zeros<arma::vec>(tmp_N_test);
  arma::mat beta_train_mean= arma::zeros<arma::mat>(N_train,p);
  arma::mat beta_test_mean = arma::zeros<arma::mat>(tmp_N_test,p);

  arma::mat  fit_train = arma::zeros<arma::mat>(1,1);
  arma::mat  fit_test  = arma::zeros<arma::mat>(1,1);
  arma::cube beta_train= arma::zeros<arma::cube>(1,1,1);
  arma::cube beta_test = arma::zeros<arma::cube>(1,1,1);

  if(save_samples){
    fit_train.zeros(nd, N_train);
    beta_train.zeros(nd, N_train, p);
    if(N_test > 0){
      fit_test.zeros(nd, N_test);
      beta_test.zeros(nd, N_test, p);
    }
  }

  arma::vec sigma_samples(total_draws);
  arma::mat total_accept_samples(nd,p);
  arma::cube theta_samples(1,1,1);
  if(sparse) theta_samples.zeros(total_draws, R, p);
  arma::cube var_count_samples(total_draws, R, p);

  // GL/RHS outputs
  arma::vec tau_samples(total_draws);
  arma::vec c2_samples(total_draws);
  arma::mat lambda_samples(total_draws, p);

  Rcpp::List tree_draws(nd);

  // -------------------------
  // Main MCMC
  // -------------------------
  for(int iter = 0; iter < total_draws; iter++){
    if(verbose){
      if( (iter < burn) && (iter % print_every == 0)){
        Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws << "; Warmup" << std::endl;
        Rcpp::checkUserInterrupt();
      } else if(((iter> burn) && (iter%print_every == 0)) || (iter == burn) ){
        Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws << "; Sampling" << std::endl;
        Rcpp::checkUserInterrupt();
      } else if( iter == total_draws-1){
        Rcpp::Rcout << "  MCMC Iteration: " << iter+1 << " of " << total_draws << "; Sampling" << std::endl;
        Rcpp::checkUserInterrupt();
      }
    }

    if( (iter == burn) && (N_test > 0)){
      for(int j = 0; j < p; j++){
        for(int m = 0; m < M; m++){
          for(suff_stat_it l_it = ss_test_vec[j][m].begin(); l_it != ss_test_vec[j][m].end(); ++l_it){
            double mu_tmp = tree_vec[j][m].get_ptr(l_it->first)->get_mu();
            for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
              int ii = *it;
              allfit_test[ii] += tX_test(j,ii) * mu_tmp;
              beta_fit_test[j + ii*p] += mu_tmp;
            }
          }
        }
      }
    }

    // ----- Loop over ensembles (tree updates; independent errors) -----
    for(int j = 0; j < p; j++){
      total_accept[j] = 0;
      for(int m = 0; m < M; m++){
        // remove fit of tree (j,m)
        for(suff_stat_it l_it = ss_train_vec[j][m].begin(); l_it != ss_train_vec[j][m].end(); ++l_it){
          double mu_tmp = tree_vec[j][m].get_ptr(l_it->first)->get_mu();
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
            int ii = *it;
            allfit_train[ii] -= tX_train(j,ii) * mu_tmp;
            residual[ii]     += tX_train(j,ii) * mu_tmp;
            if(std::abs(residual[ii] - (Y_train[ii] - allfit_train[ii])) > 1e-12)
              Rcpp::stop("after removing fit of single tree, something is wrong with residual");
            beta_fit_train[j + ii*p] -= mu_tmp;
            r_sum[subj_id_train[ii]]  += tX_train(j,ii) * mu_tmp;
            r2_sum[subj_id_train[ii]] += 2.0 * tX_train(j,ii) * mu_tmp * residual[ii] - std::pow(tX_train(j,ii) * mu_tmp, 2.0);
          }
        }
        if(iter >= burn && N_test > 0){
          for(suff_stat_it l_it = ss_test_vec[j][m].begin(); l_it != ss_test_vec[j][m].end(); ++l_it){
            double mu_tmp = tree_vec[j][m].get_ptr(l_it->first)->get_mu();
            for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
              int ii = *it;
              allfit_test[ii] -= tX_test(j,ii) * mu_tmp;
              beta_fit_test[j + ii*p] -= mu_tmp;
            }
          }
        }

        // update tree (independent errors)
        update_tree_ind(tree_vec[j][m], accept, j, sigma,
                        ss_train_vec[j][m], ss_test_vec[j][m],
                        di_train, di_test, tree_pi_vec[j], gen);
        total_accept[j] += accept;

        // add fit back
        for(suff_stat_it l_it = ss_train_vec[j][m].begin(); l_it != ss_train_vec[j][m].end(); ++l_it){
          double mu_tmp = tree_vec[j][m].get_ptr(l_it->first)->get_mu();
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
            int ii = *it;
            r2_sum[subj_id_train[ii]] += std::pow(tX_train(j,ii) * mu_tmp, 2.0) - 2.0 * tX_train(j,ii) * mu_tmp * residual[ii];
            r_sum[subj_id_train[ii]]  -= tX_train(j,ii) * mu_tmp;
            allfit_train[ii] += tX_train(j,ii) * mu_tmp;
            residual[ii]     -= tX_train(j,ii) * mu_tmp;
            if(std::abs(residual[ii] - (Y_train[ii] - allfit_train[ii])) > 1e-12)
              Rcpp::stop("after restoring fit of single tree, something is wrong with residual");
            beta_fit_train[j + ii*p] += mu_tmp;
          }
        }
        if(iter >= burn && N_test > 0){
          for(suff_stat_it l_it = ss_test_vec[j][m].begin(); l_it != ss_test_vec[j][m].end(); ++l_it){
            double mu_tmp = tree_vec[j][m].get_ptr(l_it->first)->get_mu();
            for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
              int ii = *it;
              allfit_test[ii] += tX_test(j,ii) * mu_tmp;
              beta_fit_test[j + ii*p] += mu_tmp;
            }
          }
        }
      } // m over trees

      if(sparse){
        update_theta_u(theta[j], u[j], var_count[j], R, a_u, b_u, gen);
        for(int r = 0; r < R; r++){
          theta_samples(iter, r, j) = theta[j][r];
        }
      }
    } // j over ensembles

    // update var_count_samples
    for(int j = 0; j < p; j++){
      for(int r = 0; r < R; r++){
        var_count_samples(iter, r, j) = var_count[j][r];
      }
    }

    // ----- Update residual σ (independent) -----
    update_sigma_ind(sigma, nu, lambda, di_train, gen);
    sigma_samples(iter) = sigma;

    // ----- Global–local (RHS) updates via slice sampling on log-scale -----
    std::vector<double> S_j(p, 0.0);
    std::vector<int>    L_j(p, 0);
    compute_leaf_summaries(S_j, L_j, tree_vec);

    //  λ_j
    for(int j = 0; j < p; ++j){
      auto logf_lam = [&](double loglam){
        return logpost_loglambda_rhs(loglam, /*tau=*/tau_gl, /*c2=*/c2,
                                     /*sumsq_mu_j=*/S_j[j], /*n_mu_j=*/L_j[j], /*M=*/M);
      };
      double loglam0 = std::log(lambda_loc[j]);
      double loglam1 = slice1d_on_log(logf_lam, loglam0, gen,
                                      /*w=*/0.5, /*m=*/50, /*lower=*/-12.0, /*upper=*/12.0);
      lambda_loc[j] = std::exp(loglam1);
    }

    //  τ
    {
      auto logf_tau = [&](double logtau){
        return logpost_logtau_rhs(logtau, S_j, L_j, lambda_loc,
                                  /*c2=*/c2, /*M=*/M, /*tau0=*/tau0);
      };
      double logtau0 = std::log(tau_gl);
      double logtau1 = slice1d_on_log(logf_tau, logtau0, gen,
                                      /*w=*/0.5, /*m=*/50, /*lower=*/-12.0, /*upper=*/12.0);
      tau_gl = std::exp(logtau1);
    }

    //  c2 with IG(ν/2, (ν s^2)/2) prior
    {
      const double a_slab = 0.5 * slab_df;
      const double b_slab = 0.5 * slab_df * slab_scale * slab_scale;
      auto logf_c2 = [&](double logc2){
        return logpost_logc2_rhs(logc2, S_j, L_j, lambda_loc,
                                 /*tau=*/tau_gl, /*M=*/M, a_slab, b_slab);
      };
      double logc20 = std::log(c2);
      double logc21 = slice1d_on_log(logf_c2, logc20, gen,
                                     /*w=*/0.5, /*m=*/50, /*lower=*/-12.0, /*upper=*/12.0);
      c2 = std::exp(logc21);
    }

    
    for(int j = 0; j < p; ++j){
      tree_pi_vec[j].tau = leaf_tau_from_gl(tau_gl, lambda_loc[j], c2, M);
    }

    // Save GL params
    tau_samples(iter) = tau_gl;
    c2_samples(iter)  = c2;
    for(int j = 0; j < p; ++j) lambda_samples(iter, j) = lambda_loc[j];

    // ----- Save posterior draws of fits/betas -----
    if( (iter >= burn) && ( (iter - burn) % thin == 0)){
      sample_index = (int) ((iter - burn)/thin);
      for(int j = 0; j < p; j++) total_accept_samples(sample_index,j) = total_accept[j];

      if(save_trees){
        Rcpp::List tmp_tree_draws(p);
        for(int j = 0; j < p; j++){
          Rcpp::CharacterVector tree_string_vec(M);
          for(int m = 0; m < M; m++){
            tree_string_vec[m] = write_tree(tree_vec[j][m], tree_pi_vec[j], set_str);
          }
          tmp_tree_draws[j] = tree_string_vec;
        }
        tree_draws[sample_index] = tmp_tree_draws;
      }

      if(save_samples){
        for(int ii = 0; ii < N_train; ii++){
          fit_train(sample_index,ii) = allfit_train[ii];
          fit_train_mean(ii)        += allfit_train[ii];
          for(int j = 0; j < p; j++){
            beta_train(sample_index,ii,j) = beta_fit_train[j + ii * p];
            beta_train_mean(ii,j)        += beta_fit_train[j + ii*p];
          }
        }
        if(N_test > 0){
          for(int ii = 0; ii < N_test; ii++){
            fit_test(sample_index,ii) = allfit_test[ii];
            fit_test_mean(ii)        += allfit_test[ii];
            for(int j = 0; j < p; j++){
              beta_test(sample_index,ii,j) = beta_fit_test[j + ii*p];
              beta_test_mean(ii,j)        += beta_fit_test[j + ii*p];
            }
          }
        }
      } else{
        for(int ii = 0; ii < N_train; ii++){
          fit_train_mean(ii) += allfit_train[ii];
          for(int j = 0; j < p; j++){
            beta_train_mean(ii,j) += beta_fit_train[j + ii*p];
          }
        }
        if(N_test > 0){
          for(int ii = 0; ii < N_test; ii++){
            fit_test_mean(ii) += allfit_test[ii];
            for(int j = 0; j < p; j++){
              beta_test_mean(ii,j) += beta_fit_test[j + ii*p];
            }
          }
        }
      }
    } // save block
  } // end MCMC

  fit_train_mean /= ( (double) nd);
  beta_train_mean/= ( (double) nd);
  if(N_test > 0){
    fit_test_mean  /= ( (double) nd);
    beta_test_mean /= ( (double) nd);
  }

  Rcpp::List results;
  results["fit_train_mean"] = fit_train_mean;
  results["beta_train_mean"] = beta_train_mean;
  if(save_samples){
    results["fit_train"] = fit_train;
    results["beta_train"] = beta_train;
  }
  if(N_test > 0){
    results["fit_test_mean"] = fit_test_mean;
    results["beta_test_mean"] = beta_test_mean;
    if(save_samples){
      results["fit_test"] = fit_test;
      results["beta_test"] = beta_test;
    }
  }
  results["sigma"] = sigma_samples;
  results["total_accept"] = total_accept_samples;
  results["var_count"]    = var_count_samples;
  if(save_trees) results["trees"] = tree_draws;
  if(sparse)     results["theta"] = theta_samples;

  results["tau"]    = tau_samples;
  results["c2"]     = c2_samples;
  results["lambda"] = lambda_samples;

  delete[] allfit_train;
  delete[] allfit_test;
  delete[] beta_fit_train;
  delete[] beta_fit_test;
  delete[] residual;
  delete[] r_sum;
  delete[] r2_sum;
  delete[] total_accept;
  delete[] rule_count;
  delete[] rc_rule_count;
  delete[] rc_var_count;
  delete[] theta_rc;
  delete[] u;

  return results;
}
