#include <Rcpp.h>
using namespace Rcpp;

//' Calculate a Sample Vector
//'
//' This function takes a vector of probabilities and an integer value `m`, and optionally a seed,
//' and returns a sampled vector based on the provided probabilities.
//'
//' @param probs A numeric vector of probabilities.
//' @param m An integer indicating the number of samples to return.
//' @param seed An optional integer that sets the random seed for reproducibility.
//' @return A numeric vector of sampled values.
//' @examples
//' samplev_cpp(c(0.1, 0.9), 5)
//' @export
// [[Rcpp::export]]
SEXP samplev_cpp(SEXP probs, int m, int seed=NA_INTEGER) {
  if (seed != NA_INTEGER) {
    Environment base_env("package:base");
    Function set_seed = base_env["set.seed"];
    set_seed(seed);
  }
  
  if (Rf_isMatrix(probs)) {
    NumericMatrix matProbs(probs);
    int n = matProbs.nrow();
    int k = matProbs.ncol();
    CharacterMatrix ran(n, m);
    CharacterVector lev = colnames(matProbs);
    
    for (int i = 0; i < n; ++i) {
      NumericVector U = cumsum(matProbs(i, _));
      if (std::abs(U[k-1] - 1) > 1e-5) {
        stop("error in multinom; probabilities do not sum to 1");
      }
      for (int j = 0; j < m; ++j) {
        double un = R::runif(0, 1);
        int idx = std::upper_bound(U.begin(), U.end(), un) - U.begin();
        ran(i, j) = lev[idx];
      }
    }
    return ran;
  } else {
    NumericVector vecProbs(probs);
    int k = vecProbs.size();
    NumericVector U = cumsum(vecProbs);
    CharacterVector lev = as<CharacterVector>(vecProbs.names());
    if (lev.size() == 0) {
      lev = CharacterVector(k);
      for (int i = 0; i < k; ++i) lev[i] = std::to_string(i + 1);
    }
    CharacterVector ran(m);
    
    if (std::abs(U[k-1] - 1) > 1e-5) {
      stop("error in multinom; probabilities do not sum to 1");
    }
    for (int j = 0; j < m; ++j) {
      double un = R::runif(0, 1);
      int idx = std::upper_bound(U.begin(), U.end(), un) - U.begin();
      ran[j] = lev[idx];
    }
    return ran;
  }
}