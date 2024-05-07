//--------------------------------------------------------------
// Header (header)
//--------------------------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;


//--------------------------------------------------------------
// Functions (Functions_cpp)
//--------------------------------------------------------------
// [[Rcpp::export]]
List construct_bags(NumericMatrix xy_receiver, NumericMatrix xy_sender, CharacterVector cell_ids, double dist_cutoff2, int min_instance) {
  int n_receiver = xy_receiver.nrow();
  int n_sender = xy_sender.nrow();
  List pos_sender;
  List exp_sender;
  CharacterVector pos_names;
  CharacterVector exp_names;

  for (int i = 0; i < n_receiver; ++i) {
    LogicalVector keep(n_sender);
    NumericVector dists(n_sender);
    for (int j = 0; j < n_sender; ++j) {
      double dist = pow(xy_receiver(i, 0) - xy_sender(j, 0), 2) +
                    pow(xy_receiver(i, 1) - xy_sender(j, 1), 2);
      keep[j] = dist < dist_cutoff2;
      dists[j] = dist;
    }
    
    if (sum(keep) < min_instance) continue; // Skip if there aren't enough instances
    
    NumericVector valid_dists = sqrt(dists[keep]); // Only compute sqrt for valid distances
    pos_sender.push_back(log(valid_dists));        // Store the logarithm of valid distances
    exp_sender.push_back(keep);                    // Store the indices of valid senders
    pos_names.push_back(cell_ids[i]);              // Assign cell IDs
    exp_names.push_back(cell_ids[i]);
  }

  // Setting names to the lists after the loop
  pos_sender.attr("names") = pos_names;
  exp_sender.attr("names") = exp_names;

  return List::create(Named("pos_sender") = pos_sender, Named("exp_sender") = exp_sender);
}
