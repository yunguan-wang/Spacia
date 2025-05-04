//--------------------------------------------------------------
// Header (header)
//--------------------------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <cmath>
#include <algorithm>
#include <random>
using namespace Rcpp;

//--------------------------------------------------------------
// Functions (Functions_cpp)
//--------------------------------------------------------------
// [[Rcpp::export]]
List construct_bags(NumericMatrix xy_receiver, NumericMatrix xy_sender, CharacterVector cell_ids, 
                    double dist_cutoff2, int min_instance, int subSample) {
  int n_receiver = xy_receiver.nrow();
  int n_sender = xy_sender.nrow();
  List pos_sender;
  List exp_sender;
  CharacterVector pos_names;
  CharacterVector exp_names;
  
  // Determine the target number of bags
  int target_bags = (subSample > 0 && subSample < n_receiver) ? subSample : n_receiver;
  
  // Warning if subSample is greater than n_receiver
  if (subSample > n_receiver) {
    Rcpp::warning("subSample (%d) is greater than the number of receivers (%d). Using all receivers.", subSample, n_receiver);
  }
  
  // Create a randomly shuffled index vector
  std::vector<int> indices(n_receiver);
  for (int i = 0; i < n_receiver; ++i) indices[i] = i;
  std::random_shuffle(indices.begin(), indices.end());
  
  Rcout << "Attempting to construct up to " << target_bags << " bags.\n";
  
  // Setup progress bar
  Progress p(n_receiver, true);
  
  int bags_constructed = 0;
  int receivers_checked = 0;
  
  while (bags_constructed < target_bags && receivers_checked < n_receiver) {
    // Check for user interrupt and update progress
    if (Progress::check_abort()) {
      throw Rcpp::exception("User interrupted the bag construction process.");
    }
    p.increment();
    
    int i = indices[receivers_checked];
    receivers_checked++;
    
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
    
    bags_constructed++;
  }
  
  // Setting names to the lists after the loop
  pos_sender.attr("names") = pos_names;
  exp_sender.attr("names") = exp_names;
  
  // Notify the user about the number of bags constructed and receivers checked
  Rcout << "Successfully constructed " << bags_constructed << " bags.\n";
  Rcout << "Checked " << receivers_checked << " out of " << n_receiver << " receivers.\n";
  
  return List::create(Named("pos_sender") = pos_sender, Named("exp_sender") = exp_sender);
}
