#include "DPBWT.hpp"
#include <math.h>
#include <vector>
#include <memory>
#include <cassert>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>


using std::cout;
using std::cerr;
using std::endl;

int END_ALLELE = 0;


DPBWT::DPBWT(std::vector<double> _physical_positions, 
             std::vector<double> _genetic_positions, 
             double _mutation_rate, 
             std::vector<double> ne, 
             std::vector<double> ne_times, 
             bool _sparse_sites,
             int _n_prune) :
  physical_positions(_physical_positions), genetic_positions(_genetic_positions), mutation_rate(_mutation_rate), demography(Demography(ne, ne_times)), sparse_sites(_sparse_sites), n_prune(_n_prune)
{
  if (physical_positions.size() != genetic_positions.size()) {
    cerr << "Map lengths don't match.\n";
    exit(1);
  } else if (physical_positions.size() <= 2) {
    cerr << "Need at least 3 sites, found " << physical_positions.size() << endl;
    exit(1);
  } else {
    cout << "Found " << physical_positions.size() << " sites.\n";
  }
  if (mutation_rate <= 0) {
    cerr << "Need a strictly positive mutation rate.\n";
    exit(1);
  }
  num_sites = physical_positions.size();
  num_samples = 0;

  for (int i = 0; i < num_sites - 1; i ++) {
    if (physical_positions[i + 1] <= physical_positions[i]) {
      cerr << "Physical positions must be strictly increasing, found ";
      cerr << physical_positions[i + 1] << " after " << physical_positions[i] << endl;
      exit(1);
    }
    if (genetic_positions[i + 1] <= genetic_positions[i]) {
      cerr << "Genetic coordinates must be strictly increasing, found ";
      cerr << genetic_positions[i + 1] << " after " << genetic_positions[i] << endl;
      exit(1);
    }
  }

  // Demography:
  cout << demography << endl;

  // Initialize both ends of the linked-list columns
  for (int i = 0; i < num_sites + 1; i++) {
    tops.emplace_back(std::make_unique<Node>(-1, 0));
    bottoms.emplace_back(std::make_unique<Node>(-1, 1));
    Node* top_i = tops.back().get();
    Node* bottom_i = bottoms.back().get();

    top_i->below = bottom_i;
    bottom_i->above = top_i;

    if (i > 0) {
      // tops[i - 1]->w[0] = {top_i, top_i};
      tops[i - 1]->w[0] = top_i;
      tops[i - 1]->w[1] = top_i;
      // bottoms[i - 1]->w = {bottom_i, bottom_i};
      bottoms[i - 1]->w[0] = bottom_i;
      bottoms[i - 1]->w[1] = bottom_i;
    }
  }

  std::tie(bp_boundaries, bp_sizes) = site_sizes(physical_positions);
  std::tie(cm_boundaries, cm_sizes) = site_sizes(genetic_positions);
}

std::tuple<std::vector<double>, std::vector<double>> DPBWT::site_sizes(std::vector<double> positions) {
  // Find mid-points between sites
  std::vector<double> pos_means(num_sites - 1);
  for (int i = 0; i < num_sites - 1; i++) {
    pos_means[i] = (positions[i] + positions[i + 1]) / 2.;
  }
  // Find the mean size of mid-point differences
  std::vector<double> site_sizes(num_sites);
  // Mid-point deltas tell us about the area around each site 
  for (int i = 1; i < num_sites - 1; i++) {
    site_sizes[i] = (pos_means[i] - pos_means[i - 1]);
  }
  double mean_size = (pos_means[num_sites - 2] - pos_means[0]) / double(num_sites - 2);
  site_sizes[0] = mean_size;
  site_sizes[num_sites - 1] = mean_size;
  for (double s : site_sizes) {
    if (s < 0) {
      cerr << "Found negative site size " << s << endl;
      exit(1);
    }
  }
  std::vector<double> boundaries(num_sites + 1);
  boundaries[0] = positions[0];
  boundaries[num_sites] = positions[num_sites - 1];
  for (int i = 1; i < num_sites; i++) {
    boundaries[i] = pos_means[i - 1];
  }
  return std::tuple(boundaries, site_sizes);
}

/**
 * @brief Find the next insert position
 * 
 * @param t Pointer to node below the sequence being inserted at site i 
 * @param g Allele at site i+1
 * @param i The site 
 * @return Node* Pointer to node below the sequence being inserted at site i+1
 */
Node* DPBWT::extend_node(Node* t, bool g, int i) {
  Node* t_next;
  if (!g && t->w[g]->sample_ID == -1) {
    // Logic is that if we're at the bottom and have a "0" genotype we jump up to last
    // 0-spot in next column
    t_next = tops[i]->below->w[1];
    // t_next = tops[t->site]->below->w[1];
  } else {
    t_next = t->w[g];
  }
  return t_next;
}

void DPBWT::insert(std::vector<bool> genotype) {
  insert(num_samples, genotype);
}
/**
 * @brief Insert a new sequence into the dynamic panel
 * 
 * @param genotype 
 */
void DPBWT::insert(int ID, std::vector<bool> genotype) {

  if (ID_map.find(ID) != ID_map.end()) {
    cerr << "ID " << ID << " is already in the panel.\n";
    exit(1);
  }
  if (genotype.size() != num_sites) {
    cerr << "Number of input markers does not match map.\n";
    exit(1);
  }
  if ((num_samples + 1) % 100 == 0) {
    cout << "Inserting haplotype number " << num_samples + 1 << endl;
  }
  int insert_index = num_samples;
  ID_map[ID] = insert_index;

  panel.emplace_back(std::vector<std::unique_ptr<Node>>(num_sites + 1));

  Node* t0 = bottoms[0].get();
  panel[ID_map[ID]][0] = std::move(std::make_unique<Node>(ID, genotype[0]));
  Node* z0 = panel[ID_map[ID]][0].get();

  // Inserts z0 above t0
  t0->insert_above(z0);

  Node* z_k = z0;
  Node* z_next;
  Node* tmp;
  Node* t_k;
  Node* t_next;
  // Insert new sequence into panel
  for (int k = 0; k < num_sites; k++) {
    bool g_k = genotype[k];
    bool next_genotype = (k == num_sites - 1) ? END_ALLELE : genotype[k + 1];
    // Add current thingy to panel
    // panel[ID_map[ID]].emplace_back(std::make_unique<Node>(ID, next_genotype));
    panel[ID_map[ID]][k + 1] = std::move(std::make_unique<Node>(ID, next_genotype));
    z_next = panel[ID_map[ID]][k + 1].get(); //.back().get();
    tmp = z_k->above;
    while (tmp->sample_ID != -1 && tmp->genotype != g_k) {
      tmp->w[g_k] = z_next;
      tmp = tmp->above;
    }

    // Figure out where to insert next node
    z_k->w[g_k] = z_next;
    t_k = z_k->below;
    z_k->w[!g_k] = t_k->w[!g_k];

    t_next = extend_node(t_k, g_k, k);
    z_next->sample_ID = ID;
    t_next->insert_above(z_next);
    z_k = z_next;
  }

  num_samples++;
  // Compute divergence values.
  // int z_dtmp = num_sites;
  // int b_dtmp = num_sites;
  // for (int k = num_sites; k >= 0; k--) {
  //   z_dtmp = std::min(z_dtmp, k);
  //   b_dtmp = std::min(b_dtmp, k);
  //   z_k = panel[current_ID][k].get();
  //   int above_ID = z_k->above->sample_ID;
  //   int below_ID = z_k->below->sample_ID;
  //   if (above_ID != -1) {
  //     while (z_dtmp > 0 && genotype[z_dtmp - 1] == panel[above_ID][z_dtmp - 1]->genotype) {
  //       z_dtmp--;
  //     }
  //   }
  //   // Update divergence value for z_k
  //   z_k->divergence = z_dtmp;
  //   if (k > 0 && below_ID != -1) {
  //     while (b_dtmp > 0 && genotype[b_dtmp - 1] == panel[below_ID][b_dtmp - 1]->genotype) {
  //       b_dtmp--;
  //     }
  //     // Only set this for real nodes
  //     z_k->below->divergence = b_dtmp;
  //   }
  // }

}

/**
 * @brief Deletes sequence ID from the dynamic panel. This moves the last sequence in the panel
 * to the position ID held. See alg 5 from d-PBWT paper.
 * 
 * @param ID 
 */
void DPBWT::delete_ID(int ID) {
  Node* s = panel[ID_map[ID]][0].get();
  // The last sequence in the panel
  int last_ID = panel[num_samples - 1][0]->sample_ID;

  for (int k = 0; k < num_sites; k++) {
    // Get allele
    bool s_allele = panel[ID_map[ID]][k]->genotype;

    // Update w-values
    Node* tmp = s->above;
    while (tmp != nullptr && tmp->genotype != s->genotype) {
      tmp->w[s_allele] = s->below->w[s_allele];
      tmp = tmp->above;
    }

    // Unlink node
    Node* tmp_above = s->above;
    Node* tmp_below = s->below;
    tmp_above->below = tmp_below;
    tmp_below->above = tmp_above;

    // Find next node in sequence
    s = extend_node(s, s_allele, k);
  }
  // Unlink node
  Node* tmp_above = s->above;
  Node* tmp_below = s->below;
  tmp_above->below = tmp_below;
  tmp_below->above = tmp_above;

  // If needed, move last sequence to replace the one we just deleted.
  if (ID_map[ID] != num_samples - 1) {
    for (int k = 0; k <= num_sites; k++) {
      // Hope this is correct
      panel[ID_map[ID]][k] = std::move(panel[num_samples - 1][k]);
    }
    ID_map[last_ID] = ID_map[ID];
    ID_map.erase(ID);
  }
  panel.pop_back();
  num_samples--;
}

/**
 * @brief For debugging: print the sample-IDs of the arrayified panel.
 * 
 */
void DPBWT::print_sorting() {
  for (int j = 0; j < num_sites + 1; ++j) {
    Node* node = tops[j].get();
    while (node != nullptr) {
      cout << node->sample_ID << " ";
      node = node->below;
    }
    cout << endl;
  }
}

/**
 * @brief Run Li-Stephens on input haplotype *without* inserting into the dynamic panel.
 *        See also Algorithm 4 of Lunter (2018), Bioinformatics.
 * 
 * @param genotype 
 * @return std::vector<std::tuple<int, int>> a pair containing path segments (start_pos, id) 
 */
std::vector<std::tuple<int, int>> DPBWT::fastLS(std::vector<bool> genotype) {
  // Get mutation/recombination penalties;
  std::vector<double> mu;
  std::vector<double> mu_c;
  std::tie(mu, mu_c) = mutation_penalties_correct();// correct_penalties ? mutation_penalties_correct() : mutation_penalties();
  // NB these still rely on padding around sites (like mutations), rather than distance between sites
  std::vector<double> rho;
  std::vector<double> rho_c; 
  std::tie(rho, rho_c) = sparse_sites ? recombination_penalties() : recombination_penalties_correct();

  std::vector<std::unique_ptr<TracebackState>> traceback_states;
  traceback_states.emplace_back(std::make_unique<TracebackState>(0, -1, nullptr));
  std::vector<State> current_states;
  current_states.emplace_back(bottoms[0].get(), 0, traceback_states.back().get());

  // Just like with insertion, we start at the bottom of the first column.
  // gm is the best score for the current level of the stack
  double gm = 0;
  // gm_next holds temporary values for updates of gm
  double gm_next;
  int max_states = 0;

  bool allele;
  for (int i = 0; i < num_sites; i++) {
    bool allele = genotype[i];
    int n_states = current_states.size();
    max_states = std::max(n_states, max_states);
    if (n_states == 0) {
      cerr << "No states left on stack, something is messed up in the algorithm.\n";
      exit(1);
    }

    // cout << "\nDoing site " << i << " with genotype " << allele << " and " << current_states.size() << " states on the stack. ";

    // We can get wonky losses on the very fist sequences. This quickly evens out.
    double recomb_score = std::max(gm + mu_c[i] + rho[i], gm + mu_c[i] + rho_c[i]);
    double mut_score = std::max(gm + mu[i] + rho_c[i], gm + mu_c[i] + rho_c[i]);
    double rho_delta = std::abs(rho[i] - rho_c[i]);
    gm_next = mut_score;

    std::vector<State> new_states;
    for (State s : current_states) {
      // cout << s << endl;

      // If the current sequence is worse than taking the best previous sequence (with score gm)
      // and recombining, we don't bother.
      if (s.score + rho_c[i] + mu_c[i] < gm_next + rho_delta) {
        Node* t_next = extend_node(s.below, allele, i);
        if (extensible_by(s, t_next, allele, i)) {
          new_states.emplace_back(t_next, s.score + rho_c[i] + mu_c[i], s.traceback);
          gm_next = std::min(gm_next, s.score + rho_c[i] + mu_c[i]);
        }
      }

      // If the current sequence with mismatch is worse than taking another sequence
      // (with score gm_next) and recombining, we don't bother.
      if (s.score + rho_c[i] + mu[i] < gm_next + rho_delta) {
        Node* t_next = extend_node(s.below, !allele, i);
        if (extensible_by(s, t_next, !allele, i)) {
          new_states.emplace_back(t_next, s.score + mu[i] + rho_c[i], s.traceback);
          gm_next = std::min(gm_next, s.score + mu[i] + rho_c[i]);
          // cout << "Mutated to " << new_states.back() << endl;
        }
      }
    }

    if (new_states.size() == 0) {
      cerr << "The algorithm is in an illegal state because no new_states were created.";
      exit(1);
    }

    // Find a best state in the current layer and recombine.
    State best_extension = *(std::min_element(new_states.begin(), new_states.end(),
      [](const auto& s1, const auto& s2) { return s1.score < s2.score; }));

    if (best_extension.score < gm_next - 0.001 || best_extension.score > gm_next + 0.001) {
      cerr << "The algorithm is in an illegal state because gm_next != best_extension.score, found ";
      cerr << "best_extension.score=" << best_extension.score << " and gm_next=" << gm_next << endl;
      exit(1);
    }

    // Add the new recombinant state to the stack (we never enter this clause on the first iteration)
    traceback_states.emplace_back(std::make_unique<TracebackState>(
      i + 1, best_extension.below->above->sample_ID, best_extension.traceback));
    new_states.emplace_back(bottoms[i + 1].get(), gm_next - rho_c[i] + rho[i], traceback_states.back().get());
    gm_next = std::min(gm_next, gm_next - rho_c[i] + rho[i]);


    // Need to consider how best to use this threshold, 100 might be too low
    // if (new_states.size() > 100) {
    // if (n_prune >= 0 && new_states.size() >= n_prune) {
    if (i % 100 == 0 && n_prune >= 0 && new_states.size() >= n_prune) {
      int old_size = new_states.size();
      StateTree tree = StateTree(new_states);
      tree.prune();
      current_states = tree.dump();
      // cout << "Site " << i << " pruning: " << current_states.size() << " states left out of " << old_size << endl;
    } else {
      current_states = new_states;
    }
    gm = gm_next;
    new_states.clear();
  }

  // Trace back best final state
  // cout << current_states.size() << " states in the end.\n";
  State min_state = *(std::min_element(current_states.begin(), current_states.end(),
        [](const auto& s1, const auto& s2) { return s1.score < s2.score; }));
  
  if ((num_samples + 1) % 100 == 0) {
    cout << "Found best path with score " << min_state.score << " for sequence " << num_samples + 1;
    cout << ", using a maximum of " << max_states << " states.\n";
  }
  
  std::vector<std::tuple<int, int>> best_path;
  TracebackState* traceback = min_state.traceback;
  int match_id = min_state.below->above->sample_ID;
  while (traceback != nullptr) {
    int segment_start = traceback->site;
    best_path.emplace_back(segment_start, match_id);
    match_id = traceback->best_prev_ID;
    traceback = traceback->prev;
  }

  std::reverse(best_path.begin(), best_path.end());

  // int mismatches = 0;
  // for (int k = 0; k < best_path.size(); k++) {
  //   if (panel[best_path[k]][k]->genotype != genotype[k]) {
  //     mismatches++;
  //   }
  // }
  // cout << "Best path has " << mismatches << " mismatches.\n";
  return best_path;
}


/**
 * @brief Determine whether state can be extended through panel by appending g.
 *        May alter 
 * @param s State at site i
 * @param t_next Node at site i+1 
 * @param g Candidate genotype for s at i+1
 * @param i The site index
 * @return true
 * @return false
 */
bool DPBWT::extensible_by(State& s, const Node* t_next, const bool g, const int i) {
  int next_above_candidate = t_next->above->sample_ID;

  if (next_above_candidate == -1 || g != panel[ID_map[next_above_candidate]][i]->genotype) {
    // This case take care of non-segregating sites with the opposite allele in the panel 
    return false;
  } else if (s.below->above->sample_ID == next_above_candidate) {
    // In this case we're just extending the same sequence again
    return true;
  } else if (genotype_interval_match(s.below->above->sample_ID, next_above_candidate, s.traceback->site, i)) {
    return true;
  }
  return false;
}

/**
 * @brief This relies on panel size, mutation rate and physical map
 * 
 * @return double 
 */
std::tuple<std::vector<double>, std::vector<double>> DPBWT::mutation_penalties() {
  std::vector<double> mu(num_sites);
  std::vector<double> mu_c(num_sites);

  // The expected branch length
  // double hap_ne = 2. * ne; //this may be wrong
  // double t = 2. * hap_ne / double(num_samples); 
  double t = num_samples == 1 ? demography.std_to_gen(1. / double(num_samples)) : demography.std_to_gen(2. / double(num_samples)) ;

  for (int i = 0; i < num_sites; i++) {
    // l is in bp units
    double k = 2. * mutation_rate * bp_sizes[i] * t;
    mu_c[i] = k;
    mu[i] = -std::log1p(-std::exp(-k));
  }
  return std::tuple(mu, mu_c);
}

std::tuple<std::vector<double>, std::vector<double>> DPBWT::mutation_penalties_correct() {
  std::vector<double> mu(num_sites);
  std::vector<double> mu_c(num_sites);

  // The expected branch length
  const double t = demography.expected_branch_length(num_samples + 1);

  for (int i = 0; i < num_sites; i++) {
    // l is in bp units
    double k = 2. * mutation_rate * bp_sizes[i] * t;
    mu_c[i] = k;
    mu[i] = -std::log1p(-std::exp(-k));
  }
  return std::tuple(mu, mu_c);
}

// std::tuple<std::vector<double>, std::vector<double>> DPBWT::mutation_penalties_impute5() {
//   std::vector<double> mu(num_sites);
//   std::vector<double> mu_c(num_sites);
  
//   double hom_loss = -std::log(0.9999);
//   double het_loss = -std::log(0.0001);
//   for (int i = 0; i < num_sites; i++) {
//     mu_c[i] = hom_loss;
//     mu[i] = het_loss;
//   }
//   return std::tuple(mu, mu_c);
// }

/**
 * @brief This relies on panel size and the recombination map
 * 
 * @return double 
 */
std::tuple<std::vector<double>, std::vector<double>> DPBWT::recombination_penalties() {
  // Recall: 1cM means the expected average number of intervening 
  // chromosomal crossovers in a single generation is 0.01
  std::vector<double> rho(num_sites);
  std::vector<double> rho_c(num_sites);
  
  // The expected branch length
  double t = num_samples == 1 ? demography.std_to_gen(1. / double(num_samples)) : demography.std_to_gen(2. / double(num_samples)) ;

  for (int i = 0; i < num_sites; i++) {
    // l is in cM units
    const double k = 2. * 0.01 * cm_sizes[i] * t;
    rho_c[i] = k;
    // This is for the binary mutation model, to use {a, c, g, t} we add log(3)
    rho[i] = -std::log1p(-std::exp(-k));
  }
  return std::tuple(rho, rho_c);
}

/**
 * @brief This relies on panel size and the recombination map
 * 
 * @return double 
 */
std::tuple<std::vector<double>, std::vector<double>> DPBWT::recombination_penalties_correct() {
  // Recall: 1cM means the expected average number of intervening 
  // chromosomal crossovers in a single generation is 0.01
  std::vector<double> rho(num_sites);
  std::vector<double> rho_c(num_sites);
  
  // The expected branch length
  const double t = demography.expected_branch_length(num_samples + 1);

  for (int i = 0; i < num_sites; i++) {
    // l is in cM units
    const double k = 2. * 0.01 * cm_sizes[i] * t;
    rho_c[i] = k;
    rho[i] = -(std::log1p(-std::exp(-k)) - std::log(num_samples));
  }
  return std::tuple(rho, rho_c);
}

/**
 * @brief Date the segment based on length and n_mismatches using maximum likelihood. (No demography)
 * 
 * @param id1 
 * @param id2 
 * @param start inclusive
 * @param end exclusive
 * @return double 
 */
double DPBWT::date_segment(const int id1, const int id2, const int start, const int end) {
  if (start > end) {
    cerr << "Can't date a segment with length <= 0\n";
    exit(1);
  }
  double m = 0;
  double bp_size = 0;
  double cm_size = 0;
  for (int i = start; i < end; i++) {
    if (panel[ID_map[id1]][i]->genotype != panel[ID_map[id2]][i]->genotype) {
      m++;
    }
    bp_size += bp_sizes[i];
    cm_size += cm_sizes[i];
  }
  double mu = 2. * mutation_rate * bp_size;
  double rho = 2. * 0.01 * cm_size;
  if (sparse_sites) {
    if (m > 15) {
      cout << "Warning: very many heterozygous sites, defaulting to const-demography method.\n";
      double gamma = 1. / demography.expected_time;
      return 2. / (gamma + rho);
    }
    double numerator = 0;
    double denominator = 0;
    int K = demography.times.size();
    for (int k = 0; k < K; k++) {
      double T1 = demography.times[k];
      double gamma_k = 1. / demography.sizes[k];
      double lambda_k = gamma_k + rho;
      double coal_fac = gamma_k * std::exp(T1 * gamma_k - demography.std_times[k]);

      if (k < K - 1) {
        double T2 = demography.times[k + 1];
        numerator += coal_fac * (2. / std::pow(lambda_k, 3))
          * (boost::math::gamma_q(3, lambda_k * T1) - boost::math::gamma_q(3, lambda_k * T2));
        denominator += coal_fac * (1. / std::pow(lambda_k, 2)) 
          * (boost::math::gamma_q(2, lambda_k * T1) - boost::math::gamma_q(2, lambda_k * T2));
      } else {
        numerator += coal_fac * (2. / std::pow(lambda_k, 3)) * boost::math::gamma_q(3, lambda_k * T1);
        denominator += coal_fac * (1. / std::pow(lambda_k, 2)) * boost::math::gamma_q(2, lambda_k * T1);
      }
    }
    return numerator / denominator;
  } else {
    if (m > 15) {
      cout << "Warning: very many heterozygous sites, defaulting to const-demography method.\n";
      double gamma = 1. / demography.expected_time;
      return (m + 2) / (gamma + rho + mu);
    }
    double numerator = 0;
    double denominator = 0;
    int K = demography.times.size();
    for (int k = 0; k < K; k++) {
      double T1 = demography.times[k];
      double gamma_k = 1. / demography.sizes[k];
      double lambda_k = gamma_k + rho + mu;
      double coal_fac = gamma_k * std::exp(T1 * gamma_k - demography.std_times[k]);
      double data_fac = std::pow(mu / lambda_k, m);

      if (k < K - 1) {
        double T2 = demography.times[k + 1];
        numerator += coal_fac * data_fac * ((m + 2) / std::pow(lambda_k, 3))
          * (boost::math::gamma_q(m + 3, lambda_k * T1) - boost::math::gamma_q(m + 3, lambda_k * T2));
        denominator += coal_fac * data_fac * (1. / std::pow(lambda_k, 2)) 
          * (boost::math::gamma_q(m + 2, lambda_k * T1) - boost::math::gamma_q(m + 2, lambda_k * T2));
      } else {
        numerator += coal_fac * data_fac * ((m + 2) / std::pow(lambda_k, 3)) * boost::math::gamma_q(m + 3, lambda_k * T1);
        denominator += coal_fac * data_fac * (1. / std::pow(lambda_k, 2)) * boost::math::gamma_q(m + 2, lambda_k * T1);
      }
    }
    return numerator / denominator;
  }
}

std::tuple<std::vector<double>, std::vector<int>, std::vector<double>> DPBWT::thread(std::vector<bool> genotype) {
  return thread(num_samples, genotype);
}

std::tuple<std::vector<double>, std::vector<int>, std::vector<double>> DPBWT::thread(int new_sample_ID, std::vector<bool> genotype) {
  // Compute LS path
  std::vector<std::tuple<int, int>> best_path = fastLS(genotype);

  // Insert new genotype. NB, it's important we insert before we date, 
  // because we need the new genotype in the panel to look for het-sites
  insert(new_sample_ID, genotype);

  std::vector<double> bp_starts;
  std::vector<int> target_IDs;
  std::vector<double> segment_ages;

  // Date segments
  for (int i = 0; i < best_path.size(); i++) {
    int segment_start = std::get<0>(best_path[i]);
    int segment_end = (i == best_path.size() - 1) ? num_sites : std::get<0>(best_path[i + 1]);
    bp_starts.push_back(ceil(bp_boundaries[segment_start])); // ceil bc should be int-like
    int target_ID = std::get<1>(best_path[i]);
    target_IDs.push_back(target_ID);
    segment_ages.push_back(date_segment(new_sample_ID, target_ID, segment_start, segment_end));
    // double age_bayes = date_segment_bayes(new_sample_ID, target_ID, segment_start, segment_end);
    // threading_instructions.emplace_back(bp_start, target_ID, age_ML);
  }
  return std::tuple(bp_starts, target_IDs, segment_ages);
}

std::vector<std::tuple<std::vector<double>, std::vector<int>, std::vector<double>>> DPBWT::thread_from_file(const std::string file_path, const int n_cycle) {
  std::vector<std::tuple<std::vector<double>, std::vector<int>, std::vector<double>>> all_threads;
  std::ifstream file(file_path, std::ios_base::in | std::ios_base::binary);
  try {
    boost::iostreams::filtering_istream in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(file);
    int i = 0;
    for(std::string line; std::getline(in, line); )
    {
      std::vector<bool> genotype;
      boost::char_separator<char> sep(" ");
      boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
      BOOST_FOREACH(std::string t, tokens)
      {
        genotype.push_back(boost::lexical_cast<bool>(t));
      }
      if (genotype.size() != num_sites) {
        cerr << "Invalid genotype for haplotype " << i << ", found " << genotype.size();
        cerr << " sites, expected " << num_sites << endl;
        exit(1);
      }
      if (i == 0) {
        insert(genotype);
        all_threads.push_back(std::tuple<std::vector<double>, std::vector<int>, std::vector<double>>({}, {}, {}));
      } else {
        all_threads.push_back(thread(i, genotype));
      }
      i++;
    }
  }
  catch(const boost::iostreams::gzip_error& e) {
    std::cout << e.what() << '\n';
  }
  file.close();

  //cycling pass
  std::ifstream file2(file_path, std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_istream in2;
  in2.push(boost::iostreams::gzip_decompressor());
  in2.push(file2);
  int i = 0;
  for(std::string line; std::getline(in2, line); )
  {
    if (i >= n_cycle) {
      return all_threads;
    }
    std::vector<bool> genotype(num_sites);
    boost::char_separator<char> sep(" ");
    boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
    int m = 0;
    BOOST_FOREACH(std::string t, tokens)
    {
      genotype[m] = boost::lexical_cast<bool>(t);
      m++;
    }
    for (bool b : genotype) {
      cout << " " << static_cast<int>(b);
    }
    cout << endl;
    if (genotype.size() != num_sites) {
      cerr << "Invalid genotype for haplotype " << i << ", found " << genotype.size();
      cerr << " sites, expected " << num_sites << endl;
      exit(1);
    }

    delete_ID(i);
    all_threads.push_back(thread(i, genotype));
    i++;
  }
  return all_threads;
  // // Compute LS path
  // std::vector<std::tuple<int, int>> best_path = fastLS(genotype);
  // // Insert new genotype. NB, it's important we insert before we date, 
  // // because we need the new genotype in the panel to look for het-sites
  // insert(new_sample_ID, genotype);

  // std::vector<double> bp_starts;
  // std::vector<int> target_IDs;
  // std::vector<double> segment_ages;

  // // Date segments
  // for (int i = 0; i < best_path.size(); i++) {
  //   int segment_start = std::get<0>(best_path[i]);
  //   int segment_end = (i == best_path.size() - 1) ? num_sites : std::get<0>(best_path[i + 1]);
  //   bp_starts.push_back(ceil(bp_boundaries[segment_start])); // ceil bc should be int-like
  //   int target_ID = std::get<1>(best_path[i]);
  //   target_IDs.push_back(target_ID);
  //   segment_ages.push_back(date_segment(new_sample_ID, target_ID, segment_start, segment_end));
  //   // double age_bayes = date_segment_bayes(new_sample_ID, target_ID, segment_start, segment_end);
  //   // threading_instructions.emplace_back(bp_start, target_ID, age_ML);
  // }
  // return std::tuple(bp_starts, target_IDs, segment_ages);
}

/**
 * @brief 
 * 
 * @param id1 
 * @param id2 
 * @param start inclusive!
 * @param end exclusive!
 * @return true 
 * @return false 
 */
bool DPBWT::genotype_interval_match(const int id1, const int id2, const int start, const int end) {
  if (id1 == id2) {
    return true;
  }
  if (id1 == -1 || id2 == -1) {
    return false;
  }
  for (int i = start; i < end; i++) {
    if (panel[ID_map[id1]][i]->genotype != panel[ID_map[id2]][i]->genotype) {
      return false;
    } 
  }
  return true;
}
