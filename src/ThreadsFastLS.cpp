#include "ThreadsFastLS.hpp"
#include <algorithm>
#include <boost/math/special_functions/gamma.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <math.h>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;

int END_ALLELE = 0;

ThreadsFastLS::ThreadsFastLS(std::vector<double> _physical_positions, std::vector<double> _genetic_positions,
                 double _mutation_rate, std::vector<double> ne, std::vector<double> ne_times,
                 bool _sparse_sites,
                 int _n_prune, // Threshold for pruning states in the tdpbwt algorithm
                 bool _use_hmm, int _burn_in_left, int _burn_in_right)
    : physical_positions(_physical_positions), genetic_positions(_genetic_positions),
      mutation_rate(_mutation_rate), demography(Demography(ne, ne_times)),
      sparse_sites(_sparse_sites), n_prune(_n_prune), use_hmm(_use_hmm),
      burn_in_left(_burn_in_left), burn_in_right(_burn_in_right) {
  if (physical_positions.size() != genetic_positions.size()) {
    cerr << "Map lengths don't match.\n";
    exit(1);
  }
  else if (physical_positions.size() <= 2) {
    cerr << "Need at least 3 sites, found " << physical_positions.size() << endl;
    exit(1);
  }
  if (mutation_rate <= 0) {
    cerr << "Need a strictly positive mutation rate.\n";
    exit(1);
  }
  num_sites = physical_positions.size();
  num_samples = 0;

  // Check maps are strictly increasing
  // for (int i = 0; i < num_sites - 1; i++) {
  //   if (physical_positions[i + 1] <= physical_positions[i]) {
  //     cerr << "Physical positions must be strictly increasing, found ";
  //     cerr << physical_positions[i + 1] << " after " << physical_positions[i] << endl;
  //     exit(1);
  //   }
  //   if (genetic_positions[i + 1] <= genetic_positions[i]) {
  //     cerr << "Genetic coordinates must be strictly increasing, found ";
  //     cerr << genetic_positions[i + 1] << " after " << genetic_positions[i] << endl;
  //     exit(1);
  //   }
  // }

  // Initialize map burn-in
  threading_start = physical_positions.front() + burn_in_left;
  // threading_end = physical_positions.back() - burn_in_right + 1;
  threading_end = physical_positions.back() - burn_in_right;
  trim_pos_start_idx = 0;
  for (int i = 0; i < num_sites; i++) {
    if (threading_start <= physical_positions[i]) {
      break;
    }
    else {
      trim_pos_start_idx++;
    }
  }
  // Keep segments that end before trim_pos_end
  trim_pos_end_idx = num_sites;
  for (int i = num_sites - 1; i >= 0; i--) {
    if (physical_positions[i] <= threading_end) {
      break;
    }
    else {
      trim_pos_end_idx--;
    }
  }
  if (trim_pos_start_idx >= trim_pos_end_idx - 3) {
    cerr << "Too few positions left after applying burn-in, need at least 3. Aborting." << endl;
    exit(1);
  }

  // Initialize both ends of the linked-list columns
  for (int i = 0; i < num_sites + 1; i++) {
    tops.emplace_back(std::make_unique<Node>(-1, i, 0));
    bottoms.emplace_back(std::make_unique<Node>(-1, i, 1));
    Node* top_i = tops.back().get();
    Node* bottom_i = bottoms.back().get();

    top_i->below = bottom_i;
    bottom_i->above = top_i;

    if (i > 0) {
      tops[i - 1]->w[0] = top_i;
      tops[i - 1]->w[1] = top_i;
      bottoms[i - 1]->w[0] = bottom_i;
      bottoms[i - 1]->w[1] = bottom_i;
    }
  }

  std::tie(bp_boundaries, bp_sizes) = site_sizes(physical_positions);
  std::tie(cm_boundaries, cm_sizes) = site_sizes(genetic_positions);

  if (use_hmm) {
    hmm = new HMM(demography, bp_sizes, cm_sizes, mutation_rate, 64);
  }
  else {
    hmm = nullptr;
  }

  // Initialize RNG
  std::random_device rd;
  rng = std::mt19937(rd());
}

std::tuple<std::vector<double>, std::vector<double>>
ThreadsFastLS::site_sizes(std::vector<double> positions) {
  // Find mid-points between sites
  int M = positions.size();
  std::vector<double> pos_means(M - 1);
  for (int i = 0; i < M - 1; i++) {
    pos_means[i] = (positions[i] + positions[i + 1]) / 2.;
  }
  // Find the mean size of mid-point differences
  std::vector<double> site_sizes(M);
  // Mid-point deltas tell us about the area around each site
  for (int i = 1; i < M - 1; i++) {
    site_sizes[i] = (pos_means[i] - pos_means[i - 1]);
  }
  double mean_size = (pos_means[M - 2] - pos_means[0]) / double(M - 2);
  site_sizes[0] = mean_size;
  site_sizes[M - 1] = mean_size;
  for (double s : site_sizes) {
    if (s < 0) {
      cerr << "Found negative site size " << s << endl;
      exit(1);
    }
  }
  std::vector<double> boundaries(M + 1);
  boundaries[0] = positions[0];
  boundaries[M] = positions[M - 1];
  for (int i = 1; i < M; i++) {
    boundaries[i] = pos_means[i - 1];
  }
  return std::tuple(boundaries, site_sizes);
}

/**
 *
 */
std::vector<double> ThreadsFastLS::trimmed_positions() {
  std::vector<double> trim_pos = {physical_positions.cbegin() + trim_pos_start_idx,
                                  physical_positions.cbegin() + trim_pos_end_idx};
  return trim_pos;
}

void ThreadsFastLS::delete_hmm() {
  if (use_hmm) {
    delete hmm;
    use_hmm = false;
  }
}

/**
 * @brief Find the next insert position
 *
 * @param t Pointer to node below the sequence being inserted at site i
 * @param g Allele at site i+1
 * @param i The site
 * @return Node* Pointer to node below the sequence being inserted at site i+1
 */
Node* ThreadsFastLS::extend_node(Node* t, bool g, int i) {
  Node* t_next;
  if (!g && t->w[g]->sample_ID == -1) {
    // Logic is that if we're at the bottom and have a "0" genotype we jump up to last
    // 0-spot in next column
    t_next = tops[i]->below->w[1];
  }
  else {
    t_next = t->w[g];
  }
  return t_next;
}

void ThreadsFastLS::insert(const std::vector<bool>& genotype) {
  insert(num_samples, genotype);
}
/**
 * @brief Insert a new sequence into the dynamic panel
 *
 * @param genotype
 */
void ThreadsFastLS::insert(const int ID, const std::vector<bool>& genotype) {

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
  panel[ID_map.at(ID)][0] = std::move(std::make_unique<Node>(ID, 0, genotype[0]));
  Node* z0 = panel[ID_map.at(ID)][0].get();

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
    panel[ID_map.at(ID)][k + 1] = std::move(std::make_unique<Node>(ID, k + 1, next_genotype));
    z_next = panel[ID_map.at(ID)][k + 1].get();
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
  int z_dtmp = num_sites;
  int b_dtmp = num_sites;
  for (int k = num_sites; k >= 0; k--) {
    z_dtmp = std::min(z_dtmp, k);
    b_dtmp = std::min(b_dtmp, k);
    z_k = panel[ID_map.at(ID)][k].get();
    int above_ID = z_k->above->sample_ID;
    int below_ID = z_k->below->sample_ID;
    if (above_ID != -1) {
      while (z_dtmp > 0 &&
             genotype[z_dtmp - 1] == panel[ID_map.at(above_ID)][z_dtmp - 1]->genotype) {
        z_dtmp--;
      }
    }
    // Update divergence value for z_k
    z_k->divergence = z_dtmp;
    if (k > 0 && below_ID != -1) {
      while (b_dtmp > 0 &&
             genotype[b_dtmp - 1] == panel[ID_map.at(below_ID)][b_dtmp - 1]->genotype) {
        b_dtmp--;
      }
      // Only set this for real nodes
      z_k->below->divergence = b_dtmp;
    }
  }
}

/**
 * @brief Deletes sequence ID from the dynamic panel. This moves the last sequence in the panel
 * to the position ID held. See alg 5 from d-PBWT paper.
 *
 * @param ID
 */
void ThreadsFastLS::remove(int ID) {
  Node* s = panel[ID_map.at(ID)][0].get();
  // The last sequence in the panel
  int last_ID = panel[num_samples - 1][0]->sample_ID;

  for (int k = 0; k < num_sites; k++) {
    // Get allele
    bool s_allele = panel[ID_map.at(ID)][k]->genotype;

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
  if (ID_map.at(ID) != num_samples - 1) {
    for (int k = 0; k <= num_sites; k++) {
      // Hope this is correct
      panel[ID_map.at(ID)][k] = std::move(panel[num_samples - 1][k]);
    }
    ID_map.at(last_ID) = ID_map.at(ID);
  }
  ID_map.erase(ID);
  panel.pop_back();
  num_samples--;
}

/**
 * @brief For debugging: print the sample-IDs of the arrayified panel.
 *
 */
void ThreadsFastLS::print_sorting() {
  for (int j = 0; j < num_sites + 1; ++j) {
    Node* node = tops[j].get();
    while (node != nullptr) {
      if (node->sample_ID >= 0) {
        cout << node->sample_ID << " ";
      }
      node = node->below;
    }
    cout << endl;
  }
}

/**
 * Run Li-Stephens on input haplotype *without* inserting into the dynamic panel.
 * See also Algorithm 4 of Lunter (2018), Bioinformatics.
 * For imputation we use the IMPUTE/Beagle mutation penalties
 */
std::pair<TracebackState*, Node*> ThreadsFastLS::fastLS(const std::vector<bool>& genotype,
                                                  bool imputation) {
  // Get mutation/recombination penalties;
  std::vector<double> mu;
  std::vector<double> mu_c;
  std::tie(mu, mu_c) = imputation ? mutation_penalties_impute5() : mutation_penalties();

  // NB these still rely on padding around sites (like mutations), rather than distance between
  // them
  std::vector<double> rho;
  std::vector<double> rho_c;
  if (imputation) {
    std::tie(rho, rho_c) = recombination_penalties_correct();
  }
  else if (sparse_sites) {
    std::tie(rho, rho_c) = recombination_penalties();
  }
  else {
    std::tie(rho, rho_c) = recombination_penalties_correct();
  }
  // traceback_states.emplace_back(std::make_unique<TracebackState>(0, -1, nullptr));
  traceback_states.emplace_back(std::make_unique<TracebackState>(0, nullptr, nullptr));
  std::vector<State> current_states;
  current_states.emplace_back(bottoms[0].get(), 0, traceback_states.back().get());

  // Just like with insertion, we start at the bottom of the first column.
  // z holds the best current score
  double z;
  int max_states = 0;
  bool allele;
  State best_extension = current_states.back();

  for (int i = 0; i < num_sites; i++) {
    bool allele = genotype[i];
    int n_states = current_states.size();
    max_states = std::max(n_states, max_states);
    if (n_states == 0) {
      cerr << "No states left on stack, something is messed up in the algorithm.\n";
      exit(1);
    }

    // Heuristically get a bound on states we want to add
    double extension_cost = rho_c[i] + mu_c[i];
    double mutation_cost = rho_c[i] + mu[i];
    double rho_delta = std::max(0.0, rho[i] - rho_c[i]);

    // Find the cheapest extension of the best previous state, if it can be
    if (extensible_by(best_extension, extend_node(best_extension.below, allele, i), allele, i)) {
      z = best_extension.score + extension_cost;
    }
    else {
      z = best_extension.score + mutation_cost;
    }

    std::vector<State> new_states;
    for (State& s : current_states) {
      // If extending the current sequence is worse than recombining on the best extension
      // we've found so far, we don't bother.
      if (s.score + extension_cost <= z + rho_delta) {
        Node* t_next = extend_node(s.below, allele, i);
        if (extensible_by(s, t_next, allele, i)) {
          new_states.emplace_back(t_next, s.score + extension_cost, s.traceback);
          z = std::min(z, s.score + extension_cost);
        }
      }

      // If mutating the current sequence is worse than recombining on the best extension
      // we've found so far, we don't bother.
      if (s.score + mutation_cost <= z + rho_delta) {
        Node* t_next = extend_node(s.below, !allele, i);
        if (extensible_by(s, t_next, !allele, i)) {
          new_states.emplace_back(t_next, s.score + mutation_cost, s.traceback);
          z = std::min(z, s.score + mutation_cost);
        }
      }
    }

    if (new_states.size() == 0) {
      cerr << "The algorithm is in an illegal state because no new_states were created.\n";
      exit(1);
    }

    // Find a best state in the current layer and recombine.
    best_extension =
        *(std::min_element(new_states.begin(), new_states.end(),
                           [](const auto& s1, const auto& s2) { return s1.score < s2.score; }));

    if (best_extension.score < z - 0.001 || best_extension.score > z + 0.001) {
      cerr << "The algorithm is in an illegal state because z != best_extension.score, found ";
      cerr << "best_extension.score=" << best_extension.score << " and z=" << z << endl;
      exit(1);
    }

    // Add the new recombinant state to the stack (we never enter this clause on the first
    // iteration)
    double recombinant_score = z - rho_c[i] + rho[i];
    traceback_states.emplace_back(std::make_unique<TracebackState>(
        i + 1, best_extension.below->above, best_extension.traceback));
    new_states.emplace_back(bottoms[i + 1].get(), recombinant_score, traceback_states.back().get());
    if (recombinant_score < z) {
      best_extension = new_states.back();
      z = recombinant_score;
    }

    // Pruning is turned off by default
    if (n_prune >= 0 && i % 100 == 0 && new_states.size() >= n_prune) {
      int old_size = new_states.size();
      StateTree tree = StateTree(new_states);
      tree.prune();
      current_states = tree.dump();
    }
    else {
      current_states = new_states;
    }
    new_states.clear();
  }

  // Trace back best final state
  State min_state =
      *(std::min_element(current_states.begin(), current_states.end(),
                         [](const auto& s1, const auto& s2) { return s1.score < s2.score; }));

  if ((num_samples + 1) % 100 == 0) {
    cout << "Found best path with score " << min_state.score << " for sequence " << num_samples + 1;
    cout << ", using a maximum of " << max_states << " states.\n";
  }

  return std::pair<TracebackState*, Node*>(min_state.traceback, min_state.below->above);
}

/**
 * Run Li-Stephens on input diplotype *without* inserting into the dynamic panel.
 * See also Algorithm 4 of Lunter (2018), Bioinformatics.
 * Warning: The code here is more verbose than it has to be
 **/
std::array<std::pair<TracebackState*, Node*>, 2>
ThreadsFastLS::fastLS_diploid(const std::vector<int>& genotype) {
  // Get mutation/recombination penalties;
  std::vector<double> mu;
  std::vector<double> mu_c;
  std::tie(mu, mu_c) = mutation_penalties();

  // NB these still rely on padding around sites (like mutations), rather than distance between
  // them, which is slightly wrong but probably not super-bad
  std::vector<double> rho;
  std::vector<double> rho_c;
  // Long story here...
  std::tie(rho, rho_c) =
      sparse_sites ? recombination_penalties() : recombination_penalties_correct();

  // Initialize traceback states
  traceback_states.emplace_back(std::make_unique<TracebackState>(0, nullptr, nullptr));
  // Initialize the list of current state-pairs
  // For haploids we have a list of states, here we have pairs
  std::vector<StatePair> current_pairs;
  current_pairs.emplace_back(bottoms[0].get(), bottoms[0].get(), 0, traceback_states.back().get(),
                             traceback_states.back().get());

  // Just like with insertion, we start at the bottom of the first column.
  // z holds the best current score for pairs and individual sequences
  double z;
  int max_state_pairs = 0;
  bool allele;
  StatePair best_pair = current_pairs.back();
  bool extensible_a0;
  bool extensible_a1;
  bool extensible_b0;
  bool extensible_b1;

  // Just like haploid, we iterate through sites
  for (int i = 0; i < num_sites; i++) {
    int allele = genotype[i];
    int n_state_pairs = current_pairs.size();
    std::vector<StatePair> new_pairs;
    max_state_pairs = std::max(n_state_pairs, max_state_pairs);
    if (n_state_pairs == 0) {
      cerr << "No state pairs left on stack, something is messed up in the algorithm.\n";
      exit(1);
    }

    // Heuristically get a bound on states we want to add,
    // this is identical to haploid
    double extension_cost = rho_c[i] + mu_c[i];
    double mutation_cost = rho_c[i] + mu[i];
    double rho_delta = std::max(0.0, rho[i] - rho_c[i]);

    // Find the cheapest extension of the best previous pair, if it can be extended
    // this helps with filtering out which extensions to process.
    // This is the "global" minimum
    std::tie(extensible_a0, extensible_b0) = pair_extensible_by(best_pair, 0, i);
    std::tie(extensible_a1, extensible_b1) = pair_extensible_by(best_pair, 1, i);
    if (allele == 0) {
      if (extensible_a0 && extensible_b0) {
        z = best_pair.score + 2 * extension_cost;
      }
      else if (extensible_a0 || extensible_b0) {
        z = best_pair.score + extension_cost + mutation_cost;
      }
      else {
        z = best_pair.score + 2 * mutation_cost;
      }
    }
    else if (allele == 1) {
      if ((extensible_a0 && extensible_b1) || (extensible_a1 && extensible_b0)) {
        z = best_pair.score + 2 * extension_cost;
      }
      else {
        z = best_pair.score + extension_cost + mutation_cost;
      }
    }
    else if (allele == 2) {
      if (extensible_a1 && extensible_b1) {
        z = best_pair.score + 2 * extension_cost;
      }
      else if (extensible_a1 || extensible_b1) {
        z = best_pair.score + extension_cost + mutation_cost;
      }
      else {
        z = best_pair.score + 2 * mutation_cost;
      }
    }
    else {
      cerr << "Only 0, 1, 2-alleles allowed." << endl;
      exit(1);
    }

    // Set local minima, this maps (anchor, traceback) to a score
    std::unordered_map<size_t, double> local_min;
    std::unordered_map<size_t, bool> extmap_0;
    std::unordered_map<size_t, bool> extmap_1;

    for (StatePair& p : current_pairs) {
      size_t key_a = pair_key(p.below_a->above->sample_ID, p.traceback_a->site);
      size_t key_b = pair_key(p.below_b->above->sample_ID, p.traceback_b->site);
      double z_pair;

      // Set/get extensibility
      if (extmap_0.count(key_a)) {
        extensible_a0 = extmap_0.at(key_a);
      }
      else {
        State tmp_a = State(p.below_a, 0., p.traceback_a);
        Node* t_next_a = extend_node(tmp_a.below, 0, i);
        extensible_a0 = extensible_by(tmp_a, t_next_a, 0, i);
        extmap_0[key_a] = extensible_a0;
      }
      if (extmap_1.count(key_a)) {
        extensible_a1 = extmap_1.at(key_a);
      }
      else {
        State tmp_a = State(p.below_a, 1., p.traceback_a);
        Node* t_next_a = extend_node(tmp_a.below, 1, i);
        extensible_a1 = extensible_by(tmp_a, t_next_a, 1, i);
        extmap_1[key_a] = extensible_a1;
      }
      if (extmap_0.count(key_b)) {
        extensible_b0 = extmap_0.at(key_b);
      }
      else {
        State tmp_b = State(p.below_b, 0., p.traceback_b);
        Node* t_next_b = extend_node(tmp_b.below, 0, i);
        extensible_b0 = extensible_by(tmp_b, t_next_b, 0, i);
        extmap_0[key_b] = extensible_b0;
      }
      if (extmap_1.count(key_b)) {
        extensible_b1 = extmap_1.at(key_b);
      }
      else {
        State tmp_b = State(p.below_b, 1., p.traceback_b);
        Node* t_next_b = extend_node(tmp_b.below, 1, i);
        extensible_b1 = extensible_by(tmp_b, t_next_b, 1, i);
        extmap_1[key_b] = extensible_b1;
      }

      // Compute scores
      if (allele == 0) {
        if (extensible_a0 && extensible_b0) {
          z_pair = p.score + 2 * extension_cost;
        }
        else if (extensible_a0 || extensible_b0) {
          z_pair = p.score + extension_cost + mutation_cost;
        }
        else {
          z_pair = p.score + 2 * mutation_cost;
        }
      }
      else if (allele == 1) {
        if ((extensible_a0 && extensible_b1) || (extensible_a1 && extensible_b0)) {
          z_pair = p.score + 2 * extension_cost;
        }
        else {
          z_pair = p.score + extension_cost + mutation_cost;
        }
      }
      else if (allele == 2) {
        if (extensible_a1 && extensible_b1) {
          z_pair = p.score + 2 * extension_cost;
        }
        else if (extensible_a1 || extensible_b1) {
          z_pair = p.score + extension_cost + mutation_cost;
        }
        else {
          z_pair = p.score + 2 * mutation_cost;
        }
      }
      if (!local_min.count(key_a)) {
        // Set local minima
        local_min[key_a] = z_pair;
      }
      else {
        local_min[key_a] = std::min(z_pair, local_min.at(key_a));
      }
      if (!local_min.count(key_b)) {
        local_min[key_b] = z_pair;
      }
      else {
        local_min[key_b] = std::min(z_pair, local_min.at(key_b));
      }
    }

    // START OF MAIN EXTENSION LOOP
    std::unordered_map<size_t, double> new_local_min;
    for (StatePair& p : current_pairs) {
      bool extended = false;
      size_t key_a = pair_key(p.below_a->above->sample_ID, p.traceback_a->site);
      size_t key_b = pair_key(p.below_b->above->sample_ID, p.traceback_b->site);
      double z_a = local_min.at(key_a);
      double z_b = local_min.at(key_b);

      Node* a_next;
      Node* b_next;
      extensible_a0 = extmap_0.at(key_a);
      extensible_a1 = extmap_1.at(key_a);
      extensible_b0 = extmap_0.at(key_b);
      extensible_b1 = extmap_1.at(key_b);

      // First case: we extend both sequences.
      double case1_cost = p.score + 2 * extension_cost;
      std::vector<Node*> added_1a;
      std::vector<Node*> added_1b;
      if (case1_cost <= std::min({z_a + rho_delta, z_b + rho_delta, z + 2 * rho_delta})) {
        if (allele == 0) {
          if (extensible_a0 && extensible_b0) {
            a_next = extend_node(p.below_a, 0, i);
            b_next = extend_node(p.below_b, 0, i);
            new_pairs.emplace_back(a_next, b_next, case1_cost, p.traceback_a, p.traceback_b);
            added_1a.push_back(a_next);
            added_1b.push_back(b_next);
            extended = true;
          }
        }
        else if (allele == 1) {
          if (extensible_a0 && extensible_b1) {
            a_next = extend_node(p.below_a, 0, i);
            b_next = extend_node(p.below_b, 1, i);
            new_pairs.emplace_back(a_next, b_next, case1_cost, p.traceback_a, p.traceback_b);
            added_1a.push_back(a_next);
            added_1b.push_back(b_next);
            extended = true;
          }
          if (extensible_a1 && extensible_b0) {
            a_next = extend_node(p.below_a, 1, i);
            b_next = extend_node(p.below_b, 0, i);
            new_pairs.emplace_back(a_next, b_next, case1_cost, p.traceback_a, p.traceback_b);
            added_1a.push_back(a_next);
            added_1b.push_back(b_next);
            extended = true;
          }
        }
        else {
          if (extensible_a1 && extensible_b1) {
            a_next = extend_node(p.below_a, 1, i);
            b_next = extend_node(p.below_b, 1, i);
            new_pairs.emplace_back(a_next, b_next, case1_cost, p.traceback_a, p.traceback_b);
            added_1a.push_back(a_next);
            added_1b.push_back(b_next);
            extended = true;
          }
        }
        // update local and global minima
        if (extended) {
          z = std::min(z, case1_cost);
          local_min[key_a] = std::min(z_a, case1_cost);
          local_min[key_b] = std::min(z_b, case1_cost);
          for (auto a1 : added_1a) {
            size_t new_key_a = pair_key(a1->above->sample_ID, p.traceback_a->site);
            if (new_local_min.count(new_key_a)) {
              new_local_min[new_key_a] = std::min(new_local_min.at(new_key_a), case1_cost);
            }
            else {
              new_local_min[new_key_a] = case1_cost;
            }
          }
          for (auto b1 : added_1b) {
            size_t new_key_b = pair_key(b1->above->sample_ID, p.traceback_b->site);
            if (new_local_min.count(new_key_b)) {
              new_local_min[new_key_b] = std::min(new_local_min.at(new_key_b), case1_cost);
            }
            else {
              new_local_min[new_key_b] = case1_cost;
            }
          }
          extended = false;
        }
      }

      // Second case: we mutate one sequence
      double case2_cost = p.score + extension_cost + mutation_cost;
      std::vector<Node*> added_2a;
      std::vector<Node*> added_2b;
      if (case2_cost <= std::min({z_a + rho_delta, z_b + rho_delta, z + 2 * rho_delta})) {
        if (allele == 0 || allele == 2) {
          if (extensible_a0 && extensible_b1) {
            Node* a_next = extend_node(p.below_a, 0, i);
            Node* b_next = extend_node(p.below_b, 1, i);
            new_pairs.emplace_back(a_next, b_next, case2_cost, p.traceback_a, p.traceback_b);
            added_2a.push_back(a_next);
            added_2b.push_back(b_next);
            extended = true;
          }
          if (extensible_a1 && extensible_b0) {
            Node* a_next = extend_node(p.below_a, 1, i);
            Node* b_next = extend_node(p.below_b, 0, i);
            new_pairs.emplace_back(a_next, b_next, case2_cost, p.traceback_a, p.traceback_b);
            added_2a.push_back(a_next);
            added_2b.push_back(b_next);
            extended = true;
          }
        }
        else {
          if (extensible_a0 && extensible_b0) {
            Node* a_next = extend_node(p.below_a, 0, i);
            Node* b_next = extend_node(p.below_b, 0, i);
            new_pairs.emplace_back(a_next, b_next, case2_cost, p.traceback_a, p.traceback_b);
            added_2a.push_back(a_next);
            added_2b.push_back(b_next);
            extended = true;
          }
          if (extensible_a1 && extensible_b1) {
            Node* a_next = extend_node(p.below_a, 1, i);
            Node* b_next = extend_node(p.below_b, 1, i);
            new_pairs.emplace_back(a_next, b_next, case2_cost, p.traceback_a, p.traceback_b);
            added_2a.push_back(a_next);
            added_2b.push_back(b_next);
            extended = true;
          }
        }
        // update local and global minima
        if (extended) {
          z = std::min(z, case2_cost);
          local_min[key_a] = std::min(local_min.at(key_a), case2_cost);
          local_min[key_b] = std::min(local_min.at(key_b), case2_cost);
          for (auto a2 : added_2a) {
            size_t new_key_a = pair_key(a2->above->sample_ID, p.traceback_a->site);
            if (new_local_min.count(new_key_a)) {
              new_local_min[new_key_a] = std::min(new_local_min.at(new_key_a), case2_cost);
            }
            else {
              new_local_min[new_key_a] = case2_cost;
            }
          }
          for (auto b2 : added_2b) {
            size_t new_key_b = pair_key(b2->above->sample_ID, p.traceback_b->site);
            if (new_local_min.count(new_key_b)) {
              new_local_min[new_key_b] = std::min(new_local_min.at(new_key_b), case2_cost);
            }
            else {
              new_local_min[new_key_b] = case2_cost;
            }
            extended = false;
          }
        }
      }

      // Third case: we mutate both sequences
      double case3_cost = p.score + 2 * mutation_cost;
      std::vector<Node*> added_3a;
      std::vector<Node*> added_3b;
      if (case3_cost <= std::min({z_a + rho_delta, z_b + rho_delta, z + 2 * rho_delta})) {
        if (allele == 0) {
          if (extensible_a1 && extensible_b1) {
            Node* a_next = extend_node(p.below_a, 1, i);
            Node* b_next = extend_node(p.below_b, 1, i);
            new_pairs.emplace_back(a_next, b_next, case3_cost, p.traceback_a, p.traceback_b);
            added_3a.push_back(a_next);
            added_3b.push_back(b_next);
            extended = true;
          }
        }
        else if (allele == 1) {
          if (extensible_a0 && extensible_b1) {
            Node* a_next = extend_node(p.below_a, 0, i);
            Node* b_next = extend_node(p.below_b, 1, i);
            new_pairs.emplace_back(a_next, b_next, case3_cost, p.traceback_a, p.traceback_b);
            added_3a.push_back(a_next);
            added_3b.push_back(b_next);
            extended = true;
          }
          if (extensible_a1 && extensible_b0) {
            Node* a_next = extend_node(p.below_a, 1, i);
            Node* b_next = extend_node(p.below_b, 0, i);
            new_pairs.emplace_back(a_next, b_next, case3_cost, p.traceback_a, p.traceback_b);
            added_3a.push_back(a_next);
            added_3b.push_back(b_next);
            extended = true;
          }
        }
        else {
          if (extensible_a0 && extensible_b0) {
            Node* a_next = extend_node(p.below_a, 0, i);
            Node* b_next = extend_node(p.below_b, 0, i);
            new_pairs.emplace_back(a_next, b_next, case3_cost, p.traceback_a, p.traceback_b);
            added_3a.push_back(a_next);
            added_3b.push_back(b_next);
            extended = true;
          }
        }
        // update local and global minima
        if (extended) {
          z = std::min(z, case3_cost);
          local_min[key_a] = std::min(local_min.at(key_a), case3_cost);
          local_min[key_b] = std::min(local_min.at(key_b), case3_cost);
          for (auto a3 : added_3a) {
            size_t new_key_a = pair_key(a3->above->sample_ID, p.traceback_a->site);
            if (new_local_min.count(new_key_a)) {
              new_local_min[new_key_a] = std::min(new_local_min.at(new_key_a), case3_cost);
            }
            else {
              new_local_min[new_key_a] = case3_cost;
            }
          }
          for (auto b3 : added_3b) {
            size_t new_key_b = pair_key(b3->above->sample_ID, p.traceback_b->site);
            if (new_local_min.count(new_key_b)) {
              new_local_min[new_key_b] = std::min(new_local_min.at(new_key_b), case3_cost);
            }
            else {
              new_local_min[new_key_b] = case3_cost;
            }
          }
          extended = false;
        }
      }
    }
    // END OF EXTENSION LOOP

    if (new_pairs.size() == 0) {
      cerr << "The algorithm is in an illegal state because no new_states were created.\n";
      exit(1);
    }

    // SINGLE RECOMBINATION EVENTS
    std::unordered_set<size_t> already_recombined;
    std::vector<StatePair> rec_pairs;

    for (StatePair& p : new_pairs) {
      size_t key_a = pair_key(p.below_a->above->sample_ID, p.traceback_a->site);
      size_t key_b = pair_key(p.below_b->above->sample_ID, p.traceback_b->site);
      double recombinant_score = p.score - rho_c[i] + rho[i];
      if (!already_recombined.count(key_a) &&
          std::abs(new_local_min.at(key_a) - p.score) < 0.0001) {
        already_recombined.insert(key_a);
        // Best a-state, so we recombine b
        traceback_states.emplace_back(
            std::make_unique<TracebackState>(i + 1, p.below_b->above, p.traceback_b));
        rec_pairs.emplace_back(p.below_a, bottoms[i + 1].get(), recombinant_score, p.traceback_a,
                               traceback_states.back().get());
      }
      if (!already_recombined.count(key_b) &&
          std::abs(new_local_min.at(key_b) - p.score) < 0.0001) {
        already_recombined.insert(key_b);
        // Best b-state, so we recombine a
        traceback_states.emplace_back(
            std::make_unique<TracebackState>(i + 1, p.below_a->above, p.traceback_a));
        rec_pairs.emplace_back(bottoms[i + 1].get(), p.below_b, recombinant_score,
                               traceback_states.back().get(), p.traceback_b);
      }
    }

    // DOUBLE RECOMBINATION EVENT
    best_pair =
        *(std::min_element(new_pairs.begin(), new_pairs.end(),
                           [](const auto& p1, const auto& p2) { return p1.score < p2.score; }));
    double double_recombinant_score = z - 2 * rho_c[i] + 2 * rho[i];
    traceback_states.emplace_back(
        std::make_unique<TracebackState>(i + 1, best_pair.below_a->above, best_pair.traceback_a));
    traceback_states.emplace_back(
        std::make_unique<TracebackState>(i + 1, best_pair.below_b->above, best_pair.traceback_b));

    // Insert single recombination states
    for (auto p : rec_pairs) {
      new_pairs.push_back(p);
    }
    // Insert double recombination state
    new_pairs.emplace_back(bottoms[i + 1].get(), bottoms[i + 1].get(), double_recombinant_score,
                           traceback_states[traceback_states.size() - 2].get(),
                           traceback_states.back().get());

    if (double_recombinant_score < z) {
      best_pair = new_pairs.back();
      z = double_recombinant_score;
    }
    if (std::abs(best_pair.score - z) > 0.0001) {
      cerr << "The algorithm is in an illegal state because z != best_pair.score, found ";
      cerr << "best_pair.score=" << best_pair.score << " and z=" << z << endl;
      exit(1);
    }
    current_pairs = new_pairs;
    new_pairs.clear();
    local_min.clear();
    new_local_min.clear();
  }

  // TRACEBACK PART
  // Find best final state
  StatePair min_pair =
      *(std::min_element(current_pairs.begin(), current_pairs.end(),
                         [](const auto& p1, const auto& p2) { return p1.score < p2.score; }));

  return {std::pair<TracebackState*, Node*>(min_pair.traceback_a, min_pair.below_a->above),
          std::pair<TracebackState*, Node*>(min_pair.traceback_b, min_pair.below_b->above)};
}

/**
 * This is a basic traceback that samples a random haplotype per from all matches per segment and
 * also stores the lowest-numbered match (for compression)
 */
std::vector<std::tuple<int, std::vector<int>>> ThreadsFastLS::traceback(TracebackState* tb, Node* match,
                                                                  bool return_all) {
  std::vector<std::tuple<int, std::vector<int>>> best_path;
  while (tb != nullptr) {
    int n_matches = 1;
    int segment_start = tb->site;
    int match_id = match->sample_ID;
    std::vector<int> div_states = {match_id};
    Node* div_node = match;
    // We also keep track of the min_id, this is useful for genotype parsing
    int min_id = match_id;
    // Find all samples that match on the segment
    while (div_node->above != nullptr && div_node->divergence <= segment_start) {
      div_node = div_node->above;
      // this is a bit awkward
      if (div_node->sample_ID == -1) {
        break;
      }
      div_states.push_back(div_node->sample_ID);
      n_matches += 1;
      if (div_node->sample_ID < min_id) {
        min_id = div_node->sample_ID;
      }
    }
    std::vector<int> sampled_states;
    if (return_all) {
      sampled_states = std::vector<int>(div_states.begin(), div_states.end());
    }
    else {
      std::uniform_int_distribution<> distrib(0, div_states.size() - 1);
      sampled_states.reserve(2);
      std::sample(div_states.begin(), div_states.end(), std::back_inserter(sampled_states), 1, rng);
      // Add the min-state as well
      sampled_states.push_back(min_id);
    }

    best_path.emplace_back(segment_start, sampled_states);
    match = tb->best_prev_node;
    tb = tb->prev;
  }
  std::reverse(best_path.begin(), best_path.end());

  return best_path;
}

/**
 * Similar to normal traceback, but picks the up-to L best matches and stores their overlap with the
 * input sequence returns a list of a tuple of lists :-P I.e., a list of segments, and each segment
 * is a tuple (sample_IDs, starts, ends) all of equal length <= L
 */
std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>>
ThreadsFastLS::traceback_impute(std::vector<bool>& genotypes, TracebackState* tb, Node* match, int L) {
  std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>> imputation_path;
  int prev_end = num_sites;
  while (tb != nullptr) {
    int n_matches = 1;
    int segment_start = tb->site;
    int segment_end = prev_end;
    prev_end = segment_start;
    int match_id = match->sample_ID;
    std::vector<int> div_states = {match_id};
    Node* div_node = match;
    // We also keep track of the min_id, this is useful for genotype parsing
    int min_id = match_id;

    // Find all samples that match on the segment
    // First traverse upwards
    while (div_node->above != nullptr && div_node->divergence <= segment_start) {
      div_node = div_node->above;
      // this is a bit awkward
      if (div_node->sample_ID == -1) {
        break;
      }
      div_states.push_back(div_node->sample_ID);
      n_matches += 1;
      if (div_node->sample_ID < min_id) {
        min_id = div_node->sample_ID;
      }
    }
    // Then traverse downwards
    div_node = match->below;
    while (div_node != nullptr && div_node->divergence <= segment_start) {
      div_states.push_back(div_node->sample_ID);
      n_matches += 1;
      if (div_node->sample_ID < min_id) {
        min_id = div_node->sample_ID;
      }
      div_node = div_node->below;
    }

    // Overlaps measured in matching sites
    std::vector<std::pair<int, int>> overlaps;
    overlaps.reserve(div_states.size());
    // NB can replace this with array
    for (const int s_id : div_states) {
      overlaps.push_back(overflow_region(genotypes, s_id, segment_start, segment_end));
    }

    // initialize original index locations
    std::vector<size_t> idx(overlaps.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    std::stable_sort(idx.begin(), idx.end(), [&overlaps](size_t i1, size_t i2) {
      return overlaps[i1].second - overlaps[i1].first < overlaps[i2].second - overlaps[i2].first;
    });

    std::vector<int> segment_starts;
    std::vector<int> sample_ids;
    std::vector<int> segment_ends;
    for (int j = idx.size() - 1; j >= std::max(0, (int) (idx.size() - L)); j--) {
      segment_starts.push_back(segment_start);
      segment_ends.push_back(segment_end);
      sample_ids.push_back(div_states[idx[j]]);
    }
    imputation_path.emplace_back(sample_ids, segment_starts, segment_ends);
    match = tb->best_prev_node;
    tb = tb->prev;
  }
  std::reverse(imputation_path.begin(), imputation_path.end());

  return imputation_path;
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
bool ThreadsFastLS::extensible_by(State& s, const Node* t_next, const bool g, const int i) {
  int next_above_candidate = t_next->above->sample_ID;

  if (next_above_candidate == -1 || g != panel[ID_map.at(next_above_candidate)][i]->genotype) {
    // This case take care of non-segregating sites with the opposite allele in the panel
    return false;
  }
  else if (s.below->above->sample_ID == next_above_candidate) {
    // In this case we're just extending the same sequence again
    return true;
  }
  else if (genotype_interval_match(
               s.below->above->sample_ID, next_above_candidate, s.traceback->site, i)) {
    return true;
  }
  return false;
}

std::pair<bool, bool> ThreadsFastLS::pair_extensible_by(StatePair& p, const bool g, const int i) {
  State tmp_a = State(p.below_a, 0., p.traceback_a);
  Node* t_next_a = extend_node(tmp_a.below, g, i);
  State tmp_b = State(p.below_b, 0., p.traceback_b);
  Node* t_next_b = extend_node(tmp_b.below, g, i);
  bool extensible_a = extensible_by(tmp_a, t_next_a, g, i);
  bool extensible_b = extensible_by(tmp_b, t_next_b, g, i);
  return std::pair<bool, bool>(extensible_a, extensible_b);
}

std::tuple<std::vector<double>, std::vector<double>> ThreadsFastLS::mutation_penalties() {
  std::vector<double> mu(num_sites);
  std::vector<double> mu_c(num_sites);

  double mean_bp_size =
      (double) (physical_positions.back() - physical_positions[0]) / (double) num_sites;

  // The expected branch length
  const double t = demography.expected_branch_length(num_samples + 1);

  for (int i = 0; i < num_sites; i++) {
    // l is in bp units
    double k = 2. * mutation_rate * mean_bp_size * t;
    mu_c[i] = k;
    mu[i] = -std::log1p(-std::exp(-k));
  }
  return std::tuple(mu, mu_c);
}

/**
 * @brief This gives the IMPUTE5 (and Beagle) recombination penalties
 *
 * @return tuple of penalty vectors
 */
std::tuple<std::vector<double>, std::vector<double>> ThreadsFastLS::mutation_penalties_impute5() {
  std::vector<double> mu(num_sites);
  std::vector<double> mu_c(num_sites);

  double hom_loss = -std::log(0.9999);
  double het_loss = -std::log(0.0001);
  for (int i = 0; i < num_sites; i++) {
    mu_c[i] = hom_loss;
    mu[i] = het_loss;
  }
  return std::tuple(mu, mu_c);
}

/**
 * @brief This gives the *sparse* recombination penalties
 *
 * @return tuple of penalty vectors
 */
std::tuple<std::vector<double>, std::vector<double>> ThreadsFastLS::recombination_penalties() {
  // Recall: 1cM means the expected average number of intervening
  // chromosomal crossovers in a single generation is 0.01
  std::vector<double> rho(num_sites);
  std::vector<double> rho_c(num_sites);

  // The expected branch length
  const double t = demography.expected_branch_length(num_samples + 1);
  for (int i = 0; i < num_sites; i++) {
    // l is in cM units
    double k = 2. * 0.01 * cm_sizes[i] * t;
    if (k == 0) {
      k = 0.000001;
    }
    rho_c[i] = k;
    // This is for the binary mutation model, to use {a, c, g, t} we add log(3)
    rho[i] = -std::log1p(-std::exp(-k));
  }
  return std::tuple(rho, rho_c);
}

/**
 * @brief This gives the correct, dense, recombination penalties
 *
 * @return tuple of penalty vectors
 */
std::tuple<std::vector<double>, std::vector<double>> ThreadsFastLS::recombination_penalties_correct() {
  // Recall: 1cM means the expected average number of intervening
  // chromosomal crossovers in a single generation is 0.01
  std::vector<double> rho(num_sites);
  std::vector<double> rho_c(num_sites);

  // The expected branch length
  const double t = demography.expected_branch_length(num_samples + 1);

  for (int i = 0; i < num_sites; i++) {
    // l is in cM units
    double k = 2. * 0.01 * cm_sizes[i] * t;
    if (k == 0) {
      k = 0.000001;
    }
    rho_c[i] = k;
    rho[i] = -(std::log1p(-std::exp(-k)) - std::log(num_samples));
  }
  return std::tuple(rho, rho_c);
}

/**
 * @brief Date the segment based on length and n_mismatches using maximum likelihood. (No
 * demography)
 *
 * @param id1
 * @param id2
 * @param start inclusive
 * @param end exclusive
 * @return double
 */
double ThreadsFastLS::date_segment(const int num_het_sites, const int start, const int end) {
  if (start > end) {
    cerr << "Can't date a segment with length <= 0\n";
    exit(1);
  }
  double m = (double) num_het_sites;
  double bp_size = 0;
  double cm_size = 0;
  for (int i = start; i < end; i++) {
    bp_size += bp_sizes[i];
    cm_size += cm_sizes[i];
  }
  double mu = 2. * mutation_rate * bp_size;
  double rho = 2. * 0.01 * cm_size;
  if (sparse_sites) {
    return ThreadsFastLS::date_segment_sparse(num_het_sites, cm_size, demography);
  }
  else {
    return ThreadsFastLS::date_segment(m, cm_size, bp_size, mutation_rate, demography);
  }
}

double ThreadsFastLS::date_segment(int num_het_sites, double cm_size, double bp_size,
                             double mutation_rate, Demography& demography) {
  int m = num_het_sites;
  double mu = 2. * mutation_rate * bp_size;
  double rho = 2. * 0.01 * cm_size;
  if (m > 15) {
    // cout << "Warning: very many heterozygous sites, defaulting to const-demography method.\n";
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
      numerator +=
          coal_fac * data_fac * ((m + 2) / std::pow(lambda_k, 3)) *
          (boost::math::gamma_q(m + 3, lambda_k * T1) - boost::math::gamma_q(m + 3, lambda_k * T2));
      denominator +=
          coal_fac * data_fac * (1. / std::pow(lambda_k, 2)) *
          (boost::math::gamma_q(m + 2, lambda_k * T1) - boost::math::gamma_q(m + 2, lambda_k * T2));
    }
    else {
      numerator += coal_fac * data_fac * ((m + 2) / std::pow(lambda_k, 3)) *
                   boost::math::gamma_q(m + 3, lambda_k * T1);
      denominator += coal_fac * data_fac * (1. / std::pow(lambda_k, 2)) *
                     boost::math::gamma_q(m + 2, lambda_k * T1);
    }
  }
  return numerator / denominator;
}

double ThreadsFastLS::date_segment_sparse(int num_het_sites, double cm_size, Demography& demography) {
  int m = num_het_sites;
  double rho = 2. * 0.01 * cm_size;
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
      numerator +=
          coal_fac * (2. / std::pow(lambda_k, 3)) *
          (boost::math::gamma_q(3, lambda_k * T1) - boost::math::gamma_q(3, lambda_k * T2));
      denominator +=
          coal_fac * (1. / std::pow(lambda_k, 2)) *
          (boost::math::gamma_q(2, lambda_k * T1) - boost::math::gamma_q(2, lambda_k * T2));
    }
    else {
      numerator += coal_fac * (2. / std::pow(lambda_k, 3)) * boost::math::gamma_q(3, lambda_k * T1);
      denominator +=
          coal_fac * (1. / std::pow(lambda_k, 2)) * boost::math::gamma_q(2, lambda_k * T1);
    }
  }
  return numerator / denominator;
}

std::tuple<std::vector<int>, std::vector<std::vector<int>>, std::vector<double>, std::vector<int>>
ThreadsFastLS::thread(const std::vector<bool>& genotype) {
  return thread(num_samples, genotype);
}

std::vector<std::tuple<int, std::vector<int>>>
ThreadsFastLS::threads_ls(const std::vector<bool>& genotype) {
  std::vector<std::tuple<int, std::vector<int>>> best_path;
  if (num_samples > 0) {
    TracebackState* tb;
    Node* match;
    std::tie(tb, match) = fastLS(genotype);
    best_path = traceback(tb, match);
    traceback_states.clear();
  }
  return best_path;
}

// todo: L is no longer required here
std::tuple<std::vector<int>, std::vector<std::vector<int>>, std::vector<double>, std::vector<int>>
ThreadsFastLS::thread(const int new_sample_ID, const std::vector<bool>& genotype) {
  // Compute LS path
  std::vector<std::tuple<int, std::vector<int>>> best_path;
  if (num_samples > 0) {
    TracebackState* tb;
    Node* match;
    std::tie(tb, match) = fastLS(genotype);
    best_path = traceback(tb, match);
    traceback_states.clear();
  }

  // Insert new genotype. NB, it's important we insert before we date,
  // because we need the new genotype in the panel to look for het-sites
  insert(new_sample_ID, genotype);

  // TODO: delete the HMM once num_samples > 100
  std::vector<int> bp_starts;
  std::vector<std::vector<int>> target_IDs;
  std::vector<double> segment_ages;
  int total_num_het_sites = 0;
  // Date segments
  for (int i = 0; i < best_path.size(); i++) {
    int segment_start = std::get<0>(best_path[i]);
    int segment_end = (i == best_path.size() - 1) ? num_sites : std::get<0>(best_path[i + 1]);
    std::vector<int> target_ID_L = std::get<1>(best_path[i]);

    std::vector<bool> het_hom_sites =
        fetch_het_hom_sites(new_sample_ID, target_ID_L[0], segment_start, segment_end);

    int num_het_sites = 0;
    for (auto s : het_hom_sites) {
      if (s) {
        num_het_sites++;
      }
    }
    total_num_het_sites += num_het_sites;
    // is it ok to have 100 here?
    if (use_hmm && num_samples < 1000) {
      // is it ok to have 10 here?
      if (num_het_sites > 5) {
        std::vector<int> breakpoints = hmm->breakpoints(het_hom_sites, segment_start);

        for (int i = 0; i < breakpoints.size(); i++) {
          int breakpoint_start = breakpoints[i];
          int breakpoint_end = (i == breakpoints.size() - 1) ? segment_end : breakpoints[i + 1];
          target_IDs.push_back(target_ID_L);
          bp_starts.push_back(static_cast<int>(ceil(bp_boundaries[breakpoint_start])));
          // This is wrong!!! need to actually re-do num het_sites!!!!
          segment_ages.push_back(date_segment(num_het_sites, breakpoint_start, breakpoint_end));
        }
      }
      else {
        target_IDs.push_back(target_ID_L);
        bp_starts.push_back(
            static_cast<int>(ceil(bp_boundaries[segment_start]))); // ceil bc should be int-like
        segment_ages.push_back(date_segment(num_het_sites, segment_start, segment_end));
      }
    }
    else {
      target_IDs.push_back(target_ID_L);
      bp_starts.push_back(
          static_cast<int>(ceil(bp_boundaries[segment_start]))); // ceil bc should be int-like
      segment_ages.push_back(date_segment(num_het_sites, segment_start, segment_end));
    }
  }
  std::vector<int> het_sites = het_sites_from_thread(new_sample_ID, bp_starts, target_IDs);
  return remove_burn_in(bp_starts, target_IDs, segment_ages, het_sites);
}

std::vector<ImputationSegment> ThreadsFastLS::impute(std::vector<bool>& genotype, int L) {
  // vector of sample_ids, seg_starts, seg_ends (buffered)
  std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>> best_path;
  Node* match;
  TracebackState* traceback;
  std::tie(traceback, match) = fastLS(genotype);
  best_path = traceback_impute(genotype, traceback, match, L);
  traceback_states.clear();

  std::vector<double> seg_ages;
  for (const std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>& segment :
       best_path) {
    // date anti-conservatively, using longest region (is this bad?)
    int seg_start = *std::min_element(std::get<1>(segment).begin(), std::get<1>(segment).end());
    int seg_end = *std::max_element(std::get<2>(segment).begin(), std::get<2>(segment).end());
    // NB num_het_sites is not used for sparse_sites, which we assume for imputation

    seg_ages.push_back(date_segment(0, seg_start, seg_end));
  }

  int num_segs = best_path.size();
  std::vector<ImputationSegment> imputation_segments;

  // Special case for the first segment
  // initialize 1/K weights
  for (int i = 0; i < num_segs; i++) {
    const std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>& segment = best_path[i];
    const std::vector<int>& samples = std::get<0>(segment);
    ImputationSegment imp_seg;
    imp_seg.seg_start = physical_positions[std::get<1>(segment)[0]];
    imp_seg.ids = samples;
    std::vector<double> wts;
    std::vector<double> ages(samples.size(), seg_ages[i]);
    imp_seg.ages = ages;
    for (auto s : samples) {
      wts.push_back(1. / samples.size());
    }
    imp_seg.weights = wts;
    imputation_segments.push_back(imp_seg);
  }
  return imputation_segments;
}

std::array<std::vector<std::tuple<int, std::vector<int>>>, 2>
ThreadsFastLS::diploid_ls(std::vector<int> unphased_genotypes) {
  std::vector<std::tuple<int, std::vector<int>>> best_path_a;
  std::vector<std::tuple<int, std::vector<int>>> best_path_b;

  auto ls_output = fastLS_diploid(unphased_genotypes);
  best_path_a = traceback(ls_output[0].first, ls_output[0].second, true);
  best_path_b = traceback(ls_output[1].first, ls_output[1].second, true);
  traceback_states.clear();
  return {best_path_a, best_path_b};
}

std::tuple<std::vector<int>, std::vector<std::vector<int>>, std::vector<double>, std::vector<int>>
ThreadsFastLS::remove_burn_in(std::vector<int>& bp_starts, std::vector<std::vector<int>>& target_IDs,
                        std::vector<double>& segment_ages, std::vector<int>& het_sites) {
  int num_segments = bp_starts.size();

  std::vector<int> trim_starts;
  std::vector<std::vector<int>> trim_IDs;
  std::vector<double> trim_ages;

  if (num_segments > 0) {
    // Keep segments that start on or before threading_start
    int seg_start_i = 0;
    for (int i = 0; i < num_segments; i++) {
      int seg_end = i == num_segments - 1 ? threading_end : bp_starts[i + 1];
      if (threading_start < seg_end) {
        break;
      }
      else {
        seg_start_i++;
      }
    }
    // Keep segments that end on or before threading_end
    int seg_end_i = num_segments;
    for (int i = num_segments - 1; i >= 0; i--) {
      if (bp_starts[i] <= threading_end) {
        break;
      }
      else {
        seg_end_i--;
      }
    }
    trim_starts = {bp_starts.begin() + seg_start_i, bp_starts.begin() + seg_end_i};
    trim_starts[0] = threading_start;
    trim_IDs = {target_IDs.begin() + seg_start_i, target_IDs.begin() + seg_end_i};
    trim_ages = {segment_ages.begin() + seg_start_i, segment_ages.begin() + seg_end_i};
  }

  std::vector<int> trim_hets;
  for (const auto site : het_sites) {
    if (threading_start <= site && site <= threading_end) {
      trim_hets.push_back(site);
    }
  }
  return std::tie(trim_starts, trim_IDs, trim_ages, trim_hets);
}

/**
 * @brief
 * @param id1
 * @param id2
 * @param start inclusive!
 * @param end exclusive!
 * @return bool whether sequences match on the interval
 */
bool ThreadsFastLS::genotype_interval_match(const int id1, const int id2, const int start,
                                      const int end) {
  if (id1 == id2) {
    return true;
  }
  if (id1 == -1 || id2 == -1) {
    return false;
  }
  for (int i = start; i < end; i++) {
    if (panel[ID_map.at(id1)][i]->genotype != panel[ID_map.at(id2)][i]->genotype) {
      return false;
    }
  }
  return true;
}

/**
 * Assuming input genotypes match sample_id on [segment_start, segment_end), how much can we extend
 * the region in either direction with out hitting a mismatch
 */
std::pair<int, int> ThreadsFastLS::overflow_region(const std::vector<bool>& genotypes,
                                             const int sample_id, const int segment_start,
                                             const int segment_end) {
  int overlap_start = segment_start;
  int overlap_end = segment_end;

  while (overlap_start > 0) {
    if (genotypes[overlap_start - 1] == panel[ID_map.at(sample_id)][overlap_start - 1]->genotype) {
      overlap_start--;
    }
    else {
      break;
    }
  }

  while (overlap_end < num_sites) {
    if (genotypes[overlap_end] == panel[ID_map.at(sample_id)][overlap_end]->genotype) {
      overlap_end++;
    }
    else {
      break;
    }
  }
  return std::pair<int, int>(overlap_start, overlap_end);
}

/**
 * Fetch het-hom status for id1 and id2 in the region specified by site indices
 */
std::vector<bool> ThreadsFastLS::fetch_het_hom_sites(const int id1, const int id2, const int start,
                                               const int end) {
  if (ID_map.find(id1) == ID_map.end()) {
    cerr << "fetch_het_hom_sites bad id1 " << id1 << endl;
    exit(1);
  }

  if (ID_map.find(id2) == ID_map.end()) {
    cerr << "fetch_het_hom_sites bad id2 " << id2 << endl;
    exit(1);
  }
  std::vector<bool> het_hom_sites(end - start);
  for (int i = start; i < end; i++) {
    bool g1 = panel[ID_map.at(id1)][i]->genotype;
    bool g2 = panel[ID_map.at(id2)][i]->genotype;
    het_hom_sites[i - start] = (g1 != g2);
  }
  return het_hom_sites;
}

// Given threading instructions, find all heterozygous sites
std::vector<int> ThreadsFastLS::het_sites_from_thread(const int focal_ID,
                                                const std::vector<int> bp_starts,
                                                const std::vector<std::vector<int>> target_IDs) {
  std::vector<int> het_sites;
  int num_segments = bp_starts.size();
  int site_i = 0;
  for (int seg_i = 0; seg_i < num_segments; seg_i++) {
    int segment_start = bp_starts[seg_i];
    int segment_end =
        seg_i == num_segments - 1 ? physical_positions.back() + 1 : bp_starts[seg_i + 1];
    int target_ID = target_IDs[seg_i][0];
    while (segment_start <= physical_positions[site_i] &&
           physical_positions[site_i] < segment_end && site_i < num_sites) {
      if (panel[ID_map.at(focal_ID)][site_i]->genotype !=
          panel[ID_map.at(target_ID)][site_i]->genotype) {
        het_sites.push_back(physical_positions[site_i]);
      }
      site_i++;
    }
  }
  if (site_i != num_sites) {
    cerr << "Found " << site_i + 1 << " sites, expected " << num_sites << endl;
    exit(1);
  }
  return het_sites;
}
