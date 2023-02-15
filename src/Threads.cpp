#include "Threads.hpp"
#include <algorithm>
#include <boost/math/special_functions/gamma.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <math.h>
#include <memory>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <vector>
// #include <boost/iostreams/filtering_stream.hpp>
// #include <boost/iostreams/filter/gzip.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

using std::cerr;
using std::cout;
using std::endl;

int END_ALLELE = 0;

Threads::Threads(std::vector<double> _physical_positions, std::vector<double> _genetic_positions,
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
  // else {
  //   cout << "Found " << physical_positions.size() << " sites.\n";
  // }
  if (mutation_rate <= 0) {
    cerr << "Need a strictly positive mutation rate.\n";
    exit(1);
  }
  num_sites = physical_positions.size();
  num_samples = 0;

  // Check maps are strictly increasing
  for (int i = 0; i < num_sites - 1; i++) {
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

  // Set imputation state
  impute = false;
}

void Threads::set_impute(bool impute_state) {
  impute = impute_state;
}

std::tuple<std::vector<double>, std::vector<double>>
Threads::site_sizes(std::vector<double> positions) {
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
 *
 */
std::vector<double> Threads::trimmed_positions() {
  std::vector<double> trim_pos = {physical_positions.cbegin() + trim_pos_start_idx,
                                  physical_positions.cbegin() + trim_pos_end_idx};
  return trim_pos;
}

void Threads::delete_hmm() {
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
Node* Threads::extend_node(Node* t, bool g, int i) {
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

void Threads::insert(const std::vector<bool>& genotype) {
  insert(num_samples, genotype);
}
/**
 * @brief Insert a new sequence into the dynamic panel
 *
 * @param genotype
 */
void Threads::insert(const int ID, const std::vector<bool>& genotype) {

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
void Threads::remove(int ID) {
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
void Threads::print_sorting() {
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
 * @brief Run Li-Stephens on input haplotype *without* inserting into the dynamic panel.
 *        See also Algorithm 4 of Lunter (2018), Bioinformatics.
 *
 * @param genotype
 * @return std::vector<std::tuple<int, int>> a pair containing path segments (start_pos, id)
 */
// std::vector<std::tuple<int, int>> Threads::fastLS(const std::vector<bool>& genotype, const int L)
// {
std::vector<std::tuple<int, std::vector<int>>> Threads::fastLS(const std::vector<bool>& genotype,
                                                               const int L) {
  // Get mutation/recombination penalties;
  std::vector<double> mu;
  std::vector<double> mu_c;
  std::tie(mu, mu_c) = impute ? mutation_penalties_impute5()
                              : mutation_penalties_correct(); // mutation_penalties_impute5();

  // NB these still rely on padding around sites (like mutations), rather than distance between
  // them
  std::vector<double> rho;
  std::vector<double> rho_c;
  // std::tie(rho, rho_c) = recombination_penalties_correct();
  std::tie(rho, rho_c) = // recombination_penalties_correct();
      sparse_sites ? recombination_penalties() : recombination_penalties_correct();

  std::vector<std::unique_ptr<TracebackState>> traceback_states;
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
    // cout << "\nDoing site " << i << " with genotype " << allele << " and " <<
    // current_states.size() << " states on the stack. ";

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

      // If extending the current sequence is worse than recombining on the best extension
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
    // i + 1, best_extension.below->above->sample_ID, best_extension.traceback));
    new_states.emplace_back(bottoms[i + 1].get(), recombinant_score, traceback_states.back().get());
    if (recombinant_score < z) {
      best_extension = new_states.back();
      z = recombinant_score;
    }

    // Need to consider how best to use this threshold, 100 might be too high
    // Not currently using, just adds overhead
    // ...but may be useful in real data
    if (n_prune >= 0 && i % 100 == 0 && new_states.size() >= n_prune) {
      int old_size = new_states.size();
      StateTree tree = StateTree(new_states);
      tree.prune();
      current_states = tree.dump();
      // cout << "Site " << i << " pruning: " << current_states.size() << " states left out of " <<
      // old_size << endl;
    }
    else {
      current_states = new_states;
    }
    new_states.clear();
  }

  // Trace back best final state
  // cout << current_states.size() << " states in the end.\n";
  State min_state =
      *(std::min_element(current_states.begin(), current_states.end(),
                         [](const auto& s1, const auto& s2) { return s1.score < s2.score; }));

  if ((num_samples + 1) % 100 == 0) {
    cout << "Found best path with score " << min_state.score << " for sequence " << num_samples + 1;
    cout << ", using a maximum of " << max_states << " states.\n";
  }

  std::vector<std::tuple<int, std::vector<int>>> best_path;
  TracebackState* traceback = min_state.traceback;
  Node* match = min_state.below->above;
  // int match_id = min_state.below->above->sample_ID;
  while (traceback != nullptr) {
    int n_matches = 1;
    int segment_start = traceback->site;
    int match_id = match->sample_ID;
    std::vector<int> div_states = {match_id};
    Node* div_node = match;
    while (div_node->above != nullptr && div_node->divergence <= segment_start) {
      n_matches += 1;
      // cout << div_node->above->sample_ID << endl;
      div_node = div_node->above;
      div_states.push_back(div_node->sample_ID);
    }
    std::uniform_int_distribution<> distrib(0, div_states.size() - 1);
    std::vector<int> sampled_states; //(L);
    sampled_states.reserve(L);
    if (div_states.size() < L) {
      for (int i = 0; i < L; i++) {
        sampled_states.push_back(div_states[distrib(rng)]);
      }
    }
    else {
      std::sample(div_states.begin(), div_states.end(), std::back_inserter(sampled_states), L, rng);
    }

    best_path.emplace_back(segment_start, sampled_states);
    // best_path.emplace_back(segment_start, match_id);
    // match_id = traceback->best_prev_ID;
    match = traceback->best_prev_node;
    traceback = traceback->prev;
  }
  std::reverse(best_path.begin(), best_path.end());

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
bool Threads::extensible_by(State& s, const Node* t_next, const bool g, const int i) {
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

/**
 * @brief This relies on panel size, mutation rate and physical map
 *
 * @return double
 */
std::tuple<std::vector<double>, std::vector<double>> Threads::mutation_penalties() {
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

std::tuple<std::vector<double>, std::vector<double>> Threads::mutation_penalties_correct() {
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

std::tuple<std::vector<double>, std::vector<double>> Threads::mutation_penalties_impute5() {
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
 * @brief This relies on panel size and the recombination map
 *
 * @return double
 */
std::tuple<std::vector<double>, std::vector<double>> Threads::recombination_penalties() {
  // Recall: 1cM means the expected average number of intervening
  // chromosomal crossovers in a single generation is 0.01
  std::vector<double> rho(num_sites);
  std::vector<double> rho_c(num_sites);

  // The expected branch length
  const double t = demography.expected_branch_length(num_samples + 1);
  // double t = num_samples == 1 ? demography.std_to_gen(1. / double(num_samples)) :
  // demography.std_to_gen(2. / double(num_samples)) ;

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
std::tuple<std::vector<double>, std::vector<double>> Threads::recombination_penalties_correct() {
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
 * @brief Date the segment based on length and n_mismatches using maximum likelihood. (No
 * demography)
 *
 * @param id1
 * @param id2
 * @param start inclusive
 * @param end exclusive
 * @return double
 */
double Threads::date_segment(const int id1, const int id2, const int start, const int end) {
  if (start > end) {
    cerr << "Can't date a segment with length <= 0\n";
    exit(1);
  }
  // this check should become unnecessary
  if (ID_map.find(id1) == ID_map.end()) {
    cerr << "date_segment bad id1 " << id1 << endl;
    exit(1);
  }

  if (ID_map.find(id2) == ID_map.end()) {
    cerr << "date_segment bad id2 " << id2 << endl;
    exit(1);
  }
  double m = 0;
  double bp_size = 0;
  double cm_size = 0;
  for (int i = start; i < end; i++) {
    if (panel[ID_map.at(id1)][i]->genotype != panel[ID_map.at(id2)][i]->genotype) {
      m++;
    }
    bp_size += bp_sizes[i];
    cm_size += cm_sizes[i];
  }
  double mu = 2. * mutation_rate * bp_size;
  double rho = 2. * 0.01 * cm_size;
  if (sparse_sites) {
    if (m > 15) {
      // cout << "Warning: very many heterozygous sites, defaulting to const-demography method.\n";
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
        numerator +=
            coal_fac * (2. / std::pow(lambda_k, 3)) *
            (boost::math::gamma_q(3, lambda_k * T1) - boost::math::gamma_q(3, lambda_k * T2));
        denominator +=
            coal_fac * (1. / std::pow(lambda_k, 2)) *
            (boost::math::gamma_q(2, lambda_k * T1) - boost::math::gamma_q(2, lambda_k * T2));
      }
      else {
        numerator +=
            coal_fac * (2. / std::pow(lambda_k, 3)) * boost::math::gamma_q(3, lambda_k * T1);
        denominator +=
            coal_fac * (1. / std::pow(lambda_k, 2)) * boost::math::gamma_q(2, lambda_k * T1);
      }
    }
    return numerator / denominator;
  }
  else {
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
        numerator += coal_fac * data_fac * ((m + 2) / std::pow(lambda_k, 3)) *
                     (boost::math::gamma_q(m + 3, lambda_k * T1) -
                      boost::math::gamma_q(m + 3, lambda_k * T2));
        denominator += coal_fac * data_fac * (1. / std::pow(lambda_k, 2)) *
                       (boost::math::gamma_q(m + 2, lambda_k * T1) -
                        boost::math::gamma_q(m + 2, lambda_k * T2));
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
}

std::tuple<std::vector<int>, std::vector<std::vector<int>>, std::vector<double>, std::vector<int>>
Threads::thread(const std::vector<bool>& genotype, const int L) {
  return thread(num_samples, genotype, L);
}

std::tuple<std::vector<int>, std::vector<std::vector<int>>, std::vector<double>, std::vector<int>>
Threads::thread(const int new_sample_ID, const std::vector<bool>& genotype, const int L) {
  // Compute LS path
  // std::vector<std::tuple<int, int>> best_path;

  std::vector<std::tuple<int, std::vector<int>>> best_path;
  if (num_samples > 0) {
    best_path = fastLS(genotype, L);
  }

  // Insert new genotype. NB, it's important we insert before we date,
  // because we need the new genotype in the panel to look for het-sites
  insert(new_sample_ID, genotype);

  // TODO: delete the HMM once num_samples > 100
  std::vector<int> bp_starts;
  std::vector<std::vector<int>> target_IDs;
  std::vector<double> segment_ages;

  // Date segments
  for (int i = 0; i < best_path.size(); i++) {
    int segment_start = std::get<0>(best_path[i]);
    int segment_end = (i == best_path.size() - 1) ? num_sites : std::get<0>(best_path[i + 1]);
    std::vector<int> target_ID_L = std::get<1>(best_path[i]);

    if (use_hmm && num_samples <= 100) {
      std::vector<bool> het_hom_sites =
          fetch_het_hom_sites(new_sample_ID, target_ID_L[0], segment_start, segment_end);

      int num_het_sites = 0;
      for (auto s : het_hom_sites) {
        if (s) {
          num_het_sites++;
        }
      }
      if (num_het_sites > 10) {
        // cout << "found " << num_het_sites << " het sites...\n";
        std::vector<int> breakpoints = hmm->breakpoints(het_hom_sites, segment_start);
        // cout << "broke [" << segment_start << ", " << segment_end << ") up into\n";

        // cout << "\t... found " << breakpoints.size() << " breakpoints\n";
        for (int i = 0; i < breakpoints.size(); i++) {
          int breakpoint_start = breakpoints[i];
          int breakpoint_end = (i == breakpoints.size() - 1) ? segment_end : breakpoints[i + 1];
          // cout<< "[" << breakpoint_start << ", " << breakpoint_end << ") ";
          target_IDs.push_back(target_ID_L);
          bp_starts.push_back(static_cast<int>(ceil(bp_boundaries[breakpoint_start])));
          segment_ages.push_back(
              date_segment(new_sample_ID, target_ID_L[0], breakpoint_start, breakpoint_end));
        }
        // cout << "\n";
      }
      else {
        target_IDs.push_back(target_ID_L);
        bp_starts.push_back(
            static_cast<int>(ceil(bp_boundaries[segment_start]))); // ceil bc should be int-like
        segment_ages.push_back(
            date_segment(new_sample_ID, target_ID_L[0], segment_start, segment_end));
      }
    }
    else {
      target_IDs.push_back(target_ID_L);
      bp_starts.push_back(
          static_cast<int>(ceil(bp_boundaries[segment_start]))); // ceil bc should be int-like
      segment_ages.push_back(
          date_segment(new_sample_ID, target_ID_L[0], segment_start, segment_end));
    }
  }
  std::vector<int> het_sites = het_sites_from_thread(new_sample_ID, bp_starts, target_IDs);
  return remove_burn_in(bp_starts, target_IDs, segment_ages, het_sites);
  // return std::tuple(bp_starts, target_IDs, segment_ages, het_sites);
}

std::tuple<std::vector<int>, std::vector<std::vector<int>>, std::vector<double>, std::vector<int>>
Threads::remove_burn_in(std::vector<int>& bp_starts, std::vector<std::vector<int>>& target_IDs,
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
 * @return true
 * @return false
 */
bool Threads::genotype_interval_match(const int id1, const int id2, const int start,
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
 *
 */
std::vector<bool> Threads::fetch_het_hom_sites(const int id1, const int id2, const int start,
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

// This function... does... what exactly?
std::vector<int> Threads::het_sites_from_thread(const int focal_ID,
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
