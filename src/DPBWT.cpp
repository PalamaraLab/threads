#include "DPBWT.hpp"
#include <math.h>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>

using std::cout;
using std::cerr;
using std::endl;

int END_ALLELE = 0;


DPBWT::DPBWT(std::vector<double> _physical_positions, std::vector<double> _genetic_positions, double _mutation_rate, double _ne) :
  physical_positions(_physical_positions), genetic_positions(_genetic_positions), mutation_rate(_mutation_rate), ne(_ne)
{
  if (physical_positions.size() != genetic_positions.size()) {
    cerr << "Map lengths don't match.\n";
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

  // Initialize both ends of the linked-list columns
  for (int i = 0; i < num_sites + 1; i++) {
    tops.emplace_back(std::make_unique<Node>(-1, 0));
    bottoms.emplace_back(std::make_unique<Node>(-1, 1));
    Node* top_i = tops.back().get();
    Node* bottom_i = bottoms.back().get();

    top_i->below = bottom_i;
    bottom_i->above = top_i;

    if (i > 0) {
      tops[i - 1]->w = {top_i, top_i};
      bottoms[i - 1]->w = {bottom_i, bottom_i};
    }
  }

  bp_sizes = site_sizes(physical_positions);
  cm_sizes = site_sizes(genetic_positions);
}

std::vector<double> DPBWT::site_sizes(std::vector<double> positions) {
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
  return site_sizes;
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
  insert(genotype, num_samples);
}
/**
 * @brief Insert a new sequence into the dynamic panel
 * 
 * @param genotype 
 */
void DPBWT::insert(std::vector<bool> genotype, int ID) {
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
    exit(1);
  }
  int insert_index = num_samples;
  ID_map[ID] = insert_index;

  panel.emplace_back(std::vector<std::unique_ptr<Node>>());

  Node* t0 = bottoms[0].get();
  panel[ID_map[ID]].emplace_back(std::make_unique<Node>(ID, genotype[0]));
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
    panel[ID_map[ID]].emplace_back(std::make_unique<Node>(ID, next_genotype));
    z_next = panel[ID_map[ID]].back().get();
    tmp = z_k->above;
    while (tmp->sample_ID != -1 && tmp->genotype != g_k) {
      tmp->w[g_k] = z_next;
      tmp = tmp->above;
    }

    // Figure out where to insert next node
    z_k->w[g_k] = z_next;
    t_k = z_k->below;
    z_k->w[!g_k] = t_k->w[!g_k];
    // This is my own addition b-c else everything is stuck at the bottom :-(
    t_next = extend_node(t_k, g_k, k);
    z_next->sample_ID = ID;
    t_next->insert_above(z_next);
    z_k = z_next;
  }

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

  num_samples++;
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
    while (tmp->genotype != s->genotype) {
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

  // If needed, move last sequence to replace the one we just deleted.
  if (ID_map[ID] != num_samples - 1) {
    for (int k = 0; k < num_sites; k++) {
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
 * @brief Get longest prefix matching input haplotype *without* inserting into the dynamic panel
 * 
 * @param genotype 
 * @return std::vector<int> 
 */
std::vector<int> DPBWT::longest_prefix(std::vector<bool> genotype) {
  cerr << "This is no longer implemented.\n";
  exit(1);
  // if (physical_positions.size() != genotype.size()) {
  //   cerr << "Map lengths don't match.";
  //   exit(1);
  // }
  std::vector<int> prefix;
  // int k = 0;
  // Node* t = (genotype[0]) ? bottoms[1].get() : tops[0]->below->w[1];
  // while (k < num_sites && t->above->sample_ID != -1 && t->above->divergence == 0) {
  //   k++;
  //   t = extend_node(t, genotype[k], k);
  // }
  // // We only need to consider the longest set of sequences with div=0,
  // // or the single first above-sequence after we first find div>0.
  // int last_match = t->above->sample_ID;
  // // todo: better name than k2
  // int k2 = 0;
  // while (panel[last_match][k2]->genotype == genotype[k2]) {
  //   k2++;
  // }

  // k = std::max(k, k2);

  // for (int i = 0; i < k; i++) {
  //   prefix.push_back(genotype[i]);
  // }
  return prefix;
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
  std::tie(mu, mu_c) = mutation_penalties();
  std::vector<double> rho;
  std::vector<double> rho_c; 
  std::tie(rho, rho_c) = recombination_penalties();


  std::vector<std::unique_ptr<TracebackState>> traceback_states;
  traceback_states.emplace_back(std::make_unique<TracebackState>(0, -1, nullptr));
  std::vector<State> current_states;
  current_states.emplace_back(bottoms[0].get(), 0, traceback_states.back().get());
  // current_states.emplace_back(bottoms[0].get(), 0, traceback_states.back().get(), -1);

  // Just like with insertion, we start at the bottom of the first column.
  // gm is the best score for the current level of the stack
  double gm = 0;
  // gm_next holds temporary values for updates of gm
  double gm_next;

  // double rho = recombination_penalty();
  // Whether a sequence has been successfully extended at a given site 
  bool extensible;
  Node* t_next;
  // Collect states for the next iteration
  bool allele;
  for (int i = 0; i < num_sites; i++) {
    bool allele = genotype[i];
    // cout << "\nDoing site " << i << " with genotype " << allele << " and " << current_states.size() << " states on the stack.\n";
    gm_next = gm + mu[i] + rho_c[i];
    // Whether a sequence can be extended
    bool extended = false;
    std::vector<State> new_states;
    for (State s : current_states) {
      // cout << s << endl;
      Node* t_i = s.below;
      double score = s.score;
      // If the current sequence is worse than taking the best previous sequence (with score gm)
      // and recombining, we don't bother.
      if (score + rho_c[i] + mu_c[i] < gm + rho[i] + mu_c[i]) {
        t_next = extend_node(t_i, allele, i);
        if (extensible_by(s, t_next, allele, i)) {
          new_states.emplace_back(t_next, score + rho_c[i] + mu_c[i], s.traceback);
          // cout << "Extended to " << new_states.back() << endl;
          gm_next = std::min(gm_next, score + rho_c[i] + mu_c[i]);
          if (score == gm) {
            // If the *best* state was extended, no need to recombine
            extended = true;
          }
        }
      }

      // If the current sequence with mismatch is worse than taking another sequence
      // (with score gm_next) and recombining, we don't bother.
      if (score + mu[i] + rho_c[i] < gm_next - rho_c[i] + rho[i]) {
        t_next = extend_node(t_i, !allele, i);
        if (extensible_by(s, t_next, !allele, i)) {
          new_states.emplace_back(t_next, score + mu[i] + rho_c[i], s.traceback);
        }
      }
    }

    // No state was extended, so we're forced to recombine. TODO: always recombine when possible
    if (i > 0) {
      // Short sanity check
      if (i == 0) {
        cerr << "The algorithm is in an illegal state because nothing was extended on the first site.\n";
        exit(1);
      }

      // Find a best state in the previous layer to add as parent to the recombinant state.
      State best_prev = *(std::min_element(current_states.begin(), current_states.end(),
        [](const auto& s1, const auto& s2) { return s1.score < s2.score; }));
      if (best_prev.score < gm - 0.00001 || best_prev.score > gm + 0.00001) {
        cerr << "The algorithm is in an illegal state because gm != best_prev.score, found ";
        cerr << "best_prev.score=" << best_prev.score << " and gm=" << gm << endl;
        exit(1);
      }

      // Add the new recombinant state to the stack (we never enter this clause on the first iteration)
      t_next = extend_node(bottoms[i].get(), allele, i);
      // Make sure the new genotype actually exists at the next site
      if (t_next->above->genotype == allele && t_next->above->sample_ID != -1) {
        traceback_states.emplace_back(std::make_unique<TracebackState>(
          i, best_prev.below->above->sample_ID, best_prev.traceback));
        new_states.emplace_back(t_next, gm + mu_c[i] + rho[i], traceback_states.back().get());
      }
    }

    // Need to consider how best to use this threshold, 100 might be too low
    // if (new_states.size() > 100) {
    if (false) {
      StateTree tree = StateTree(new_states);
      int old_size = new_states.size();
      tree.prune();
      current_states = tree.dump();
      int new_size = current_states.size();
      cout << "Site " << i << ": pruned " << old_size - new_size << " states out of " << old_size << endl;
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
  cout << "Found best path with score " << min_state.score << " ending at " << min_state << endl;
  
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
  int i = 0;
  // The expected branch length
  double hap_ne = 2. * ne; //this may be wrong
  double t = 2. * hap_ne / double(num_samples); 
  for (int i = 0; i < num_sites; i++) {
    // l is in bp units
    double k = 2. * mutation_rate * bp_sizes[i] * t;
    mu_c[i] = k;
    mu[i] = -std::log1p(-std::exp(-k));
  }
  return std::tuple(mu, mu_c);
}

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
  double hap_ne = 2. * ne; //this may be wrong
  double t = 2. * hap_ne / double(num_samples); 
  for (int i = 0; i < num_sites; i++) {
    // l is in cM units
    double k = 2. * 0.01 * cm_sizes[i] * t;
    rho_c[i] = k;
    rho[i] = -std::log1p(-std::exp(-k));
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
double DPBWT::age_segment_ML(const int id1, const int id2, const int start, const int end) {
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
  // The local max of p(l, m | t) w.r.t. t
  return (m + 1) / (mu + rho);
  // double gamma = 1 / Ne;
  // double mu = mutation_rate;

  // return 0.5 * (gamma + cm_size)^2 * (gamma + mu * bp_size)^(n_mismatch+1) / (gamma + mu * bp_size + cm_size) ^(m+3);
}

/**
 * @brief Date the segment based on length and n_mismatches using a demography prior and Bayes
 * 
 * @param id1 
 * @param id2 
 * @param start 
 * @param end 
 * @return double 
 */
double DPBWT::age_segment_bayes(const int id1, const int id2, const int start, const int end) {
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
  double gamma = 1. / ne;
  double mu = 2. * mutation_rate * bp_size;
  double rho = 2. * 0.01 * cm_size; 
  return (m + 2) / (gamma + rho + mu); 
  // return (m + 2) * (m + 1) * std::pow(gamma + mu, m + 1) * std::pow(gamma + rho, 2) / (gamma * std::pow(gamma + rho + mu, m + 3));
  // return 2. * std::pow((mu + gamma) / (mu + rho), n_mismatch + 2);
}

std::vector<std::tuple<int, int, double, double>> DPBWT::thread(std::vector<bool> genotype) {
  std::vector<std::tuple<int, int, double, double>> threading_instructions;
  int new_sample_ID = num_samples;

  cout << "Finding best path...\n";
  std::vector<std::tuple<int, int>> best_path = fastLS(genotype);

  cout << "Inserting...\n";
  insert(genotype);

  cout << "Dating segments...\n";
  for (int i = 0; i < best_path.size(); i++) {
    int segment_start = std::get<0>(best_path[i]);
    int segment_end = (i == best_path.size() - 1) ? num_sites : std::get<0>(best_path[i + 1]);
    int target_ID = std::get<1>(best_path[i]);
    double age_ML = age_segment_ML(new_sample_ID, target_ID, segment_start, segment_end);
    double age_bayes = age_segment_bayes(new_sample_ID, target_ID, segment_start, segment_end);
    threading_instructions.emplace_back(segment_start, target_ID, age_ML, age_bayes);
  }
  return threading_instructions;
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
