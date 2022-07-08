#include "DPBWT.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>

using std::cout;
using std::endl;

int END_ALLELE = 0;

#define THROW_LINE(a) (std::string(__FILE__) + ":" + std::to_string(__LINE__) + ": " + a)


DPBWT::DPBWT(std::vector<int> _physical_positions, std::vector<float> _genetic_positions, float _mutation_rate) :
  physical_positions(_physical_positions), genetic_positions(_genetic_positions), mutation_rate(_mutation_rate)
{
  if (physical_positions.size() != genetic_positions.size()) {
    throw std::invalid_argument(THROW_LINE("Map lengths don't match."));
  }
  if (mutation_rate <= 0) {
    throw std::invalid_argument(THROW_LINE("Need a strictly positive mutation rate."));
  }

  num_sites = physical_positions.size();
  num_samples = 0;

  // Initialize both ends of the linked-list columns
  for (int i = 0; i < physical_positions.size() + 1; i++) {
    tops.emplace_back(std::make_unique<Node>(-1, i, 0));
    bottoms.emplace_back(std::make_unique<Node>(-1, i, 1));
    Node* top_i = tops.back().get();
    Node* bottom_i = bottoms.back().get();

    top_i->below = bottom_i;
    bottom_i->above = top_i;

    if (i > 0) {
      tops[i - 1]->w = {top_i, top_i};
      bottoms[i - 1]->w = {bottom_i, bottom_i};
    }
  }

  std::vector<int> positions_delta;
  for (int i = 1; i < num_sites; i++) {
    positions_delta.push_back(physical_positions[i] - physical_positions[i - 1]);
  }
}

/**
 * @brief Find the next insert position
 * 
 * @param t Pointer to node below the sequence being inserted at site i 
 * @param g Allele at site i+1
 * @return Node* Pointer to node below the sequence being inserted at site i+1
 */
Node* DPBWT::extend_node(Node* t, bool g) {
  Node* t_next;
  if (!g && t->w[g]->sample_ID == -1) {
    // Logic is that if we're at the bottom and have a "0" genotype we jump up to last
    // 0-spot in next column
    t_next = tops[t->site]->below->w[1];
  } else {
    t_next = t->w[g];
  }
  return t_next;
}

/**
 * @brief Insert a new sequence into the dynamic panel
 * 
 * @param genotype 
 */
void DPBWT::insert(std::vector<bool> genotype) {
  int current_ID = num_samples;
  if (genotype.size() != num_sites) {
    throw std::invalid_argument(THROW_LINE("Number of input markers does not match map."));
  }
  cout << "Doing sample " << current_ID << endl;
  panel.emplace_back(std::vector<std::unique_ptr<Node>>());

  Node* t0 = bottoms[0].get();
  panel[current_ID].emplace_back(std::make_unique<Node>(current_ID, 0, genotype[0]));
  Node* z0 = panel[current_ID][0].get();

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
    panel[current_ID].emplace_back(std::make_unique<Node>(current_ID, k + 1, next_genotype));
    z_next = panel[current_ID].back().get();
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
    t_next = extend_node(t_k, g_k);
    z_next->sample_ID = current_ID;
    t_next->insert_above(z_next);
    z_k = z_next;
  }

  // Compute divergence values.
  int z_dtmp = num_sites;
  int b_dtmp = num_sites;
  for (int k = num_sites; k >= 0; k--) {
    z_dtmp = std::min(z_dtmp, k);
    b_dtmp = std::min(b_dtmp, k);
    z_k = panel[current_ID][k].get();
    int above_ID = z_k->above->sample_ID;
    int below_ID = z_k->below->sample_ID;
    if (above_ID != -1) {
      while (z_dtmp > 0 && genotype[z_dtmp - 1] == panel[above_ID][z_dtmp - 1]->genotype) {
        z_dtmp--;
      }
    }
    // Update divergence value for z_k
    z_k->divergence = z_dtmp;
    if (k > 0 && below_ID != -1) {
      while (b_dtmp > 0 && genotype[b_dtmp - 1] == panel[below_ID][b_dtmp - 1]->genotype) {
        b_dtmp--;
      }
      // Only set this for real nodes
      z_k->below->divergence = b_dtmp;
    }
  }

  num_samples++;
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
 * @brief For debugging: print the divergence values of the arrayified panel.
 * 
 */
void DPBWT::print_divergence() {
  for (int j = 0; j < num_sites + 1; ++j) {
    Node* node = tops[j].get();
    while (node != nullptr) {
      cout << node->divergence << " ";
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
  if (physical_positions.size() != genotype.size()) {
    throw std::invalid_argument(THROW_LINE("Map lengths don't match."));
  }
  std::vector<int> prefix;
  int k = 0;
  Node* t = (genotype[0]) ? bottoms[1].get() : tops[0]->below->w[1];
  while (k < num_sites && t->above->sample_ID != -1 && t->above->divergence == 0) {
    k++;
    t = extend_node(t, genotype[k]);
  }
  // We only need to consider the longest set of sequences with div=0,
  // or the single first above-sequence after we first find div>0.
  int last_match = t->above->sample_ID;
  // todo: better name than k2
  int k2 = 0;
  while (panel[last_match][k2]->genotype == genotype[k2]) {
    k2++;
  }

  k = std::max(k, k2);

  for (int i = 0; i < k; i++) {
    prefix.push_back(genotype[i]);
  }
  return prefix;
}

/**
 * @brief Run Li-Stephens on input haplotype *without* inserting into the dynamic panel.
 *        See also Algorithm 4 of Lunter (2018), Bioinformatics.
 * 
 * @param genotype 
 * @return std::vector<int> 
 */
std::vector<int> DPBWT::fastLS(std::vector<bool> genotype) {
  // TODO: make these unique ptrs 
  std::vector<std::unique_ptr<TracebackState>> traceback_states;
  traceback_states.emplace_back(std::make_unique<TracebackState>(0, -1, nullptr));
  std::vector<State> current_states;
  current_states.emplace_back(bottoms[0].get(), 0, traceback_states.back().get(), -1);
  // current_states.insert(*(recombinant_states.back()));

  // Just like with insertion, we start at the bottom of the first column.
  // states.emplace_back(std::make_unique<State>(bottoms[0].get(), 0, nullptr, 0, genotype[0], true, -1));
  // current_states.insert({ -1, states.back().get() });
  // gm is the best score for the current level of the stack
  double gm = 0;
  // gm_next holds temporary values for updates of gm
  double gm_next;
  double mu = mutation_penalty();
  double rho = recombination_penalty();
  // Whether a sequence has been successfully extended at a given site 
  bool extensible;
  Node* t_next;
  // The haplotype match
  int best_above;
  // Collect states for the next iteration
  bool allele;
  // int i = 0;
  for (int i = 0; i < num_sites; i++) {
    bool allele = genotype[i];

  // }
  // while (i < num_sites - 1) {
    cout << "\nDoing site " << i << " with genotype " << allele << " and " << current_states.size() << " states on the stack.\n";
    // i++;
    gm_next = gm + mu;
    // Whether a sequence can be extended
    bool extended = false;
    // bool allele = (i == num_sites) ? END_ALLELE : genotype[i];
    std::vector<State> new_states;
    for (State s : current_states) {
      cout << s << endl;
      Node* t_i = s.below;
      double score = s.score;
      // If the current sequence is worse than taking another sequence (with score gm)
      // and recombining, we don't bother.
      if (score < gm + rho) {
        t_next = extend_node(t_i, allele);
        std::tie(extensible, best_above) = extensible_by(s, t_next, allele);
        if (extensible) {
          // TODO: Massive performance boost by only making this new state if it has the lowest score of all states
          // with state.below = t_next, otherwise moving along. This can be done using a hash-table dict thingy instead of 
          // the a vector for new_states. (I think the correct structure is a "map")
          new_states.emplace_back(t_next, score, s.traceback, best_above);
          cout << "\tExtended to " << new_states.back() << endl;
          // if (new_states.find(t_next->sample_ID) == new_states.end() || new_states[t_next->sample_ID]->score > score) {
          //   states.emplace_back(std::make_unique<State>(t_next, score, s, s.traceback, allele, div_extensible, best_above));
          //   new_states[t_next->sample_ID] = states.back().get();
          // }
          // Will probably need a no-recombination penalty, it might
          // start counting more in the array setting.
          gm_next = std::min(gm_next, score);
          if (score == gm) {
            extended = true;
          }
        }
      }

      // If the current sequence with mismatch is worse than taking another sequence
      // (with score gm_next) and recombining, we don't bother.
      if (score + mu < gm_next + rho) {
        // cout << "Mutation candidate\n";
        t_next = extend_node(t_i, !allele);
        std::tie(extensible, best_above) = extensible_by(s, t_next, !allele);
        if (extensible) {
          new_states.emplace_back(t_next, score + mu, s.traceback, best_above);
          cout << "\tMutated to " << new_states.back() << endl;
          // if (new_states.find(t_next->sample_ID) == new_states.end() || new_states[t_next->sample_ID]->score > score + mu) {
          //   states.emplace_back(std::make_unique<State>(t_next, score + mu, s, s.traceback, !allele, div_extensible, best_above));
          //   new_states[t_next->sample_ID] = states.back().get();
          // }
          // states.emplace_back(std::make_unique<State>(t_next, score + mu, s, s.traceback, !allele, div_extensible, best_above));
          // new_states.push_back(states.back().get());
        }
      }
    }

    // No state was extended, so we're forced to recombine
    if (!extended) {
      // TODO: appropriately handle non-segregating sites here.
      // cout << "Nothing extended, recombining.\n";
      // Find a best state in the previous layer to add as parent to the recombinant state.
      State best_prev = *(std::min_element(current_states.begin(), current_states.end(),
        [](const auto& s1, const auto& s2) { return s1.score < s2.score; }));
      // State* best_prev = min_pair.second;

      // Short sanity check
      if (best_prev.score != gm) {
        throw std::runtime_error(THROW_LINE("The algorithm is in an illegal state because gm != best_prev.score"));
      }

      // Add the new recombinant state to the stack
      // (we never enter this clause on the first iteration)
      traceback_states.emplace_back(std::make_unique<TracebackState>(i, best_prev.below->above->sample_ID, best_prev.traceback));
      t_next = extend_node(bottoms[i - 1].get(), allele);
      new_states.emplace_back(t_next, gm + rho, traceback_states.back().get(), -1);
      // if (new_states.find(t_next->sample_ID) == new_states.end() || new_states[t_next->sample_ID]->score > gm + rho) {
      //   states.emplace_back(std::make_unique<State>(t_next, gm + rho, best_prev, i, allele, true, -1));
      //   new_states[t_next->sample_ID] = states.back().get();
      // }

      // states.emplace_back(std::make_unique<State>(t_next, gm + rho, best_prev, i, allele, true, -1));
      // new_states.push_back(states.back().get());
    }
    // This is where the magic happens
    // Need to consider how best to use this threshold, 100 might be too low
    if (new_states.size() > 100) {
      StateTree tree = StateTree(new_states);
      tree.prune();
      current_states = tree.dump();
    } else {
      current_states = new_states;
    }

    gm = gm_next;
    new_states.clear();
  }

  // Trace back best final state
  State min_state = *(std::min_element(current_states.begin(), current_states.end(),
        [](const auto& s1, const auto& s2) { return s1.score < s2.score; }));
  cout << "Found best path with score " << min_state.score << " ending at " << min_state << endl;
  std::vector<int> best_path;

  int target = num_sites;
  TracebackState* traceback = min_state.traceback;
  int match_id = min_state.below->above->sample_ID;
  int segment_length = min_state.site() - traceback->site;
  for (int i = 0; i <= std::min(segment_length, num_sites - 1); i++) {
    best_path.push_back(match_id);
  }
  TracebackState* last_traceback = traceback;
  traceback = traceback->prev;
  while (traceback != nullptr) {
    segment_length = last_traceback->site - traceback->site;
    match_id = last_traceback->best_prev_ID;
    for (int i = 0; i < segment_length; i++) {
      best_path.push_back(match_id);
    }
    last_traceback = traceback;
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
 * @return true
 * @return false
 */
std::tuple<bool, int> DPBWT::extensible_by(State& s, const Node* t_next, const bool g) {
  bool extensible = false;
  int best_above = s.best_above;
  if (s.best_above != -1 && g == panel[s.best_above][t_next->site - 1]->genotype) {
    // In this case we already know the single best match in the region
    extensible = true;
  } else if (t_next->above->divergence <= s.traceback->site) {
    // In this case we've matched with a whole set of haplotypes in the region
    extensible = true;
  } else {
    // In this case we don't have a whole set and we don't know whether there
    // is a single best match
    Node* best_prefix_candidate = t_next->above;
    if (best_prefix_candidate->sample_ID != -1 && best_prefix_candidate->genotype == g && 
        genotype_interval_match(s.below->above->sample_ID, best_prefix_candidate->sample_ID, s.traceback->site, s.site())) {
      best_above = best_prefix_candidate->sample_ID;
      extensible = true;
    }
  }
  // if (s->extended_by_div) {
  //   // If current state was extended by looking at divergence values,
  //   //we first try to do the same with the new extension.
  //   if (t_next->above->divergence <= s->start) {
  //     next_extended_by_div = true;
  //     extensible = true;
  //   } else {
  //     // If extension-by-divergence fails, find a "best-above" candidate.
  //     next_extended_by_div = false;

  //   }
  // } else {
  //   // If current state was not extended by looking at divergence 
  //   next_extended_by_div = false;
  //   if () {
  //     extensible = true;
  //   }
  // }
  return std::tuple(extensible, best_above);
}

/**
 * @brief This relies on panel size, mutation rate and physical map
 * 
 * @return double 
 */
double DPBWT::mutation_penalty() {
  // 
  // std::vector<double>(positions.size()) mut_penalties;
  // mut_penalties[0] = 0;
  // for (int i = 1; i < positions.size(); i++) {
  //   double pos_delta = static_cast<double>(positions[i] - positions[i - 1]);
    
  // }
  return 0.6;
}

/**
 * @brief This relies on panel size and the recombination map
 * 
 * @return double 
 */
double DPBWT::recombination_penalty() {
  // TODO: Recombination penalty *has to be constant* for maths to work out.
  // But mutation penalty can vary
  // This can be workarounded by chunking genome into cM bins and not by site.
  // Will make emission probs somewhat more complicated (can have multiple or no muts per bin)
  // but is do-able and should not be that bad for performance, we'll just have to use a poisson
  // instead of a bernoulli for the emissions. And we also need to rewrite extension functions and
  // disallow recombination in between two sites in the same bin.
  
  return 1.0;
}

// /**
//  * @brief  
//  * 
//  * @param state 
//  * @param sample_ID 
//  * @return true 
//  * @return false 
//  */
// bool DPBWT::state_node_prefix_match(const State state, const int sample_ID) {
//   State* s = &state;
//   int start = state.start;
//   while (s != nullptr && s->below->site >= start) {
//     if (s->genotype != panel[sample_ID][s->below->site]->genotype) {
//       return false;
//     }
//     if (s->parent != nullptr && s->below->site == s->parent->below->site) {
//       cout << "State at site " << s->below->site << " has the same position as its parent.\n";
//       throw std::runtime_error(THROW_LINE("!!!"));
//     }
//     s = s->parent;
//   }
//   return true;
// }

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
    if (panel[id1][i]->genotype != panel[id2][i]->genotype) {
      return false;
    } 
  }
  return true;
}
