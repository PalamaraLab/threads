#include "DPBWT.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <iterator>
#include <algorithm>
using std::cout;
using std::endl;

#define THROW_LINE(a) (std::string(__FILE__) + ":" + std::to_string(__LINE__) + ": " + a)

State::State(Node* _below, double _score, State* _parent, int _start, bool _genotype, bool _extended_by_div, int _best_above) :
  below(_below), score(_score), parent(_parent), start(_start), genotype(_genotype), extended_by_div(_extended_by_div), best_above(_best_above) {
}

ostream& operator<<(ostream& os, State& state) {
  os << "State above node " << state.below->sample_ID << " at site " << state.below->site;
  os << " with score " << state.score << ", start " << state.start << ", genotype " << state.genotype;
  os << ", div status " << state.extended_by_div << ", above-seq " << state.best_above;
  if (state.parent == nullptr) {
    os << ", and no parent.";
  } else {
    os << ", and parent at site " << state.parent->below->site;
  }
  return os;
}


DPBWT::DPBWT(std::vector<int> _physical_positions, std::vector<float> _genetic_positions, float _mutation_rate) :
  physical_positions(_physical_positions), genetic_positions(_genetic_positions), mutation_rate(_mutation_rate),
  virtual_top(Node(-1, -1, 0)), virtual_bottom(Node(-1, -1, 1))
{
  if (physical_positions.size() != genetic_positions.size()) {
    throw std::invalid_argument(THROW_LINE("Map lengths don't match."));
  }
  if (mutation_rate <= 0) {
    throw std::invalid_argument(THROW_LINE("Need a strictly positive mutation rate."));
  }

  assert(mutation_rate > 0);
  num_sites = physical_positions.size();
  num_samples = 0;

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
}

Node* DPBWT::extend_node(Node* t, bool g) {
  // Node* t = node->below;
  if (!g && t->w[g]->sample_ID == -1) {
    // Logic is that if we're at the bottom and have a "0" genotype we jump up to last
    // 0-spot in next column
    t = tops[t->site]->below->w[1];
  } else {
    t = t->w[g];
  }
  return t;
}

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
    bool next_genotype = (k == num_sites - 1) ? 0 : genotype[k + 1];
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
  // Here is a comment to explain this logic:
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
  if (states.size() > 0) {
    throw std::runtime_error(THROW_LINE("Still states on stack when there shouldn't be any"));
  }
  std::vector<State*> current_states;
  states.emplace_back(std::make_unique<State>(bottoms[0].get(), 0, nullptr, 0, genotype[0], true, -1));
  current_states.push_back(states.back().get());
  int i = 0;
  // gm is the best score for the current level of the stack
  double gm = 0;
  double gm_next;
  double mu = mutation_penalty();
  double rho = recombination_penalty();
  // all will be made clear
  bool extended;
  Node* t_next;
  bool extensible;
  bool div_extensible;
  int best_above;
  std::vector<State*> new_states;
  bool allele;
  while (i < num_sites) {
    cout << "Doing site " << i << " with " << current_states.size() << " states on the stack.\n";
    i++;
    gm_next = gm + mu;
    extended = false;
    allele = (i == num_sites) ? 0 : genotype[i];
    for (State* s : current_states) {
      cout << "Doing " << *s << endl;
      Node* t_i = s->below;
      double score = s->score;
      // If the current sequence is worse than taking another sequence (with score gm)
      // and recombining, we don't bother.
      if (score < gm + rho) {
        cout << "Extension candidate\n";
        cout << "Extended " << *t_i << " by " << allele << " to " << *t_next << endl;
        t_next = extend_node(t_i, allele);
        std::tie(extensible, div_extensible, best_above) = extensible_by(*s, t_next, allele);
        if (extensible) {
          cout << "Extensible, adding a new state.\n";
          // State new_state = State(t_next, score, s, s->start, allele);
          states.emplace_back(std::make_unique<State>(t_next, score, s, s->start, allele, div_extensible, best_above));
          new_states.push_back(states.back().get());
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
        cout << "Mutation candidate\n";
        t_next = extend_node(t_i, !allele);
        std::tie(extensible, div_extensible, best_above) = extensible_by(*s, t_next, !allele);
        if (extensible) {
          cout << "Extensible, adding a new state.\n";
          states.emplace_back(std::make_unique<State>(t_next, score + mu, s, s->start, !allele, div_extensible, best_above));
          new_states.push_back(states.back().get());
        }
      }
    }
    // No state was extended, so we're forced to recombine
    if (!extended) {
      cout << "Nothing extended, recombining.\n";
      // Find a best state in the previous layer to add as parent to the recombinant state.
      State* best_prev = *(std::min_element(std::begin(current_states), std::end(current_states), 
        [](State* s1, State* s2) { return s1->score < s2->score; }));
      
      // This should hold up
      if (best_prev->score != gm) {
        throw std::runtime_error(THROW_LINE("The algorithm is in an illegal state because gm != best_prev.score"));
      }

      // Add the new recombinant state to the stack
      // (we never enter this clause on the first iteration)
      t_next = extend_node(bottoms[i - 1].get(), allele); //(allele) ? bottoms[i].get() : tops[i - 1]->below->w[1];
      states.emplace_back(std::make_unique<State>(t_next, gm + rho, best_prev, i, allele, true, -1));
      cout << "New state:\n";
      new_states.push_back(states.back().get());
      cout << *(new_states.back()) << endl;
    }
    current_states = new_states;
    gm = gm_next;
    new_states.clear();
  }
  // Trace back best final state
  State* s = *(std::min_element(std::begin(current_states), 
        std::end(current_states), 
        [](const State* s1, const State* s2) { return s1->score < s2->score; }));
  cout << "FOUND BEST PATH WITH SCORE " << s->score << endl;
  std::vector<int> best_path;

  int target = num_sites;
  int best_match = -1;
  while (s != nullptr) {
    if (s->below->site < target) {
      best_match = s->below->above->sample_ID;
      target = s->start;
    }
    best_path.push_back(best_match);
    s = s->parent;
  }
  std::reverse(best_path.begin(), best_path.end());
  states.clear();
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
std::tuple<bool, bool, int> DPBWT::extensible_by(State s, Node* t_next, bool g) {
  bool next_extended_by_div;
  bool extensible = false;
  int best_above = s.best_above;
  if (s.extended_by_div) {
    // If current state was extended by looking at divergence values, 
    //we first try to do the same with the new extension.
    if (t_next->above->divergence <= s.start) {
      next_extended_by_div = true;
      extensible = true;
    } else {
      // If extension-by-divergence fails, find a "best-above" candidate.
      next_extended_by_div = false;
      Node* best_prefix_candidate = t_next->above;
      cout << "candidate above-node: " << *(t_next->above) << endl;
      if (best_prefix_candidate->sample_ID != -1 && best_prefix_candidate->genotype == g && state_node_prefix_match(s, best_prefix_candidate->sample_ID)) {
        best_above = best_prefix_candidate->sample_ID;
        extensible = true;
      }
    }
  } else {
    // If current state was not extended by looking at divergence 
    next_extended_by_div = false;
    if (s.best_above != -1 && g == panel[s.best_above][t_next->site]->genotype) {
      extensible = true;
    }
  }
  return std::tuple(extensible, next_extended_by_div, best_above);
}

/**
 * @brief This relies on panel size, mutation rate and physical map
 * 
 * @return double 
 */
double DPBWT::mutation_penalty() {
  return 0.6;
}

/**
 * @brief This relies on panel size and the recombination map
 * 
 * @return double 
 */
double DPBWT::recombination_penalty() {
  return 1.0;
}

bool DPBWT::state_node_prefix_match(State state, int sample_ID) {
  State* s = &state;
  int start = state.start;
  while (s != nullptr && s->below->site >= start) {
    if (s->genotype != panel[sample_ID][s->below->site]->genotype) {
      cout << sample_ID << " does not have a " << s->genotype << " at " << s->below->site << endl;
      return false;
    }
    if (s->parent != nullptr && s->below->site == s->parent->below->site) {
      cout << "State at site " << s->below->site << " has the same position as its parent.\n";
      throw std::runtime_error(THROW_LINE("sorry ! "));
    }
    s = s->parent;
  }
  return true;
}
