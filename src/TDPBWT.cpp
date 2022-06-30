#include "TDPBWT.hpp"

// #include <pair>
// #include <sstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cassert>
#include <memory>
#include <stdexcept>
using std::cout;
using std::endl;

#define THROW_LINE(a) (std::string(__FILE__) + ":" + std::to_string(__LINE__) + ": " + a)

TDPBWT::TDPBWT(std::vector<int> _physical_positions, std::vector<float> _genetic_positions, float _mutation_rate) :
    physical_positions(_physical_positions), genetic_positions(_genetic_positions), mutation_rate(_mutation_rate),
    virtual_top(Node(-1, -1, 0)), virtual_bottom(Node(-1, -1, 1))
{
  assert(physical_positions.size() == genetic_positions.size());
  assert(mutation_rate > 0);
  num_sites = physical_positions.size();
  num_samples = 0;
}

// void TDPBWT::insert(std::vector<bool> genotype) {
//   // If nothing in here, make a pointer to all tops,
//   // then insert genotype
//   int current_ID = num_samples;
//   if (genotype.size() != num_sites) {
//     throw std::invalid_argument(THROW_LINE("Number of input markers does not match map."));
//   }
//   cout << "Doing sample " << current_ID << endl;

//   t0 = cols[0].back();
//   z0 = Node(current_ID, 0, genotype[0]);
//   // Inserts z0 above t0
//   t0->insert_above(z0);

//   z_k = z0;
//   for (int k = 0; k < N; k++) {
//     bool g_k = genotype[k];
//     Node z_next = Node(current_ID, k + 1, genotype[k + 1]);
//     Node* tmp = z_k->above;

//     while (tmp->genotype != g_k) {
//       tmp->w[g_k] = z_next;
//       tmp = tmp->above;
//     }
//     z_k->w[g_k] = z_next;
//     Node* t_k = z_k->below;
//     z_k->w[!g_k] = t_k->w[!g_k];
//     Node* t_next = t_k->w[g_k];
//     z_next->ID = current_ID;
//     t_next->insert_above(z_next);
//   }
//   num_samples++;
// }


void TDPBWT::insert(std::vector<bool> genotype) {
  // If nothing in here, make a pointer to all tops,
  // then insert genotype
  int current_ID = num_samples;
  if (genotype.size() != num_sites) {
    throw std::invalid_argument(THROW_LINE("Number of input markers does not match map."));
  }
  cout << "Doing sample " << current_ID << endl;
  // It's important new backs have genotype '1'
  all_nodes.emplace_back(std::make_unique<Node>(current_ID, num_sites, 1));
  Node* new_back = all_nodes.back().get();
  backs.push_back(new_back); // backs!

  // On empty grid, simply initialise a doubly-linked list of genotypes.
  if (num_samples == 0) {

    for (int j = 0; j < num_sites; ++j) {
      all_nodes.emplace_back(std::make_unique<Node>(current_ID, j, genotype[j]));
      Node* node = all_nodes.back().get();

      if (genotype[j]) {
        // Set prev/next 0 to null and prev/next 1 to self
        node->next_1 = node;
        node->prev_1 = node;
        node->next_0 = &virtual_bottom;
        node->prev_0 = &virtual_top;
      } else {
        // Set prev/next 1 to null and prev/next 0 to self
        node->next_0 = node;
        node->prev_0 = node;
        node->next_1 = &virtual_bottom;
        node->prev_1 = &virtual_top;
      }
      if (j > 0) {
        Node* prev_node = tops.back();
        prev_node->insert_right(node);
      }
      node->insert_right(new_back);
      tops.push_back(node);
      bottoms.push_back(node);
    }
    ++num_samples;
    print_sorting();
    return;
  }

  // Here we do the insertion bit.
  // Remember, we do it backwards!
  Node* prev_back = backs[backs.size() - 2];
  prev_back->insert_below(new_back);
  Node* prev = new_back;
  for (int j = num_sites - 1; j >= 0; --j) {
    // Node for sample current_ID at site j
    all_nodes.emplace_back(std::make_unique<Node>(current_ID, j, genotype[j]));
    Node* node = all_nodes.back().get();
    node->right = prev;
    prev->left = node;

    // Find spot to insert new node

    if (prev->genotype) {
      // If g_prev=1, find the lowest (j-1)-parent that also has x at j and map below that.
      Node* above = prev->first_uncle(genotype[j]); 
      if (above != nullptr) {
        above->insert_below(node);
        if (bottoms[j] == above) {
          bottoms[j] = node;
        }
      } else if (genotype[j]) {
        bottoms[j]->insert_below(node);
        bottoms[j] = node;
      } else {
        tops[j]->insert_above(node);
        tops[j] = node;
      }
    } else {
      // If g_prev=0, find the highest (j-1)-descendant that also has 0 at j and map above that.
      Node* below = prev->first_nephew(genotype[j]);
      if (below != nullptr) {
        below->insert_above(node);
        if (tops[j] == below) {
          tops[j] = node;
        }
      } else if (genotype[j]) {
        bottoms[j]->insert_below(node);
        bottoms[j] = node;
      } else {
        tops[j]->insert_above(node);
        tops[j] = node;
      }
    }

    // Set/update next_0s etc (U/V arrays)
    Node* old_next_1 = node->below == nullptr ? &virtual_bottom : node->below->next_1;
    Node* old_next_0 = node->below == nullptr ? &virtual_bottom : node->below->next_0;
    Node* old_prev_1 = node->above == nullptr ? &virtual_top : node->above->prev_1;
    Node* old_prev_0 = node->above == nullptr ? &virtual_top : node->above->prev_0;
    if (genotype[j]) {
      node->next_0 = old_next_0;
      node->prev_0 = old_prev_0;
      node->next_1 = node;
      node->prev_1 = node;
      if (node->above != nullptr) {
        node->above->update_next_1_recursive(old_next_1, node);
      }
      if (node->below != nullptr) {
        node->below->update_prev_1_recursive(old_prev_1, node);
      }
    } else {
      node->next_1 = old_next_1;
      node->prev_1 = old_prev_1;
      node->next_0 = node;
      node->prev_0 = node;
      if (node->above != nullptr) {
        node->above->update_next_0_recursive(old_next_0, node);
      }
      if (node->below != nullptr) {
        node->below->update_prev_0_recursive(old_prev_0, node);
      }
    }
    prev = node;
  }

  // Update divergence arrays
  // for (int j = 0; j something; ++j) {
  //   //profit
  //   // see dpbwt insertion algorithm
  // }


  ++num_samples;
  cout << "sorting:\n";
  print_sorting();
  cout << "U0:\n";
  print_U0();
  cout << "U1:\n";
  print_U1();
  return;
}

void TDPBWT::print_sorting() {
  for (int j = 0; j < num_sites; ++j) {
    Node* node = tops[j];
    while (node != nullptr) {
      cout << node->sample_ID << " ";
      node = node->below;
    }
    cout << endl;
  }
}

void TDPBWT::print_U0() {
  for (int j = 0; j < num_sites; ++j) {
    Node* node = tops[j];
    while (node != nullptr) {
      cout << node->prev_0->sample_ID << " ";
      node = node->below;
    }
    cout << endl;
  }
}

void TDPBWT::print_U1() {
  for (int j = 0; j < num_sites; ++j) {
    Node* node = tops[j];
    while (node != nullptr) {
      cout << node->prev_1->sample_ID << " ";
      node = node->below;
    }
    cout << endl;
  }
}

/*
 * Find the longest suffix match for g.
 * Returns ID of match and start position
 */
std::pair<int, int> TDPBWT::longest_suffix(std::vector<bool> genotype) {
  assert(g.size() == num_sites);
  bool have_matched = false;

  Node* ts = tops[num_sites - 1];
  Node* te = bottoms[num_sites - 1];
  Node* next_ts;
  Node* next_te;
  for (int j = num_sites - 1; j >= 0; --j) {
    cout << "Doing site " << j << endl;
    cout << "Found genotype " << genotype[j] << endl;
    if (genotype[j]) {
      next_ts = ts->next_1;
      next_te = te->prev_1;
    } else {
      next_ts = ts->next_0;
      next_te = te->prev_0;
    }
    if (next_ts != nullptr) {
      cout << *next_ts << endl;
    } else {
      cout << "nullptr\n";
    }
    if (next_te != nullptr) {
      cout << *next_te << endl;
    } else {
      cout << "nullptr\n";
    }
    // This is going to be a bit weird and I need proof that it works
    if (!have_matched) {
      have_matched = next_ts != nullptr && next_te != nullptr && next_ts->equivalent_to(next_te);
      ts = next_ts->left;
      te = next_te->left;
    } else if (next_ts->equivalent_to(next_te)) {
      ts = next_ts->left;
      te = next_te->left;
    } else {
      return std::pair<int, int>(ts->sample_ID, j);
    }
    // if (feasible) {
    // } else {
    //   return std::pair<int, int>(ts->sample_ID, j);
    // }
  }
}

    // def longest_match(self, g):
    //     """
    //     Find longest suffix match for g in n-idp time
    //     """
    //     assert len(g) == self.M
    //     a = 0
    //     b = self.N - 1
    //     s = 0
    //     for j in range(self.M - 1, -1, -1):
    //         x = g[j]
    //         # this is the lfx part of loop
    //         ts = self.U[j, a, x]
    //         te = self.V[j, b, x]
    //         if self.feasible(ts, te):  # Test if interval is possible
    //             a = self.LF[j, ts]
    //             b = self.LF[j, te]
    //         else:
    //             s = j + 1
    //             break
    //     k = self.A[s, a]
    //     return k, s
