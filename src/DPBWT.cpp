#include "DPBWT.hpp"

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

  for (int i = 0; i < physical_positions.size(); i++) {
    tops.emplace_back(std::make_unique<Node>(-1, i, 0));
    bottoms.emplace_back(std::make_unique<Node>(-1, i, 1));
    Node* top_i = tops.back().get();
    Node* bottom_i = bottoms.back().get();

    top_i->below = bottom_i;
    bottom_i->above = top_i;

    // top_i->below = bottom_i;
    // bottom_i->above = top_i;
    if (i > 0) {
      tops[i - 1]->w = {top_i, top_i};
      bottoms[i - 1]->w = {bottom_i, bottom_i};
    }
  }
}

void DPBWT::insert(std::vector<bool> genotype) {
  // If nothing in here, make a pointer to all tops,
  // then insert genotype.
  // This is ok but ignores last genotype for some reason and instead has a useless col-0.
  int current_ID = num_samples;
  if (genotype.size() != num_sites) {
    throw std::invalid_argument(THROW_LINE("Number of input markers does not match map."));
  }
  cout << "Doing sample " << current_ID << endl;
  // std::vector<std::unique_ptr<Node>> row;
  panel.emplace_back(std::vector<std::unique_ptr<Node>>());

  Node* t0 = bottoms[0].get();
  // Node z0 = Node(current_ID, 0, genotype[0]);
  panel[current_ID].emplace_back(std::make_unique<Node>(current_ID, 0, genotype[0]));
  Node* z0 = panel[current_ID][0].get();

  // Inserts z0 above t0
  t0->insert_above(z0);

  Node* z_k = z0;
  Node* z_next;
  Node* tmp;
  Node* t_k;
  Node* t_next;
  for (int k = 0; k < num_sites - 1; k++) {
    bool g_k = genotype[k];
    panel[current_ID].emplace_back(std::make_unique<Node>(current_ID, k + 1, genotype[k + 1]));
    z_next = panel[current_ID][k + 1].get();
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
    if (t_k->sample_ID == -1 and !g_k) {
      // Logic is that if we're at the bottom and have a "0" genotype we jump up to last
      // 0-spot in next column
      t_next = tops[k]->below->w[1];
    } else {
      t_next = t_k->w[g_k];
    }
    z_next->sample_ID = current_ID;
    t_next->insert_above(z_next);
    z_k = z_next;
  }
  // Special case for the last column

  num_samples++;
}

void DPBWT::print_sorting() {
  for (int j = 0; j < num_sites; ++j) {
    Node* node = tops[j].get();
    while (node != nullptr) {
      cout << node->sample_ID << " ";
      node = node->below;
    }
    cout << endl;
  }
}
