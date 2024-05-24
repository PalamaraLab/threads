#ifndef THREADS_ARG_STATE_HPP
#define THREADS_ARG_STATE_HPP

#include "Node.hpp"

#include <vector>
#include <unordered_map>

class TracebackState {
public:
  int site = 0;
  // ID to trace back through to last traceback state
  Node* best_prev_node = nullptr;
  TracebackState* prev = nullptr;

  TracebackState(int _site, Node* _best_prev_node, TracebackState* _prev);
};

class State {
public:
  // Nothing here can be const because we want to use std::sort later on
  // Panel entry directly below
  Node* below = nullptr;
  // Score of the current state
  double score = 0.0;
  // Pointer to last recombinant state
  TracebackState* traceback = nullptr;
  // Shorthand for this.below->site
  bool genotype();

  State(Node* _below, double _score, TracebackState* _traceback);

  friend std::ostream& operator<<(std::ostream& os, State& state);
};

class StatePair {
public:
  // Panel entry directly below
  Node* below_a = nullptr;
  Node* below_b = nullptr;
  // Score of the current state
  double score = 0.0;
  // Pointer to last recombinant state
  TracebackState* traceback_a = nullptr;
  TracebackState* traceback_b = nullptr;

  StatePair(Node* _below_a, Node* _below_b, double _score, TracebackState* _traceback_a,
            TracebackState* _traceback_b);

  friend std::ostream& operator<<(std::ostream& os, StatePair& state_pair);
};

class StateBranch {
public:
  std::vector<State> states;
  void insert(const State& state);
  void prune();
};

class StateTree {
public:
  std::unordered_map<int, StateBranch> branches;

  StateTree(std::vector<State>& states);
  void prune();
  std::vector<State> dump();
};

#endif // THREADS_ARG_STATE_HPP
