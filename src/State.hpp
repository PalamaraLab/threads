#include "Node.hpp"

#include <vector>
#include <unordered_map>

class TracebackState {
public:
  int site;
  // ID to trace back through to last traceback state
  // int best_prev_ID;
  Node* best_prev_node;
  TracebackState* prev;

  TracebackState(int _site, Node* _best_prev_node, TracebackState* _prev);
};

class State {
public:
  // Nothing here can be const because we want to use std::sort later on
  // Panel entry directly below
  Node* below;
  // Score of the current state
  double score;
  // Pointer to last recombinant state
  TracebackState* traceback;
  // Shorthand for this.below->site
  bool genotype();

  State(Node* _below, double _score, TracebackState* _traceback);

  friend ostream& operator<<(ostream& os, State& state);
};

class StatePair {
public:
  // Panel entry directly below
  Node* below_a;
  Node* below_b;
  // Score of the current state
  double score;
  // Pointer to last recombinant state
  TracebackState* traceback_a;
  TracebackState* traceback_b;

  StatePair(Node* _below_a, Node* _below_b, double _score, TracebackState* _traceback_a,
            TracebackState* _traceback_b);

  friend ostream& operator<<(ostream& os, StatePair& state_pair);
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

