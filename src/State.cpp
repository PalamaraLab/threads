#include "State.hpp"

#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>
using std::cout;
using std::endl;

TracebackState::TracebackState(int _site, int _best_prev_ID, TracebackState* _prev) :
  site(_site), best_prev_ID(_best_prev_ID), prev(_prev) {
}

State::State(Node* _below, double _score, TracebackState* _traceback) :
  below(_below), score(_score), traceback(_traceback) {
}

bool State::genotype() {
  return this->below->above->genotype;
}

ostream& operator<<(ostream& os, State& state) {
  os << "State between nodes " << state.below->sample_ID;
  os << " and " << state.below->above->sample_ID; 
  os << " with score " << state.score;
  os << ", traceback to " << state.traceback->site;
  os << ", genotype " << state.genotype();
  return os;
}

void StateBranch::insert(const State& state) {
  states.push_back(state);
}

void StateBranch::prune() {
  std::sort(states.begin(), states.end(), [](const State& s1, const State& s2) {
    return (s1.traceback->site == s2.traceback->site) ? s1.score < s2.score : s1.traceback->site > s2.traceback->site;
  });

  std::vector<State> new_states;
  double running_min_score = std::numeric_limits<double>::infinity();
  for (const State& s : states) {
    if (s.score < running_min_score) {
      new_states.push_back(s);
      running_min_score = s.score;
    }
  }
  states = new_states;
}

StateTree::StateTree(std::vector<State>& states) {
  for (auto s : states) {
    int sample_ID = s.below->sample_ID;
    if (branches.find(sample_ID) == branches.end()) {
      branches[sample_ID] = StateBranch();
    }
    branches[sample_ID].insert(s);
  }
}

void StateTree::prune() {
  for (auto pair : branches) {
    branches[pair.first].prune();
  }
}

std::vector<State> StateTree::dump() {
  std::vector<State> states;
  for (auto pair : branches) {
    for (auto s : pair.second.states) {
      states.push_back(s);
    }
  }
  return states;
}
