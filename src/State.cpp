// This file is part of the Threads software suite.
// Copyright (C) 2024 Threads Developers.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "State.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

TracebackState::TracebackState(int _site, Node* _best_prev_node, TracebackState* _prev)
    : site(_site), best_prev_node(_best_prev_node), prev(_prev) {
}

State::State(Node* _below, double _score, TracebackState* _traceback)
    : below(_below), score(_score), traceback(_traceback) {
}

bool State::genotype() const {
  return this->below->above->genotype;
}

std::ostream& operator<<(std::ostream& os, State& state) {
  os << "State with nodes (" << state.below->sample_ID;
  os << ") and (" << state.below->above->sample_ID;
  os << ") with score " << state.score;
  os << ", traceback to " << state.traceback->site;
  os << ", genotype " << state.genotype();
  return os;
}

StatePair::StatePair(Node* _below_a, Node* _below_b, double _score, TracebackState* _traceback_a,
                     TracebackState* _traceback_b)
    : below_a(_below_a), below_b(_below_b), score(_score), traceback_a(_traceback_a),
      traceback_b(_traceback_b) {
}

std::ostream& operator<<(std::ostream& os, StatePair& pair) {
  os << "State-pair between a: nodes " << pair.below_a->sample_ID;
  os << ", " << pair.below_a->above->sample_ID;
  os << " and b: nodes " << pair.below_b->sample_ID;
  os << ", " << pair.below_b->above->sample_ID;
  os << " with score " << pair.score;
  os << ", traceback to " << pair.traceback_a->site;
  os << " and " << pair.traceback_b->site;
  return os;
}

void StateBranch::insert(const State& state) {
  states.push_back(state);
}

void StateBranch::prune() {
  std::sort(states.begin(), states.end(), [](const State& s1, const State& s2) {
    return (s1.traceback->site == s2.traceback->site) ? s1.score < s2.score
                                                      : s1.traceback->site > s2.traceback->site;
  });

  std::vector<State> new_states;
  double running_min_score = std::numeric_limits<double>::max();
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

std::vector<State> StateTree::dump() const {
  std::vector<State> states;
  for (auto pair : branches) {
    for (auto s : pair.second.states) {
      states.push_back(s);
    }
  }
  return states;
}
