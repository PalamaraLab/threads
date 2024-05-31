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

#ifndef THREADS_ARG_STATE_HPP
#define THREADS_ARG_STATE_HPP

#include "Node.hpp"

#include <unordered_map>
#include <vector>

class TracebackState {
public:
  TracebackState(int _site, Node* _best_prev_node, TracebackState* _prev);

public:
  int site = 0;
  Node* best_prev_node = nullptr; ///< ID to trace back through to last traceback state
  TracebackState* prev = nullptr;
};

class State {
public:
  State(Node* _below, double _score, TracebackState* _traceback);

  /// Shorthand for this.below->site
  bool genotype() const;

  friend std::ostream& operator<<(std::ostream& os, State& state);

public:
  // Nothing here can be const because we want to use std::sort later on
  Node* below = nullptr; ///< Panel entry directly below
  double score = 0.0;    ///< Score of the current state
  TracebackState* traceback = nullptr;
};

class StatePair {
public:
  StatePair(Node* _below_a, Node* _below_b, double _score, TracebackState* _traceback_a,
            TracebackState* _traceback_b);

  friend std::ostream& operator<<(std::ostream& os, StatePair& state_pair);

public:
  Node* below_a = nullptr; ///< Panel entry directly below
  Node* below_b = nullptr;
  double score = 0.0;                    ///< Score of the current state
  TracebackState* traceback_a = nullptr; ///< Pointer to last recombinant state
  TracebackState* traceback_b = nullptr;
};

class StateBranch {
public:
  void insert(const State& state);
  void prune();

public:
  std::vector<State> states;
};

class StateTree {
public:
  StateTree(std::vector<State>& states);
  void prune();
  std::vector<State> dump() const;

public:
  std::unordered_map<int, StateBranch> branches;
};

#endif // THREADS_ARG_STATE_HPP
