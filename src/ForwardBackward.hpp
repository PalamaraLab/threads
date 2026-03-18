// This file is part of the Threads software suite.
// Copyright (C) 2024-2025 Threads Developers.
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

#ifndef THREADS_ARG_FORWARD_BACKWARD_HPP
#define THREADS_ARG_FORWARD_BACKWARD_HPP

#include <vector>
#include <utility>

constexpr int FWBW_MISSING = -9;

// Li-Stephens haploid forward algorithm (normalized).
// H: reference panel, row-major (m x n)
// s: query haplotype (length m)
// e: emission probabilities, row-major (m x 2), e[l*2+0]=mismatch, e[l*2+1]=match
// r: recombination probabilities (length m)
// Returns (F, c) where F is (m x n) forward matrix and c is (m) normalization factors.
std::pair<std::vector<double>, std::vector<double>>
forwards_ls_hap(int n, int m, const double* H, const double* s,
                const double* e, const double* r);

// Li-Stephens haploid backward algorithm (normalized).
// Same inputs as forward, plus c (normalization factors from forward pass).
// Returns B: (m x n) backward matrix.
std::vector<double>
backwards_ls_hap(int n, int m, const double* H, const double* s,
                 const double* e, const double* c, const double* r);

#endif // THREADS_ARG_FORWARD_BACKWARD_HPP
