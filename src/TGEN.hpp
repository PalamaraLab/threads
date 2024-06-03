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

#ifndef THREADS_ARG_TGEN_HPP
#define THREADS_ARG_TGEN_HPP

#include <limits>
#include <memory>
#include <vector>

// Forward declaration of the implementation class
// Using PIMPL idiom to prevent libraries that #include "TGEN.hpp" needing to see boost/icl headers
class TGENImpl;

class TGEN {
public:
  TGEN(std::vector<int> _positions, std::vector<std::vector<int>> _bp_starts,
       std::vector<std::vector<int>> _target_IDs, std::vector<std::vector<int>> _het_sites);
  ~TGEN();

  std::vector<std::vector<bool>>& query(int start_pos, int end_pos,
                                        const std::vector<int>& samples);
  void clear_cache();

private:
  std::unique_ptr<TGENImpl> pimpl; ///< Pointer to implementation
};

#endif // THREADS_ARG_TGEN_HPP
