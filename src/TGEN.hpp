#ifndef THREADS_ARG_TGEN_HPP
#define THREADS_ARG_TGEN_HPP

#include <Eigen/Dense>
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
