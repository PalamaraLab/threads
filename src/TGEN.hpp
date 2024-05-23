#ifndef THREADS_ARG_TGEN_HPP
#define THREADS_ARG_TGEN_HPP

#include <Eigen/Dense>
#include <limits>
#include <memory>
#include <vector>


class TGENImpl; // Forward declaration of the implementation class
// Using PIMPL idiom to prevent libraries that #include "TGEN.hpp" needing to see boost/icl headers

class TGEN {

private:
  std::unique_ptr<TGENImpl> pimpl; // Pointer to implementation

public:
  // Constructors
  TGEN(std::vector<int> _positions, std::vector<std::vector<int>> _bp_starts,
       std::vector<std::vector<int>> _target_IDs, std::vector<std::vector<int>> _het_sites);

  ~TGEN(); // Define the destructor for proper deletion of pimpl

  std::vector<std::vector<bool>>& query(int start_pos, int end_pos, const std::vector<int>& samples);

  // Deprecated Eigen version:
  // Eigen::MatrixXi& query(const int start_pos, const int end_pos, const std::vector<int>& samples);
  void clear_cache();
};

#endif // THREADS_ARG_TGEN_HPP
