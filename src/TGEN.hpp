#include <Eigen/Dense>
#include <boost/icl/separate_interval_set.hpp>
#include <iostream>
#include <limits>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;

// Here is a typical class that may model intervals in your application.
class TgenSegment {
public:
  TgenSegment() : _first(), _past(), target(std::numeric_limits<int>::max()) {
  }
  TgenSegment(const TgenSegment& other)
      : _first(other._first), _past(other._past), het_sites(other.het_sites), target(other.target) {
  }
  TgenSegment(int lo, int up) : _first(lo), _past(up), target(std::numeric_limits<int>::max()) {
  }
  TgenSegment(int lo, int up, std::vector<int> _hets, int _target)
      : _first(lo), _past(up), het_sites(_hets), target(_target) {
  }

  std::vector<int> het_sites;
  int target;

  int lower() const {
    return _first;
  }
  int upper() const {
    return _past;
  }

  friend std::ostream& operator<<(std::ostream& os, const TgenSegment& seg) {
    os << "Segment with range [" << seg.lower() << ", " << seg.upper() << ") with target "
       << seg.target;
    os << " and het sites ";
    for (int i : seg.het_sites) {
      os << i << " ";
    }
    return os;
  }

  TgenSegment operator&(const TgenSegment& other) const {
    int new_first = std::max(lower(), other.lower());
    int new_past = std::min(upper(), other.upper());

    std::vector<int> merged_sites;
    auto this_start = std::lower_bound(het_sites.begin(), het_sites.end(), new_first);
    auto this_end = std::lower_bound(het_sites.begin(), het_sites.end(), new_past);
    auto other_start = std::lower_bound(other.het_sites.begin(), other.het_sites.end(), new_first);
    auto other_end = std::lower_bound(other.het_sites.begin(), other.het_sites.end(), new_past);
    std::merge(this_start, this_end, other_start, other_end, std::back_inserter(merged_sites));
    return TgenSegment(new_first, new_past, merged_sites, std::min(target, other.target));
  }

private:
  int _first, _past;
};

namespace boost {
namespace icl {
// See
// https://www.boost.org/doc/libs/1_81_0/libs/icl/doc/html/boost_icl/examples/custom_interval.html
// Class template interval_traits serves as adapter to register and customize your interval class
template <>
struct interval_traits<TgenSegment>      // 1.  Partially specialize interval_traits
{                                        // 2.  Define associated types
  typedef TgenSegment interval_type;     // 2.1 TgenSegment will be the interval_type
  typedef int domain_type;               // 2.2 The elements of the domain are ints
  typedef std::less<int> domain_compare; // 2.3 This is the way our element shall be ordered.
                                         // 3.  Next we define the essential functions
                                         //    of the specialisation
                                         // 3.1 Construction of intervals
  static interval_type construct(const domain_type& lo, const domain_type& up) {
    return interval_type(lo, up);
  }
  // 3.2 Selection of values
  static domain_type lower(const interval_type& inter_val) {
    return inter_val.lower();
  };
  static domain_type upper(const interval_type& inter_val) {
    return inter_val.upper();
  };
};

template <>
struct interval_bound_type<TgenSegment> // 4.  Finally we define the interval borders.
{                                       //    Choose between static_open         (lo..up)
  typedef interval_bound_type type;     //                   static_left_open    (lo..up]
  BOOST_STATIC_CONSTANT(bound_type, value = interval_bounds::static_right_open);
}; //               and static_closed       [lo..up]

} // namespace icl
} // namespace boost

typedef boost::icl::separate_interval_set<int, std::less, TgenSegment> SegmentSet;

class TGEN {
private:
  std::vector<SegmentSet> interval_sets;

  std::unordered_map<int, int> pos_idx_map;
  Eigen::VectorXi reference_genome;
  std::vector<bool> reference_genome_vec;
  std::vector<std::vector<int>> bp_starts;
  std::vector<std::vector<int>> target_IDs;
  std::vector<std::vector<int>> het_sites;
  std::vector<int> positions;
  std::vector<std::vector<bool>> genotypes;

  Eigen::MatrixXi genotype_cache;
  std::unordered_map<int, int> cached_genotypes_map;
  void read_interval(const int sample, const int current_sample, const int start_pos_idx,
                     const int end_pos_idx, const std::vector<int> het_sites_range,
                     const int start_pos_idx_offset);

public:
  // Constructors
  TGEN(std::vector<int> _positions, std::vector<std::vector<int>> _bp_starts,
       std::vector<std::vector<int>> _target_IDs, std::vector<std::vector<int>> _het_sites);

  std::vector<std::vector<bool>>& query(const int start_pos, const int end_pos, const std::vector<int>& samples);

  // Deprecated Eigen version:
  // Eigen::MatrixXi& query(const int start_pos, const int end_pos, const std::vector<int>& samples);
  void clear_cache();
};
