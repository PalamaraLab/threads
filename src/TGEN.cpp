#include "TGEN.hpp"

#include "TgenSegment.hpp"

#include <Eigen/Core>
#include <boost/icl/separate_interval_set.hpp>

#include <limits>
#include <queue>
#include <utility>
#include <vector>


namespace boost::icl {
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

} // namespace boost::icl

typedef boost::icl::separate_interval_set<int, std::less, TgenSegment> SegmentSet;

class TGENImpl {
public:
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

  // Constructor
  TGENImpl(std::vector<int> _positions, std::vector<std::vector<int>> _bp_starts,
           std::vector<std::vector<int>> _target_IDs, std::vector<std::vector<int>> _het_sites)
      : bp_starts(std::move(_bp_starts)), target_IDs(std::move(_target_IDs)),
        het_sites(std::move(_het_sites)), positions(std::move(_positions)) {
    positions.push_back(std::numeric_limits<int>::max());
    // position-to-site map
    for (int i = 0; i < positions.size(); i++) {
      pos_idx_map[positions[i]] = i;
    }

    // set reference genome here
    reference_genome = Eigen::VectorXi::Zero(positions.size());
    reference_genome_vec = std::vector<bool>(positions.size(), false);
    for (int h : het_sites[0]) {
      reference_genome[pos_idx_map.at(h)] = 1;
      reference_genome_vec[pos_idx_map.at(h)] = true;
    }

    interval_sets.reserve(bp_starts.size());
    // Initialize interval maps for each sample. This can get too slow
    for (int i = 0; i < bp_starts.size(); i++) {
      SegmentSet iset;
      int n_segs = bp_starts[i].size();
      std::vector<int>& sample_hets = het_sites[i];
      int n_hets = sample_hets.size();
      int het_site_idx = 0;
      int pos_idx = 0;
      std::vector<int> seg_hets;
      for (int j = 0; j < n_segs; j++) {
        int seg_start = bp_starts[i][j];
        int seg_end = j == n_segs - 1 ? positions.back() : bp_starts[i][j + 1];

        // Binary search through positions to get segment boundaries
        int seg_start_pos = *std::lower_bound(positions.begin(), positions.end(), seg_start);
        int seg_end_pos = j == n_segs - 1 ? positions.back()
                                          : *std::lower_bound(positions.begin(), positions.end(),
                                                              bp_starts[i][j + 1]);

        while (het_site_idx < n_hets && seg_start <= sample_hets[het_site_idx] &&
               sample_hets[het_site_idx] < seg_end) {
          seg_hets.push_back(sample_hets[het_site_idx]);
          het_site_idx++;
        }
        iset += TgenSegment(seg_start_pos, seg_end_pos, seg_hets, target_IDs[i][j]);
      }
      interval_sets.push_back(iset);
    }

    // initialize cached genotypes
    cached_genotypes_map[0] = -1;
  }

  // Eigen-based query
  Eigen::MatrixXi& query(const int bp_from, const int bp_to, const std::vector<int>& samples) {
    clear_cache();

// Deprecated eigen-based query
// Eigen::MatrixXi& TGEN::query(const int bp_from, const int bp_to, const std::vector<int>& samples) {
//   clear_cache();

//   // Find number of expected sites
//   int start_pos = *std::lower_bound(positions.begin(), positions.end(), bp_from);
//   int end_pos = *std::upper_bound(positions.begin(), positions.end(), bp_to);
//   int idx_offset = pos_idx_map[start_pos];
//   genotype_cache.resize(samples.size(), pos_idx_map[end_pos] - idx_offset);

//   TgenSegment range(start_pos, end_pos);

//   for (int i = 0; i < samples.size(); i++) {
//     cached_genotypes_map[samples[i]] = i;
//     if (samples[i] == 0) {
//       auto insert_range =
//           Eigen::seq(0, pos_idx_map[end_pos] - idx_offset - 1); // eigen seq is inclusive
//       auto copy_range = Eigen::seq(idx_offset, pos_idx_map[end_pos] - 1);
//       genotype_cache(i, insert_range) = reference_genome(
//           copy_range); //.eval(); //WARNING need .eval() here (or do we? I don't think we do)
//     }
//     else {
//       SegmentSet& segments(interval_sets[samples[i]]);

//       // Initialize the queue
//       std::queue<TgenSegment> seg_queue;
//       auto eqr = segments.equal_range(range);
//       for (SegmentSet::const_iterator iter = eqr.first; iter != eqr.second; iter++) {
//         seg_queue.push(*iter & range);
//       }

//       // Process everything in the queue
//       while (!seg_queue.empty()) {
//         TgenSegment& segment = seg_queue.front();
//         if (cached_genotypes_map.find(segment.target) != cached_genotypes_map.end()) {
//           // We've reached somewhere along the tree where we can copy from
//           int seg_start_idx = pos_idx_map[segment.lower()];
//           int seg_end_idx = pos_idx_map[segment.upper()];
//           auto insert_range = Eigen::seq(
//               seg_start_idx - idx_offset, seg_end_idx - idx_offset - 1); // eigen seq is inclusive
//           if (segment.target == 0) {
//             auto copy_range = Eigen::seq(seg_start_idx, seg_end_idx - 1);
//             // We've reached the root of the tree and copy from the "reference" genome
//             genotype_cache(i, insert_range) = reference_genome(copy_range);
//           }
//           else {
//             // We've found a cached genotype to copy from

//             genotype_cache(i, insert_range) =
//                 genotype_cache(cached_genotypes_map[segment.target], insert_range);
//           }

//           // We then flip all the het sites
//           for (int h : segment.het_sites) {
//             genotype_cache(i, pos_idx_map[h] - idx_offset) =
//                 1 - genotype_cache(i, pos_idx_map[h] - idx_offset);
//           }
//         }
//         else {
//           // We've not yet reached somewhere to copy from, so we keep traversing
//           auto new_eqr = interval_sets[segment.target].equal_range(segment);
//           for (SegmentSet::const_iterator iter = new_eqr.first; iter != new_eqr.second; iter++) {
//             seg_queue.push(*iter & segment);
//           }
//         }
//         seg_queue.pop();
//       }
//     }
//   }
//   return genotype_cache;
// }

// std::vector-based query
// Warning: This makes a copy when returned through the python interface
std::vector<std::vector<bool>>& TGEN::query(const int bp_from, const int bp_to,
                                             const std::vector<int>& samples) {
  genotypes.clear();

    // Find number of expected sites
    int start_pos = *std::lower_bound(positions.begin(), positions.end(), bp_from);
    int end_pos = *std::upper_bound(positions.begin(), positions.end(), bp_to);
    int idx_offset = pos_idx_map[start_pos];

    int n_samples = samples.size();
    int n_sites = pos_idx_map[end_pos] - idx_offset;
    genotypes.reserve(samples.size());
    for (int i = 0; i < samples.size(); i++) {
      genotypes.push_back(std::vector<bool>(n_sites));
    }
    TgenSegment range(start_pos, end_pos);

    for (int i = 0; i < samples.size(); i++) {
      std::vector<bool>& current_gt = genotypes.at(i);
      cached_genotypes_map[samples[i]] = i;
      if (samples[i] == 0) {
        std::copy(reference_genome_vec.begin() + idx_offset,
                  reference_genome_vec.begin() + pos_idx_map[end_pos], current_gt.begin());
      }
      else {
        SegmentSet& segments(interval_sets[samples[i]]);

        // Initialize the queue
        std::queue<TgenSegment> seg_queue;
        auto eqr = segments.equal_range(range);
        for (SegmentSet::const_iterator iter = eqr.first; iter != eqr.second; iter++) {
          seg_queue.push(*iter & range);
        }

        // Process everything in the queue
        while (!seg_queue.empty()) {
          TgenSegment& segment = seg_queue.front();
          if (cached_genotypes_map.find(segment.target) != cached_genotypes_map.end()) {
            // We've reached somewhere along the tree where we can copy from
            int seg_start_idx = pos_idx_map[segment.lower()];
            int seg_end_idx = pos_idx_map[segment.upper()];
            if (segment.target == 0) {
              // We've reached the root of the tree and copy from the "reference" genome
              std::copy(reference_genome_vec.begin() + seg_start_idx,
                        reference_genome_vec.begin() + seg_end_idx,
                        current_gt.begin() + seg_start_idx - idx_offset);
            }
            else {
              std::vector<bool>& target_gt = genotypes.at(cached_genotypes_map[segment.target]);
              std::copy(target_gt.begin() + seg_start_idx - idx_offset,
                        target_gt.begin() + seg_end_idx - idx_offset,
                        current_gt.begin() + seg_start_idx - idx_offset);
              // We've found a cached genotype to copy from
            }
            // We then flip all the het sites
            // then we can also delay the .eval() until the end, right?
            for (int h : segment.het_sites) {
              current_gt.at(pos_idx_map[h] - idx_offset) =
                  !current_gt.at(pos_idx_map[h] - idx_offset);
            }
          }
          else {
            // We've not yet reached somewhere to copy from, so we keep traversing
            auto new_eqr = interval_sets[segment.target].equal_range(segment);
            for (SegmentSet::const_iterator iter = new_eqr.first; iter != new_eqr.second; iter++) {
              seg_queue.push(*iter & segment);
            }
          }
          seg_queue.pop();
        }
      }
    }
    return genotypes;
  }

  void clear_cache() {
    cached_genotypes_map.clear();
    genotype_cache.resize(0, 0);
    cached_genotypes_map[0] = -1;
  }
};

TGEN::TGEN(std::vector<int> _positions, std::vector<std::vector<int>> _bp_starts,
           std::vector<std::vector<int>> _target_IDs, std::vector<std::vector<int>> _het_sites)
    : pimpl(std::make_unique<TGENImpl>(std::move(_positions), std::move(_bp_starts),
                                       std::move(_target_IDs), std::move(_het_sites))) {}

TGEN::~TGEN() = default;

Eigen::MatrixXi& TGEN::query(const int start_pos, const int end_pos, const std::vector<int>& samples) {
  // Forward the call to pimpl
  return pimpl->query(start_pos, end_pos, samples);
}

std::vector<std::vector<bool>>& TGEN::query2(const int start_pos, const int end_pos, const std::vector<int>& samples) {
  // Forward the call to pimpl
  return pimpl->query2(start_pos, end_pos, samples);
}

void TGEN::clear_cache() {
  // Forward the call to pimpl
  pimpl->clear_cache();
}
