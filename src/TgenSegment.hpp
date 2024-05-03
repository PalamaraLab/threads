#ifndef THREADS_INFER_TGEN_SEGMENT_HPP
#define THREADS_INFER_TGEN_SEGMENT_HPP

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

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

  [[nodiscard]] int lower() const {
    return _first;
  }
  [[nodiscard]] int upper() const {
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
    return {new_first, new_past, merged_sites, std::min(target, other.target)};
  }

private:
  int _first, _past;
};

#endif // THREADS_INFER_TGEN_SEGMENT_HPP
