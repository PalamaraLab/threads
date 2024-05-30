#ifndef THREADS_ARG_TGEN_SEGMENT_HPP
#define THREADS_ARG_TGEN_SEGMENT_HPP

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

// Here is a typical class that may model intervals in your application.
class TgenSegment {
public:
  TgenSegment() : target(std::numeric_limits<int>::max()), first(), past() {
  }

  TgenSegment(const TgenSegment& other)
      : het_sites(other.het_sites), target(other.target), first(other.first), past(other.past) {
  }

  TgenSegment(int lo, int up) : target(std::numeric_limits<int>::max()), first(lo), past(up) {
  }

  TgenSegment(int lo, int up, std::vector<int> _hets, int _target)
      : het_sites(_hets), target(_target), first(lo), past(up) {
  }

  [[nodiscard]] int lower() const {
    return first;
  }

  [[nodiscard]] int upper() const {
    return past;
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

  TgenSegment calc_intersection_with(const TgenSegment& other) const {
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

public:
  std::vector<int> het_sites;
  int target = std::numeric_limits<int>::max();

private:
  int first = 0;
  int past = 0;
};

#endif // THREADS_ARG_TGEN_SEGMENT_HPP
