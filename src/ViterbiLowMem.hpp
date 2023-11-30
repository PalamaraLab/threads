#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class TracebackNode {
private:
  inline size_t coord_id_key(int i, int j) {
    return (size_t) i << 32 | (unsigned int) j;
  }

public:
  int sample_id;
  int site;
  double score;
  TracebackNode* previous;
  TracebackNode(int _sample_id, int _site, TracebackNode* _previous, double _score);
  size_t key();
};

class ViterbiPath {
private:
public:
  double score;
  int target_id;
  std::vector<int> bp_starts;
  std::vector<int> segment_starts;
  std::vector<int> sample_ids;
  std::vector<double> heights;
  std::vector<int> het_sites;
  ViterbiPath(int _target_id);
  ViterbiPath(int _target_id, std::vector<int> _segment_starts, std::vector<int> _sample_ids,
              std::vector<double> _heights, std::vector<int> _het_sites);
  void append(int segment_start, int sample_id);
  void append(int segment_start, int sample_id, double height, std::vector<int>& new_het_sites);
  std::tuple<std::vector<int>, std::vector<int>, std::vector<double>>
  dump_data_in_range(int start_idx, int end_idx);
  void reverse();
  int size();
  void map_positions(std::vector<int>& positions);
};

class ViterbiState {

private:
  inline size_t coord_id_key(int i, int j) {
    return (size_t) i << 32 | (unsigned int) j;
  }
  // don't think we need unique_ptr any more here
  std::unordered_map<size_t, TracebackNode> traceback_states;
  TracebackNode* recursive_insert(std::unordered_map<size_t, TracebackNode>& state_map,
                                  TracebackNode* state);

public:
  int target_id;
  int best_match = -1;
  double best_score;
  int sites_processed = 0;
  double mutation_penalty;
  std::vector<int> sample_ids;
  std::vector<double> sample_scores;
  std::unordered_map<int, TracebackNode*> current_tracebacks;

  ViterbiState(int _target_id, std::vector<int> _sample_ids);

  void process_site(const std::vector<int>& genotype, double rho, double rho_c, double _mu,
                    double _mu_c);
  void set_samples(std::unordered_set<int> new_sample_ids);
  int count_branches();
  void prune();
  ViterbiPath traceback();
};
