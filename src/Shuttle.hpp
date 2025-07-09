/*
 * The idea here is to have an initial tree and a queue of tweaks.
 * Then we traverse the queue and periodically stop to apply tweaks and do whatever operation.
 * We don't need to keep track of coalescence times.
 * But we need mismatches. Like a set of mismatches per variant, for example.
 *
 *
 *
 *
 */

#include <unordered_set>
#include <vector>
#include <Eigen/Core>

using Eigen::ArrayXXf;

class Weft {
public:
    Weft(const int _id, const int _target, const int _position) : id(_id), target(_target), position(_position) {};
    int id;
    int target;
    int position;
};

class Shuttle {
private:
    int next_weft_idx = 0;
    std::vector<Weft> wefts;
    std::vector<int> site_first_carriers;
    std::vector<std::unordered_set<int>> site_hets;
    void weave_weft();
    std::unordered_set<int> set_current_carriers();
public:
    std::unordered_set<int> current_carriers;
    int num_sites = 0;
    int num_samples = 0;
    std::vector<int> site_positions;
    int current_site = -1;
    std::vector<int> targets;
    std::vector<std::unordered_set<int>> descendants;
    void proceed_to_next_site();
    bool is_finished();
    Shuttle(std::vector<int> init_targets,
            std::vector<int> _site_first_carriers,
            std::vector<std::unordered_set<int>>& _site_hets,
            std::vector<int>& _site_positions,
            const std::vector<int>& weft_ids,
            const std::vector<int>& weft_targets,
            const std::vector<int>& weft_positions);
    
    void GU_mult(Eigen::Ref<const ArrayXXf>& U, Eigen::Ref<ArrayXXf>& out);
};
