#include "Shuttle.hpp"
#include <unordered_set>
#include <vector>
#include <queue>
#include <Eigen/Core>
#include <iostream>

using Eigen::ArrayXXf;

Shuttle::Shuttle(std::vector<int> init_targets,
                 std::vector<int> _site_first_carriers,
                 std::vector<std::unordered_set<int>>& _site_hets,
                 std::vector<int>& _site_positions,
                 const std::vector<int>& weft_ids,
                 const std::vector<int>& weft_targets,
                 const std::vector<int>& weft_positions) :
        targets(init_targets), site_hets(_site_hets), site_first_carriers(_site_first_carriers), site_positions(_site_positions) {
    //assert here:
    //  lists are all equal size
    //  positions is increasing
    //  site positions are not empty
    //  site positions and site hets and site_first_carriers have same length
    num_sites = site_positions.size();
    num_samples = targets.size();
    // wefts.reserve(weft_ids.size());
    for (int i = 0; i < weft_ids.size(); i++) {
        wefts.emplace_back(weft_ids.at(i), weft_targets.at(i), weft_positions.at(i));
    }
    next_weft_idx = 0;
    current_site = -1;
    // descendants.reserve(num_samples);
    descendants.push_back(std::unordered_set<int>());
    for (int i = 1; i < num_samples; i++) {
        descendants.push_back(std::unordered_set<int>());
        descendants.at(targets.at(i)).insert(i);
    }

    proceed_to_next_site();
}

bool Shuttle::is_finished() {
    return current_site >= num_sites;
}

void Shuttle::proceed_to_next_site() {
    if (current_site >= num_sites - 1) {
        current_site = num_sites;
        return;
    }
    current_site++;
    int pos = site_positions.at(current_site);
    while (next_weft_idx < wefts.size() && wefts.at(next_weft_idx).position <= pos) {
        weave_weft();
    }
    set_current_carriers();
}

void Shuttle::weave_weft() {
    Weft& weft = wefts.at(next_weft_idx);
    descendants.at(targets.at(weft.id)).erase(weft.id);
    descendants.at(weft.target).insert(weft.id);
    targets.at(weft.id) = weft.target;
    next_weft_idx++;
}

std::unordered_set<int> Shuttle::set_current_carriers() {
    current_carriers.clear();
    std::unordered_set<int>& current_hets = site_hets.at(current_site);
    current_carriers.insert(site_first_carriers.at(current_site));

    std::queue<int> carrier_queue;
    for (int d : descendants.at(site_first_carriers.at(current_site))) {
        carrier_queue.push(d);
    }

    while (!carrier_queue.empty()) {
        int c = carrier_queue.front();
        // std::cout << c << "\n";
        carrier_queue.pop();
        if (current_hets.count(c) > 0) {
            continue;
        } else {
            current_carriers.insert(c);
            for (int d : descendants.at(c)) {
                carrier_queue.push(d);
            }
        }
    }
    return current_carriers;
}

void Shuttle::GU_mult(Eigen::Ref<const ArrayXXf>& U, Eigen::Ref<ArrayXXf>& out) {
    // G is the encoded (2N x M)-genotype matrix
    // U must be size (M x K)
    // out must be size (N x K), not "2N" because we assume diploid
    if (U.rows() != num_sites) {
        throw std::runtime_error("Invalid 'U' size");
    }
    if (out.cols() != U.cols() || out.rows() != num_samples / 2) {
        std::cout << "'out' should be " << num_samples/2 << "x" << U.cols() << " but is " << out.rows() << "x" << out.cols() << "\n";
        throw std::runtime_error("Invalid 'out' size");
    }
    

    const std::size_t k = U.cols();
    while (!is_finished()) {
        for (const auto c : current_carriers) {
            int dip_id = c / 2;
            for (int j = 0; j < k; j++) {
                out(dip_id, j) += U(current_site, j);
            }
        }
        proceed_to_next_site();
    }
}
