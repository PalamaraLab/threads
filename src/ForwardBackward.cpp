// This file is part of the Threads software suite.
// Copyright (C) 2024-2025 Threads Developers.
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

#include "ForwardBackward.hpp"

std::pair<std::vector<double>, std::vector<double>>
forwards_ls_hap(int n, int m, const double* H, const double* s,
                const double* e, const double* r) {
    std::vector<double> F(m * n, 0.0);
    std::vector<double> c(m, 0.0);
    const double inv_n = 1.0 / n;

    // Initialization (l = 0)
    for (int i = 0; i < n; ++i) {
        int match = (H[i] == s[0]) || (s[0] == FWBW_MISSING);
        F[i] = inv_n * e[match];
        c[0] += F[i];
    }
    double inv_c0 = 1.0 / c[0];
    for (int i = 0; i < n; ++i) {
        F[i] *= inv_c0;
    }

    // Forward pass
    for (int l = 1; l < m; ++l) {
        const double r_l = r[l];
        const double one_minus_r = 1.0 - r_l;
        const double r_n = r_l * inv_n;
        const double s_l = s[l];

        for (int i = 0; i < n; ++i) {
            double f_val = F[(l - 1) * n + i] * one_minus_r + r_n;
            int match = (H[l * n + i] == s_l) || (s_l == FWBW_MISSING);
            f_val *= e[l * 2 + match];
            F[l * n + i] = f_val;
            c[l] += f_val;
        }

        double inv_cl = 1.0 / c[l];
        for (int i = 0; i < n; ++i) {
            F[l * n + i] *= inv_cl;
        }
    }

    return {std::move(F), std::move(c)};
}

std::vector<double>
backwards_ls_hap(int n, int m, const double* H, const double* s,
                 const double* e, const double* c, const double* r) {
    std::vector<double> B(m * n, 0.0);
    const double inv_n = 1.0 / n;

    // Initialization (l = m-1)
    for (int i = 0; i < n; ++i) {
        B[(m - 1) * n + i] = 1.0;
    }

    // Backward pass
    std::vector<double> tmp_B(n);
    for (int l = m - 2; l >= 0; --l) {
        double tmp_B_sum = 0.0;
        const double s_lp1 = s[l + 1];

        for (int i = 0; i < n; ++i) {
            int match = (H[(l + 1) * n + i] == s_lp1) || (s_lp1 == FWBW_MISSING);
            tmp_B[i] = e[(l + 1) * 2 + match] * B[(l + 1) * n + i];
            tmp_B_sum += tmp_B[i];
        }

        const double r_lp1 = r[l + 1];
        const double r_n_lp1 = r_lp1 * inv_n;
        const double one_minus_r_lp1 = 1.0 - r_lp1;
        const double inv_c_lp1 = 1.0 / c[l + 1];

        for (int i = 0; i < n; ++i) {
            B[l * n + i] = (r_n_lp1 * tmp_B_sum + one_minus_r_lp1 * tmp_B[i]) * inv_c_lp1;
        }
    }

    return B;
}
