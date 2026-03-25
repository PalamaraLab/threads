
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

#include "VCFWriter.hpp"

#include <iostream>
#include <vector>

VCFWriter::VCFWriter(const ThreadingInstructions& instructions)
    : gt_iterator(GenotypeIterator(instructions)) {
}

void VCFWriter::write_vcf() {
    int i = 0;
    int num_dip_samples = gt_iterator.num_samples / 2;
    write_header();

    // Pre-size line buffer: metadata (~200 chars) + genotypes (4 chars per diploid sample)
    std::string line;
    line.reserve(256 + num_dip_samples * 4);

    while (gt_iterator.has_next_genotype()) {
        line.clear();

        // Variant metadata
        line += chrom[i]; line += '\t';
        line += pos[i]; line += '\t';
        line += id[i]; line += '\t';
        line += ref[i]; line += '\t';
        line += alt[i]; line += '\t';
        line += qual[i]; line += '\t';
        line += filter[i]; line += '\t';
        line += "NS=";
        line += std::to_string(num_dip_samples);
        line += "\tGT";

        // Genotypes: phased diploid (0|0, 0|1, etc.)
        const auto& geno = gt_iterator.next_genotype();
        int n = gt_iterator.num_samples;
        for (int j = 0; j < n; j += 2) {
            line += '\t';
            line += ('0' + geno[j]);
            line += '|';
            line += ('0' + geno[j + 1]);
        }
        line += '\n';

        std::cout.write(line.data(), line.size());
        i++;
    }
}

void VCFWriter::write_header() {
    std::cout << "##fileformat=VCFv4.2\n"
        << "##source=Threads\n"
        << "##contig=<ID=1,length=" << gt_iterator.positions.back() << ">\n"
        << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& name : sample_names) {
        std::cout << '\t' << name;
    }
    std::cout << '\n';
}
