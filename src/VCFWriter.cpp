
// This file is part of the Threads software suite.
// Copyright (C) 2024 Threads Developers.
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
    std::vector<int> gt;
    while (gt_iterator.has_next_genotype()) {
        // Output variant metadata:
        // Chromosome
        std::cout << chrom.at(i) << "\t";
        // Position
        std::cout << pos.at(i) << "\t";
        // ID
        std::cout << id.at(i) << "\t";
        // REF
        std::cout << ref.at(i) << "\t";
        // ALT
        std::cout << alt.at(i) << "\t";
        // QUAL
        std::cout << qual.at(i) << "\t";
        // FILTER
        std::cout << filter.at(i) << "\t";
        // INFO
        std::cout << "NS=" << num_dip_samples << "\t";
        // FORMAT
        std::cout << "GT\t";

        // Output genotype
        int j = 0;
        for (const auto& g : gt_iterator.next_genotype()) {
            if (j % 2) {
                std::cout << g << "\t";
            } else {
                std::cout << g << "|";
            }
            j++;
        }
        std::cout << "\n";
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
    for (auto name : sample_names) {
        std::cout << "\t" << name;
    }
    std::cout << "\n";
}
