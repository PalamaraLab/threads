
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

VCFWriter::VCFWriter(ThreadingInstructions& instructions) 
    : gt_iterator(GenotypeIterator(instructions)) {
}

void VCFWriter::write_vcf() {
    int num_dip_samples = gt_iterator.num_samples / 2;
    write_header();
    std::vector<int> gt;
    while (gt_iterator.has_next_genotype()) {
        // Output variant metadata:
        // Chromosome
        std::cout << "1\t";
        // Position
        std::cout << gt_iterator.current_position << "\t";
        // ID
        std::cout << "1:" << gt_iterator.current_position << "\t";
        // REF
        std::cout << "A\t";
        // ALT
        std::cout << "G\t";
        // QUAL
        std::cout << ".\t";
        // FILTER
        std::cout << "PASS\t";
        // INFO
        std::cout << "NS=" << num_dip_samples << "\t";
        // FORMAT
        std::cout << "GT\t";

        // Output genotype
        int i = 0;
        for (auto& gt : gt_iterator.next_genotype()) {
            if (i % 2) {
                std::cout << gt << "\t";
            } else {
                std::cout << gt << "|";
            }
            i++;
        }
        std::cout << "\n";
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
    for (int i=0; i < gt_iterator.num_samples; i = i+2) {
        std::cout << "\t" << "sample_" << i /2;
    }
    std::cout << "\n";
}
