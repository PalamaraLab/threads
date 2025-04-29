
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

#ifndef THREADS_ARG_VCF_WRITER_HPP
#define THREADS_ARG_VCF_WRITER_HPP

#include "ThreadingInstructions.hpp"
#include "GenotypeIterator.hpp"

#include <string>
#include <iostream>
#include <vector>

class VCFWriter {
private:
    GenotypeIterator gt_iterator;
    std::vector<std::string> chrom;
    std::vector<std::string> pos;
    std::vector<std::string> id;
    std::vector<std::string> ref;
    std::vector<std::string> alt;
    std::vector<std::string> qual;
    std::vector<std::string> filter;
    std::vector<std::string> sample_names;
public:
    VCFWriter(const ThreadingInstructions& instructions);

    void set_chrom(const std::vector<std::string>& _chrom) { chrom = _chrom; }
    void set_pos(const std::vector<std::string>& _pos) { pos = _pos; }
    void set_id(const std::vector<std::string>& _id) { id = _id; }
    void set_ref(const std::vector<std::string>& _ref) { ref = _ref; }
    void set_alt(const std::vector<std::string>& _alt) { alt = _alt; }
    void set_qual(const std::vector<std::string>& _qual) { qual = _qual; }
    void set_filter(const std::vector<std::string>& _filter) { filter = _filter; }
    void set_sample_names(const std::vector<std::string>& _sample_names) { sample_names = _sample_names; }

    void write_vcf();
    void write_header();
};

#endif // THREADS_ARG_VCF_WRITER_HPP
