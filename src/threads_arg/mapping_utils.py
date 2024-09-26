# This file is part of the Threads software suite.
# Copyright (C) 2024 Threads Developers.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import arg_needle_lib
import numpy as np
import logging
import time
import importlib
import os
from cyvcf2 import VCF

def mapping_string(carrier_sets, edges):
    if len(edges) == 0:
        return "NaN"
    elif len(edges) == 1:
        return f"-1,{edges[0].child.height:.4f},{edges[0].parent.height:.4f}"
    else:
        return ";".join([f"{'.'.join([str(c) for c in carrier_set])},{edge.child.height:.4f},{edge.parent.height:.4f}" for carrier_set, edge in zip(carrier_sets, edges)])

def get_leaves(arg, edge, position):
    leaves = []
    populate_leaves(arg, edge, position - arg.offset, leaves)
    return leaves

def populate_leaves(arg, edge, position, leaf_list):
    child = edge.child
    if arg.is_leaf(child.ID):
        return leaf_list.append(child.ID)
    else:
        for edge in child.child_edges_at(position):
            populate_leaves(arg, edge, position, leaf_list)
