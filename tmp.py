from time import time
from threads_arg.serialization import instructions_to_weft, load_weft, load_instructions

instructions = load_instructions("test/data/arg_consistent.threads")
# instructions_to_weft(instructions, "tmp.weft")
start_time = time()
shuttle = load_weft("tmp.weft")
end_time = time()

num_sites = len(instructions.positions)
import numpy as np
rng = np.random.default_rng()
U = rng.random((num_sites, 1000))
out = np.zeros((500, 1000), order="F").astype(np.float32)

start_time = time()
shuttle.GU_mult(U, out)
end_time = time()

shuttle_time = end_time - start_time

import pgenlib
alleles_out = np.empty((num_sites, 2 * 500), dtype=np.int32)
phasepresent_out = np.empty((num_sites, 500), dtype=np.uint8)
reader = pgenlib.PgenReader("test/data/panel.pgen".encode())
reader.read_alleles_and_phasepresent_range(0, num_sites, alleles_out, phasepresent_out)
atrans = alleles_out.transpose()

start_time = time()
res2_tmp = np.dot(atrans, U)
res2 = res2_tmp[::2] + res2_tmp[1::2]
end_time = time()

dot_time = end_time - start_time

assert np.isclose(res2, out).all()
print(dot_time)
print(shuttle_time)
