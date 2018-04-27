import sys
import os
import numpy as np
sys.path.insert(0, os.path.abspath(".."))
from analysis import misclib


hiv_shape = "hiv_data/hiv_rre/hiv_rre.shape"

shape = misclib.read_n_shape(hiv_shape, 2)

mix1 = []
for sl4, sl5 in zip(shape["4SL"], shape["5SL"]):

    if sl4 == -999 or sl5 == -999:
        mix1.append(-999.0)
        continue

    if sl4 == np.nan or sl4 == np.nan:
        mix1.append('nan')
        continue

    mix1.append(0.5*sl4 + 0.5*sl5) # Manually enter desired weight fractions

vals = " ".join([str(x) for x in mix1])
print(vals)
