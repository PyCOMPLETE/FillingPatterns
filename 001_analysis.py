import json
import matplotlib.pyplot as pl
import numpy as np

import mystyle as ms
import filling_pattern as fp

scheme_fname = '25ns_2880b_2872_2673_2688_320bpi_12inj_800ns_bs200ns_4x80_opt.json'

with open(scheme_fname, 'r') as fid:
	data = json.load(fid)

scheme = fp.Filling_Pattern(pattern_b1 = data['beam1'], 
				pattern_b2 = data['beam2'])