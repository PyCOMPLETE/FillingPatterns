import json
import matplotlib.pyplot as pl
import numpy as np

import mystyle as ms
import filling_pattern as fp

scheme_fnames = [
'25ns_2760b_2748_2495_2560_288bpi_14inj_800ns_bs200ns_STD.json',
'25ns_2844b_2832_2560_2631_288bpi_15inj_800ns_bs200ns_4x72_opt.json',
'25ns_2808b_2800_2618_2658_320bpi_14inj_800ns_bs200ns_4x80.json',
'25ns_2880b_2872_2673_2688_320bpi_12inj_800ns_bs200ns_4x80_opt.json',
'25ns_2744b_2736_2242_2370_240bpi_13inj_800ns_bs200ns_5x48.json',
'25ns_2748b_2736_2258_2378_288bpi_12inj_800ns_bs200ns_6x48.json',
]

for fname in scheme_fnames:
    with open(fname, 'r') as fid:
        data = json.load(fid)

    patt = fp.Filling_Pattern(pattern_b1 = data['beam1'], 
                    pattern_b2 = data['beam2'])

    print('\n____________________\n')
    print(fname.split('.')[0])
    print('N_bunches = \t%d'%patt.b1.n_bunches)
    print('N_coll_ATLAS = \t%d'%patt.n_coll_ATLAS)
    print('N_coll_LHCb = \t%d'%patt.n_coll_LHCb)
    print('N_coll_ALICE = \t%d'%patt.n_coll_ALICE)
    print('N_inj = \t%d'%(patt.b1.n_injections))
    print('N_unused = \t%d (%.1f'%(patt.b1.n_unused_slots,
            patt.b1.inefficiency_perc)+' %)')
    print('Injection types = %s'%repr(map(list, patt.b1.inj_composition_types)))