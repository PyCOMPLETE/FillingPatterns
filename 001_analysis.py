import json
import matplotlib.pyplot as plt
import numpy as np

import mystyle as ms
import filling_pattern as fp

study_name = 'comparison_25ns'
scheme_fnames = [
'25ns_2760b_2748_2495_2560_288bpi_14inj_800ns_bs200ns_STD.json',
'25ns_2844b_2832_2560_2631_288bpi_15inj_800ns_bs200ns_4x72_opt.json',
'25ns_2808b_2800_2618_2658_320bpi_14inj_800ns_bs200ns_4x80.json',
'25ns_2904b_2896_2656_2734_320bpi_12inj_800ns_bs200ns_4x80b_opt.json',
'25ns_2744b_2736_2242_2370_240bpi_13inj_800ns_bs200ns_5x48.json',
'25ns_2748b_2736_2258_2378_288bpi_12inj_800ns_bs200ns_6x48.json',
]

plt.close('all')
ms.mystyle_arial(fontsz=14, dist_tick_lab=4)
fig1 = plt.figure(1, figsize=(1.5*8, .9*6))
fig1.set_facecolor('w')

with open(study_name+'.tsv', 'w') as fid:
    fid.write('\t'.join([
        'Scheme name',
        'N_bunches',
        'N_coll_ATLAS',
        'N_coll_LHCb',
        'N_coll_ALICE', 
        'N_inj',
        'N_unused',
        'Injection types',
        'Injection len',
        ])+'\n')

for fname in scheme_fnames:
    with open(fname, 'r') as fid:
        data = json.load(fid)

    patt = fp.Filling_Pattern(pattern_b1 = data['beam1'], 
                    pattern_b2 = data['beam2'])

    print('\n____________________\n')
    print('Scheme name: '+fname.split('.')[0])

    print('\n   N. collisions')
    print('ATLAS/CMS:     %d'%patt.n_coll_ATLAS)
    print('LHCb:          %d'%patt.n_coll_LHCb)
    print('ALICE:         %d'%patt.n_coll_ALICE)

    print('\nN. bunches:    %d'%patt.b1.n_bunches)
    print('\nN. injections:  %d'%(patt.b1.n_injections))
    print('\nUnused:   %d slots'%(patt.b1.n_unused_slots))
    print('          (%.1f'%patt.b1.inefficiency_perc+' % LHC)')
    print('\nInjection types = %s'%repr(map(list, patt.b1.inj_composition_types)))
    print('Injection len: %s'%repr(map(len, patt.b1.inj_pattern_types)))

    with open(study_name+'.tsv', 'a') as fid:
        fid.write('\t'.join([
            fname.split('.')[0],
            '%d'%patt.b1.n_bunches,
            '%d'%patt.n_coll_ATLAS,
            '%d'%patt.n_coll_LHCb,
            '%d'%patt.n_coll_ALICE, 
            '%d'%(patt.b1.n_injections),
            '%d'%(patt.b1.n_unused_slots),
            repr(map(list, patt.b1.inj_composition_types)),
            repr(map(len, patt.b1.inj_pattern_types)),
            ])+'\n')

    fig1.clear()
    ax = plt.subplot2grid(shape=(3,2), loc=(0,0), 
            rowspan=1, colspan=2, fig=fig1)
    
    axinj_list = [
        plt.subplot2grid(shape=(3,2), loc=(1,0), 
            rowspan=1, colspan=1, fig=fig1),
        plt.subplot2grid(shape=(3,2), loc=(2,0), 
            rowspan=1, colspan=1, fig=fig1),
        plt.subplot2grid(shape=(3,2), loc=(1,1), 
            rowspan=1, colspan=1, fig=fig1),
        plt.subplot2grid(shape=(3,2), loc=(2,1), 
            rowspan=1, colspan=1, fig=fig1),
        ]

    x_pattern = np.arange(0, len(patt.b1.pattern), 0.1)
    y_pattern = 0.*x_pattern
    for ii, filled in enumerate(patt.b1.pattern):
        if filled>0.5:
            y_pattern[np.logical_and(x_pattern>=ii, x_pattern<(ii+1))] = 1.

    ax.fill_between(x = x_pattern, y1=y_pattern, alpha=0.5, linewidth=0., edgecolor='k')
    ax.set_ylim(0, 1.2)
    ax.set_yticks([])
    ax.set_xlim(0, 3564)
    ax.set_xlabel('25 ns slot')

    fig1.suptitle(fname.split('.')[0])
    fig1.subplots_adjust(left=.03, right=.97, bottom=.12, top=.9, hspace=.8)

    for i_plot, i_inj in enumerate(np.argsort(patt.b1.inj_nbun_types)[::-1]):
        
        this_patt = patt.b1.inj_pattern_types[i_inj]
        thisax = axinj_list[i_plot]

        x_pattern = np.arange(0, len(this_patt), 0.1)
        y_pattern = 0.*x_pattern
        for ii, filled in enumerate(this_patt):
            if filled>0.5:
                y_pattern[np.logical_and(x_pattern>=ii, 
                            x_pattern<(ii+1))] = 1.

        thisax.fill_between(x = x_pattern, y1=y_pattern, 
                    alpha=0.5, linewidth=0., edgecolor='k')
        thisax.set_ylim(0, 1.2)
        thisax.set_yticks([])
        thisax.set_xlim(-10, 350)
        thisax.set_xlabel('25 ns slot')

    fig1.savefig(fname.split('.')[0]+'.png', dpi=200)


plt.show()