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
'25ns_2604b_2592_2288_2396_288bpi_12inj_900ns_bs225ns_BCMS_baseline.json',
'25ns_2604b_2592_2313_2416_240bpi_13inj_800ns_bs200ns_5x48_LHCb.json',
'25ns_2744b_2736_2246_2370_240bpi_13inj_800ns_bs200ns_5x48b_opt.json',
'25ns_2748b_2736_2258_2378_288bpi_12inj_800ns_bs200ns_6x48.json',
'25ns_2372b_2360_1784_2216_256bpi_12inj_800ns_bs200ns_run3study.json',
'25ns_2492b_2480_2048_2301_240bpi_13inj_800ns_bs200ns_run3_study_corrected.json'
]

study_name = 'hl_lhc_scenarios'
scheme_fnames = [
 '25ns_2760b_2748_2492_2574_288bpi_13inj_800ns_bs200ns.json',
 '25ns_2744b_2736_2246_2370_240bpi_13inj_800ns_bs200ns_BCMS_5x48b.json',
 '8b4e_1972b_1960_1178_1886_224bpi_12inj_800ns_bs200ns.json'
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

for ifname, fname in enumerate(scheme_fnames):
    with open(fname, 'r') as fid:
        data = json.load(fid)

    patt = fp.Filling_Pattern(pattern_b1 = data['beam1'], 
                    pattern_b2 = data['beam2'])

    if ifname == 0:
        ref_ATLAS = patt.n_coll_ATLAS
        ref_LHCb = patt.n_coll_LHCb
        ref_ALICE = patt.n_coll_ALICE 


    print('\n____________________\n')
    print('Scheme name: '+fname.split('.')[0])

    print('\n   N. collisions')
    print('ATLAS/CMS:     %d (%+.1f'%(patt.n_coll_ATLAS, 100*(float(patt.n_coll_ATLAS)/ref_ATLAS-1.))+'%)')
    print('LHCb:          %d (%+.1f'%(patt.n_coll_LHCb, 100*(float(patt.n_coll_LHCb)/ref_LHCb-1.))+'%)')
    print('ALICE:         %d (%+.1f'%(patt.n_coll_ALICE, 100*(float(patt.n_coll_ALICE)/ref_ALICE-1.))+'%)')

    print('\nN. bunches:    %d'%patt.b1.n_bunches)
    print('\nN. injections: %d'%(patt.b1.n_injections))
    print('\nUnused:        %d slots'%(patt.b1.n_unused_slots))
    print('               (%.1f'%patt.b1.inefficiency_perc+'% LHC)')
    print('\nPatterns from SPS:')
    for pp in patt.b1.inj_composition_types:
        print('%s'%repr(list(pp)))
    print('\nInjection len: %s'%repr(map(len, patt.b1.inj_pattern_types)))
    print('\nGaps = %s'%repr(list(set(patt.b1.gap_lengths))))
    print('Abort gap: %s'%repr(patt.b1.agap_length))

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

        if i_plot == 3:
            break

    fig1.savefig(fname.split('.')[0]+'.png', dpi=200)


    txtlines = []
    txtlines.append('Scheme name, '+fname.split('.')[0])
    txtlines.append(','.join(['B1 filled buckets'] +
                ['%d'%(bb*10+1) for bb in np.where(patt.b1.pattern)[0]]))
    txtlines.append(','.join(['B2 filled buckets'] +
                ['%d'%(bb*10+1) for bb in np.where(patt.b2.pattern)[0]]))
    txtlines.append(','.join(['B1 injection buckets'] +
                ['%d'%(bb*10+1) for bb in patt.b1.inj_slots]))
    txtlines.append(','.join(['B2 injection buckets'] +
                ['%d'%(bb*10+1) for bb in patt.b2.inj_slots]))
    txtlines.append(','.join(['B1 injection n. bunches'] +
                ['%d'%bb for bb in patt.b1.inj_nbun]))
    txtlines.append(','.join(['B2 injection n. bunches'] +
                ['%d'%bb for bb in patt.b2.inj_nbun]))
    txtlines.append('N. coll. ATLAS/CMS, %d'%(patt.n_coll_ATLAS))
    txtlines.append('N. coll. LHCb, %d'%(patt.n_coll_LHCb))
    txtlines.append('N. coll. ALICE, %d'%(patt.n_coll_ALICE))
    txtlines.append('B1 n. bunches, %d'%(patt.b1.n_bunches))
    txtlines.append('B2 n. bunches, %d'%(patt.b2.n_bunches))
    txtlines.append('B1 n. injections, %d'%(patt.b1.n_injections))
    txtlines.append('B2 n. injections, %d'%(patt.b2.n_injections))
    txtlines.append(','.join(['B1 max. injection length [ns]', 
                '%d'%((max(map(len, patt.b1.inj_pattern_types)))*25)]))
    txtlines.append(','.join(['B2 max. injection length [ns]', 
                '%d'%((max(map(len, patt.b2.inj_pattern_types)))*25)])) 
    txtlines.append(','.join(['B1 LHC injection kicker gap [ns]', 
                '%d'%((patt.b1.actual_MKI_slots+1)*25)]))
    txtlines.append(','.join(['B2 LHC injection kicker gap [ns]', 
                '%d'%((patt.b2.actual_MKI_slots+1)*25)]))
    txtlines.append(','.join(['B1 SPS injection kicker gap [ns]', 
                '%d'%((patt.b1.actual_MKP_slots+1)*25)]))
    txtlines.append(','.join(['B2 SPS injection kicker gap [ns]', 
                '%d'%((patt.b2.actual_MKP_slots+1)*25)])) 
    txtlines.append(','.join(['B1 abort gap [ns]', 
                '%d'%((patt.b1.agap_length+1)*25)]))
    txtlines.append(','.join(['B2 abort gap [ns]', 
                '%d'%((patt.b2.agap_length+1)*25)]))

    

    fnamecsv = fname.split('.')[0]+'.csv'
    with open(fnamecsv, 'w') as fcsv:
        fcsv.write('\n'.join(txtlines))

plt.show()