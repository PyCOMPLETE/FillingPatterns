import fillingpatterns as fp

fname = '25ns_2760b_2748_2492_2574_288bpi_13inj_800ns_bs200ns.csv'
# fname = '25ns_2744b_2736_2246_2370_240bpi_13inj_800ns_bs200ns_BCMS_5x48b.csv'
# fname = '8b4e_1972b_1960_1178_1886_224bpi_12inj_800ns_bs200ns.csv'

# Load filling pattern
patt = fp.FillingPattern.from_json(fname)
# Compute bb schedule
patt.compute_beam_beam_schedule(n_lr_per_side=16, use_old_algorithm = 0)

#########
# Plots #
#########
import matplotlib.pyplot as plt
plt.close('all')

for ib, beam in zip([1, 2], [patt.b1, patt.b2]):

    bbs = beam.bb_schedule

    fig1 = plt.figure(100 + ib, figsize=(6.4*1.5, 1.6*4.8))
    ax1 = fig1.add_subplot(4, 1, 1)
    ax2 = fig1.add_subplot(4, 1, 2, sharex=ax1)
    ax3 = fig1.add_subplot(4, 1, 3, sharex=ax1)
    ax4 = fig1.add_subplot(4, 1, 4, sharex=ax1)

    ax1.plot(bbs['collides in ATLAS/CMS'], '.', color='b', label='ATLAS/CMS')
    ax1.plot(0.05 +bbs['collides in LHCB'], '.', color='r', label='LHCB')
    ax1.plot(-0.05 + bbs['collides in ALICE'],'.', color='g', label='ALICE')
    ax1.legend(ncol=3, loc='center right', fontsize='medium')

    ax1.plot(bbs['collides in ATLAS/CMS'], '.', color='b')
    ax1.plot(0.05 +bbs['collides in LHCB'], '.', color='r')
    ax1.plot(-0.05 + bbs['collides in ALICE'],'.', color='g')

    ax2.plot(bbs['# of LR in ATLAS/CMS'], '.', color='b')
    ax3.plot(bbs['# of LR in LHCB'], '.', color='r')
    ax4.plot(bbs['# of LR in ALICE'], '.', color='g')

    ax1.set_ylabel('Head-on')
    ax2.set_ylabel('N. LR in ATLAS/CMS')
    ax3.set_ylabel('N. LR in LHCB')
    ax4.set_ylabel('N. LR in ALICE')

    ax1.set_yticks([0, 1])
    ax1.set_yticklabels(['No', 'Yes'])
    ax4.set_xlim(0, 3500)
    ax4.set_xlabel('25 ns slot')

    for aa in [ax1, ax2, ax3]:
        aa.tick_params(labelbottom=False)

    for aa in [ax1, ax2, ax3, ax4]:
        aa.grid(True, linestyle=':')

    for aa in [ax2, ax3, ax4]:
        aa.set_ylim(bottom=0)

    fig1.subplots_adjust(left=.06, right=.96, top=.92)
    fig1.suptitle(patt.scheme_name + ' - ' + beam.beam_name)

    fig1.savefig(fname.split('.csv')[0] + f"_{beam.beam_name.replace(' ', '')}_bb_summary.png", dpi=200)

plt.show()

