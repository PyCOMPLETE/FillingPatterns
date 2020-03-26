import fillingpatterns as fp
fname = '25ns_2744b_2736_2246_2370_240bpi_13inj_800ns_bs200ns_BCMS_5x48b.csv'
patt = fp.FillingPattern.from_csv(fname)

patt.compute_beam_beam_schedule(n_lr_per_side=16)

# Choose beam
bbs = patt.b1.bb_schedule

import matplotlib.pyplot as plt
plt.close('all')
fig1 = plt.figure(1, figsize=(6.4*1.5, 4.8))
ax1 = fig1.add_subplot(2, 1, 1)
ax2 = fig1.add_subplot(2, 1, 2, sharex=ax1)

ax1.plot(bbs['collides in ATLAS/CMS'], '.', color='b')
ax1.plot(0.05 +bbs['collides in LHCB'], '.', color='r')
ax1.plot(-0.05 + bbs['collides in ALICE'],'.', color='g')

ax2.plot(bbs['# of LR in ATLAS/CMS'], '.', color='b')
ax2.plot(bbs['# of LR in LHCB'], '.', color='r')
ax2.plot(bbs['# of LR in ALICE'], '.', color='g')

ax2.set_xlim(0, 3500)
fig1.subplots_adjust(left=.06, right=.96)
plt.show()
