import json
import matplotlib.pyplot as pl
import numpy as np

import mystyle as ms

scheme_fname = '25ns_2880b_2872_2673_2688_320bpi_12inj_800ns_bs200ns_4x80_opt.json'
scheme_fname = 'example_mixed.json'

ring_length_slots = 3564
min_MKI_slots = 31
min_MKP_slots = 7
agap_first_slot = 3443


with open(scheme_fname, 'r') as fid:
	data = json.load(fid)

beam = 'beam2'

agap_length = ring_length_slots - agap_first_slot

pattern = np.array(data[beam])

pl.close('all')
ms.mystyle_arial(fontsz=14, dist_tick_lab=5)
fig1 = pl.figure(1, figsize=(1.5*8, .5*6))
fig1.set_facecolor('w')

ax = fig1.add_subplot(311)

x_pattern = np.arange(0, len(pattern), 0.1)
y_pattern = 0.*x_pattern

for ii, filled in enumerate(pattern):
	if filled>0.5:
		y_pattern[np.logical_and(x_pattern>=ii, x_pattern<(ii+1))] = 1.

ax.fill_between(x = x_pattern, y1=y_pattern, alpha=0.5, linewidth=0., edgecolor='k')
ax.set_ylim(0, 1.2)
ax.set_yticks([])
fig1.subplots_adjust(left=.03, right=.97, bottom=.2, top=.9)
ax.set_xlim(0, 3564)
fig1.suptitle(scheme_fname.split('.')[0])

# Identify trains
ddd = np.diff(pattern)
start_trains = np.where(ddd==1)[0]+1
end_trains = np.where(ddd==-1)[0]+1
if end_trains[0]<start_trains[0]:
	start_trains = np.array([0] + list(start_trains))

# Remove gaps that need to be generated in the PS
tobekept = (start_trains[1:] - end_trains[:-1]) >= min_MKP_slots
start_trains = np.array([start_trains[0]] + list(start_trains[1:][tobekept]))
end_trains =  np.array(list(end_trains[:-1][tobekept]) + [end_trains[-1]])

# Lenghts of trains and gaps
train_lengths = end_trains - start_trains
gap_lengths =  start_trains[1:] - end_trains[:-1]
train_nbunches = np.array([np.sum(pattern[start_trains[ii]:end_trains[ii]]) for ii in \
					range(len(start_trains))])

i_train_first = np.array([0] + list(np.where(gap_lengths>=min_MKI_slots)[0]+1))
inj_slots = start_trains[i_train_first]
n_injections = len(inj_slots)

#Find pattern for individual injections
inj_patterns = []
inj_compositions = []
for ii, _ in enumerate(i_train_first):

	i_start = inj_slots[ii]

	if ii+1 == n_injections:
		inj_composition = train_nbunches[i_train_first[ii]:-1]
		i_end = -1
	else:
		inj_composition = train_nbunches[i_train_first[ii]:i_train_first[ii+1]]
		i_end = inj_slots[ii+1]

	inj_pattern = pattern[i_start:i_end]
	# Cut away zeros at the end
	inj_last_filled = np.max(np.where(inj_pattern==1)[0])
	inj_pattern = inj_pattern[:inj_last_filled+1]

	inj_patterns.append(inj_pattern)
	inj_compositions.append(inj_composition)

# Identify injection types
inj_pattern_types = []
inj_composition_types = []
for patt, comp in zip(inj_patterns, inj_compositions):
	found = False
	for patt_type in inj_pattern_types:
		if len(patt_type) == len(patt):	
			if np.sum(patt_type-patt)==0:
				found = True
				break
	if not found:
		inj_pattern_types.append(patt)
		inj_composition_types.append(comp)

# Unused space
needed_nslots = (n_injections-1)*min_MKI_slots \
				+ np.sum(map(len, inj_patterns)) + agap_length
inefficiency_perc = 100*(1 - float(needed_nslots) / float(ring_length_slots))



# ax.set_xlim(0, n_part)
# ax.set_xticks(np.arange(*(ax.get_xlim()), step=200))

# fig1.savefig(scheme_fname.split('.')[0]+'_start.png', dpi=200)

# ax.set_xlim(3500-n_part, 3500)
# ax.set_xticks(np.arange(*(ax.get_xlim()), step=200))

# fig1.savefig(scheme_fname.split('.')[0]+'_end.png', dpi=200)

# ax.plot(np.arange(len(pattern))+0.5, pattern, '.')

# ax.set_xlim(100, 250)
# ax.set_xticks(np.arange(*(ax.get_xlim()), step=10))
# fig1.savefig(scheme_fname.split('.')[0]+'_zoom.png', dpi=200)

pl.show()