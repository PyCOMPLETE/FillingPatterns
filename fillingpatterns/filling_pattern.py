import json

import numpy as np


class FillingPatternSingleBeam(object):

    def __init__(self, pattern, ring_length_slots = 3564,
            min_MKI_slots = 31, min_MKP_slots = 7, agap_first_slot = 3443):

        agap_length = ring_length_slots - agap_first_slot
        pattern = np.int_(np.array(pattern))

        assert(len(pattern)<=ring_length_slots)

        # Add empty slots at the end if needed
        pattern = np.array(list(pattern) + (ring_length_slots-len(pattern))*[0])

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
        injection_number = 0*pattern - 1
        slot_within_injection = 0*pattern - 1
        for ii, _ in enumerate(i_train_first):

            i_start = inj_slots[ii]

            if ii+1 == n_injections:
                inj_composition = train_nbunches[i_train_first[ii]:]
                inj_pattern = pattern[i_start:]
            else:
                inj_composition = train_nbunches[i_train_first[ii]:i_train_first[ii+1]]
                i_end = inj_slots[ii+1]
                inj_pattern = pattern[i_start:i_end]

            # Cut away zeros at the end
            inj_last_filled = np.max(np.where(inj_pattern==1)[0])
            inj_pattern = inj_pattern[:inj_last_filled+1]
            inj_len = len(inj_pattern)
            injection_number[i_start:i_start+inj_len] = ii
            slot_within_injection[i_start:i_start+inj_len] = np.arange(inj_len)

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
                        + np.sum(list(map(len, inj_patterns))) + agap_length
        inefficiency_perc = 100*(1 - float(needed_nslots) / float(ring_length_slots))

        self.actual_MKP_slots = min(gap_lengths)
        self.actual_MKI_slots = min(gap_lengths[gap_lengths>min_MKP_slots])

        self.agap_length = agap_length

        self.pattern = pattern
        self.n_bunches = np.sum(pattern)

        self.train_lengths =  train_lengths
        self.gap_lengths = gap_lengths
        self.train_nbunches = train_nbunches

        self.inj_slots = inj_slots
        self.n_injections = n_injections

        self.inj_patterns = inj_patterns
        self.inj_compositions = inj_compositions
        self.inj_nbun = np.array(list(map(np.sum, inj_compositions)))

        self.inj_pattern_types = inj_pattern_types
        self.inj_composition_types = inj_composition_types
        self.inj_nbun_types = np.array(list(map(np.sum, inj_composition_types)))

        self.n_unused_slots = ring_length_slots - needed_nslots
        self.inefficiency_perc = inefficiency_perc

        self.injection_number = injection_number
        self.slot_within_injection = slot_within_injection

        self.ring_length_slots = ring_length_slots
        self.min_MKI_slots = min_MKI_slots
        self.min_MKP_slots = min_MKP_slots
        self.agap_first_slot = agap_first_slot

    def belongs_to_group(self, slots_within_injection=None,
            n_bunches_in_injection=None):

        mask = self.pattern>0

        # Check for nbun in injection
        if n_bunches_in_injection is not None:
            nbun_my_inj = np.take(self.inj_nbun, self.injection_number)
            nbun_my_inj[self.injection_number==-1]=-1
            mask = np.logical_and(mask, nbun_my_inj == n_bunches_in_injection)

        if slots_within_injection is not None:
            mask = np.logical_and(mask, np.array(
                [(ss in slots_within_injection) for ss in self.slot_within_injection]))

        return mask


class FillingPattern(object):

    @classmethod
    def from_json(cls, fname):
        with open(fname, 'r') as fid:
            data = json.load(fid)
        patt = cls(pattern_b1 = data['beam1'],
                   pattern_b2 = data['beam2'],
                   scheme_name=fname.split('.json')[0])
        return patt

    @classmethod
    def from_csv(cls, fname):

        with open(fname, 'r') as fid:
            lines = fid.readlines()

        patt_dict = {ll.split(',', 1)[0]: ll.split(',', 1)[1].replace('\n', '') for ll in lines}

        for bb in [1, 2]:
            # Build slot array
            patt_arr = np.array(3564*[False])
            filled_buckets = np.fromstring(patt_dict[f'B{bb} filled buckets'],
                       sep=',', dtype=np.int)
            for fbk in filled_buckets:
                assert(np.mod(fbk-1, 10) == 0)
                slot = (fbk-1)//10
                patt_arr[slot] = True
            patt_dict[f'pattern_array_b{bb}'] = patt_arr

            # Get gaps in slots
            mki_gap_in_slots = int(np.round(
                float(patt_dict[f'B{bb} LHC injection kicker gap [ns]'])/25 + - 1))
            patt_dict[f'mki_gap_in_slots_b{bb}'] = mki_gap_in_slots
            mkp_gap_in_slots = int(np.round(
                float(patt_dict[f'B{bb} SPS injection kicker gap [ns]'])/25 + - 1))
            patt_dict[f'mkp_gap_in_slots_b{bb}'] = mkp_gap_in_slots
            abort_gap_in_slots = int(np.round(
                float(patt_dict[f'B{bb} abort gap [ns]'])/25 + - 1))
            patt_dict[f'abort_gap_in_slots_b{bb}'] = abort_gap_in_slots

        patt_dict['mki_min'] = np.min([patt_dict[f'mki_gap_in_slots_b{bb}'] for bb in [1, 2]])
        patt_dict['mkp_min']= np.min([patt_dict[f'mkp_gap_in_slots_b{bb}'] for bb in [1, 2]])
        patt_dict['agap_min'] = np.min([patt_dict[f'abort_gap_in_slots_b{bb}'] for bb in [1, 2]])
        patt_dict['agap_first_slot'] = 3564 - patt_dict['agap_min']

        patt = cls(pattern_b1=patt_dict['pattern_array_b1'],
            pattern_b2=patt_dict['pattern_array_b2'],
            ring_length_slots = 3564,
            min_MKI_slots=patt_dict['mki_min'],
            min_MKP_slots=patt_dict['mkp_min'],
            agap_first_slot=patt_dict['agap_first_slot'],
            scheme_name=patt_dict['Scheme name'])

        return patt

    def __init__(self, pattern_b1, pattern_b2, ring_length_slots = 3564,
            min_MKI_slots = 31, min_MKP_slots = 7, agap_first_slot = 3443,
            scheme_name='no_name'):

        self.scheme_name = scheme_name
        self.b1 = FillingPatternSingleBeam(pattern_b1, ring_length_slots,
            min_MKI_slots, min_MKP_slots, agap_first_slot)

        self.b2 = FillingPatternSingleBeam(pattern_b2, ring_length_slots,
            min_MKI_slots, min_MKP_slots, agap_first_slot)

        self.b1.beam_name = 'beam 1'
        self.b2.beam_name = 'beam 2'

        self.n_coll_ATLAS = np.sum(self.b1.pattern*self.b2.pattern)
        self.n_coll_CMS = self.n_coll_ATLAS
        self.n_coll_ALICE = np.sum(self.b2.pattern * \
                                np.roll(self.b1.pattern, 891))
        self.n_coll_LHCb = np.sum(self.b1.pattern * \
                                np.roll(self.b2.pattern, 894))

    def to_csv(self, fname='auto'):

        txtlines = []
        txtlines.append('Scheme name,' + self.scheme_name)
        txtlines.append(','.join(['B1 filled buckets'] +
                    ['%d'%(bb*10+1) for bb in np.where(self.b1.pattern)[0]]))
        txtlines.append(','.join(['B2 filled buckets'] +
                    ['%d'%(bb*10+1) for bb in np.where(self.b2.pattern)[0]]))
        txtlines.append(','.join(['B1 injection buckets'] +
                    ['%d'%(bb*10+1) for bb in self.b1.inj_slots]))
        txtlines.append(','.join(['B2 injection buckets'] +
                    ['%d'%(bb*10+1) for bb in self.b2.inj_slots]))
        txtlines.append(','.join(['B1 injection n. bunches'] +
                    ['%d'%bb for bb in self.b1.inj_nbun]))
        txtlines.append(','.join(['B2 injection n. bunches'] +
                    ['%d'%bb for bb in self.b2.inj_nbun]))
        txtlines.append('N. coll. ATLAS/CMS, %d'%(self.n_coll_ATLAS))
        txtlines.append('N. coll. LHCb, %d'%(self.n_coll_LHCb))
        txtlines.append('N. coll. ALICE, %d'%(self.n_coll_ALICE))
        txtlines.append('B1 n. bunches, %d'%(self.b1.n_bunches))
        txtlines.append('B2 n. bunches, %d'%(self.b2.n_bunches))
        txtlines.append('B1 n. injections, %d'%(self.b1.n_injections))
        txtlines.append('B2 n. injections, %d'%(self.b2.n_injections))
        txtlines.append(','.join(['B1 max. injection length [ns]',
                    '%d'%((max(map(len, self.b1.inj_pattern_types)))*25)]))
        txtlines.append(','.join(['B2 max. injection length [ns]',
                    '%d'%((max(map(len, self.b2.inj_pattern_types)))*25)]))
        txtlines.append(','.join(['B1 LHC injection kicker gap [ns]',
                    '%d'%((self.b1.actual_MKI_slots+1)*25)]))
        txtlines.append(','.join(['B2 LHC injection kicker gap [ns]',
                    '%d'%((self.b2.actual_MKI_slots+1)*25)]))
        txtlines.append(','.join(['B1 SPS injection kicker gap [ns]',
                    '%d'%((self.b1.actual_MKP_slots+1)*25)]))
        txtlines.append(','.join(['B2 SPS injection kicker gap [ns]',
                    '%d'%((self.b2.actual_MKP_slots+1)*25)]))
        txtlines.append(','.join(['B1 abort gap [ns]',
                    '%d'%((self.b1.agap_length+1)*25)]))
        txtlines.append(','.join(['B2 abort gap [ns]',
                    '%d'%((self.b2.agap_length+1)*25)]))

        if fname == 'auto':
            fnamecsv = self.scheme_name + '.csv'
        else:
            fnamecsv = fname

        with open(fnamecsv, 'w') as fcsv:
            fcsv.write('\n'.join(txtlines))

    def compute_beam_beam_schedule(self, n_lr_per_side):
        from . import bbFunctions
        print('Computing collision schedules...')
        self.b1.bb_schedule = bbFunctions.B1CollisionScheduleDF(
            self.b1.pattern, self.b2.pattern,  n_lr_per_side)
        print('Done Beam 1')
        self.b2.bb_schedule = bbFunctions.B2CollisionScheduleDF(
            self.b1.pattern, self.b2.pattern,  n_lr_per_side)
        print('Done Beam 2')

        for bb in [self.b1, self.b2]:
            for ee in ['ATLAS/CMS', 'ALICE', 'LHCB']:
                 bb.bb_schedule[f'collides in {ee}'] = ~(np.isnan(bb.bb_schedule[f'HO partner in {ee}']))

