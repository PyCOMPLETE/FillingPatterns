
import numpy as np


class Filling_Pattern_Single_Beam(object):

    def __init__(self, pattern, ring_length_slots = 3564, 
            min_MKI_slots = 31, min_MKP_slots = 7, agap_first_slot = 3443):

        agap_length = ring_length_slots - agap_first_slot
        pattern = np.array(pattern)

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
                        + np.sum(map(len, inj_patterns)) + agap_length
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
        self.inj_nbun = np.array(map(np.sum, inj_compositions))

        self.inj_pattern_types = inj_pattern_types
        self.inj_composition_types = inj_composition_types
        self.inj_nbun_types = np.array(map(np.sum, inj_composition_types))

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



class Filling_Pattern(object):

    def __init__(self, pattern_b1, pattern_b2, ring_length_slots = 3564, 
            min_MKI_slots = 31, min_MKP_slots = 7, agap_first_slot = 3443):

        self.b1 = Filling_Pattern_Single_Beam(pattern_b1, ring_length_slots,
            min_MKI_slots, min_MKP_slots, agap_first_slot)

        self.b2 = Filling_Pattern_Single_Beam(pattern_b2, ring_length_slots,
            min_MKI_slots, min_MKP_slots, agap_first_slot)

        self.n_coll_ATLAS = np.sum(self.b1.pattern*self.b2.pattern)
        self.n_coll_CMS = self.n_coll_ATLAS
        self.n_coll_ALICE = np.sum(self.b2.pattern * \
                                np.roll(self.b1.pattern, 891))
        self.n_coll_LHCb = np.sum(self.b1.pattern * \
                                np.roll(self.b2.pattern, 894))
