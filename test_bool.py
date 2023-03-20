#%%

import fillingpatterns as fp
import numpy as np
from fillingpatterns import bbFunctions as bbF


# %%
def test_bool0():
    [df_B1,df_B2] = fp.FillingPattern.compute_beam_beam_schedule(np.zeros(3564),np.zeros(3564))
    assert df_B1.shape[0] == 0, "the data_frame_B1 should be empty"
    assert df_B2.shape[0] == 0, "the data_frame_B2 should be empty"

#%%
def test_trivial_bb():
    df = bbF.CollisionsInTheGivenSlot_vec([1,0],[1,0],0,0)[0]
    assert df["BB_collision"] == 0, "should collide with the first slot"
    assert df["pos_collision"] == 0, "should collide in slot 0"
    
    


#%%