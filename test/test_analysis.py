
import numpy as np
import sys
import os
print(os.getcwd())
file_list = os.listdir()

# Print the list
for item in file_list:
    print(item)
import bbFunctions as bbF


def test_bool0():
    [df_B1,df_B2] = bbF.bbschedule(np.zeros(3564),np.zeros(3564),20)
    assert df_B1.shape[0] == 0, "the data_frame_B1 should be empty"
    assert df_B2.shape[0] == 0, "the data_frame_B2 should be empty"

def test_bool1():
    v = bbF.events_in_slots([0,1,0],[0,0,1],1)[0]
    assert v.iloc[0].name == 1, "wrong to which bunch referred the collision"
    assert v.iloc[0][0] == 2, "wrong bb partner"

def test_bool2():
    [df_B1,df_B2] = bbF.bbschedule(np.ones(3564),np.ones(3564),1)
    assert np.array([~np.isnan(i) for i in df_B1['HO partner in ATLAS/CMS']]).all(), "the data_frame_B1 should have all H-O collision in ATLAS/CMS"
    assert np.array([~np.isnan(i) for i in df_B2['HO partner in ATLAS/CMS']]).all(), "the data_frame_B2 should have all H-O collision in ATLAS/CMS"
    assert np.array([~np.isnan(i) for i in df_B1['HO partner in ALICE']]).all(), "the data_frame_B1 should have all H-O collision in ALICE"
    assert np.array([~np.isnan(i) for i in df_B2['HO partner in ALICE']]).all(), "the data_frame_B2 should have all H-O collision in ALICE"
    assert np.array([~np.isnan(i) for i in df_B1['HO partner in LHCB']]).all(), "the data_frame_B1 should have all H-O collision in LHCB"
    assert np.array([~np.isnan(i) for i in df_B2['HO partner in LHCB']]).all(), "the data_frame_B2 should have all H-O collision in LHCB"
    assert np.array([df_B1['BB partners in ATLAS/CMS'][i] == [(i-1)%3564,(i+1)%3564] for i in range(3564)]).all(), "the data_frame_B1 should have 2 LR for each bunch in ATLAS/CMS"
    assert np.array([df_B2['BB partners in ATLAS/CMS'][i] == [(i-1)%3564,(i+1)%3564] for i in range(3564)]).all(), "the data_frame_B2 should have 2 LR for each bunch in ATLAS/CMS"
    assert np.array([df_B1['BB partners in ALICE'][i] == [(i+891-1)%3564,(i+891+1)%3564] for i in range(3564)]).all(), "the data_frame_B1 should have 2 LR for each bunch in ALICE"
    assert np.array([df_B2['BB partners in ALICE'][i] == [(i-891-1)%3564,(i-891+1)%3564] for i in range(3564)]).all(), "the data_frame_B2 should have 2 LR for each bunch in ALICE"
    assert np.array([df_B1['BB partners in LHCB'][i] == [(i-894-1)%3564,(i-894+1)%3564] for i in range(3564)]).all(), "the data_frame_B1 should have 2 LR for each bunch in LHCB"
    assert np.array([df_B2['BB partners in LHCB'][i] == [(i+894-1)%3564,(i+894+1)%3564] for i in range(3564)]).all(), "the data_frame_B2 should have 2 LR for each bunch in LHCB"



