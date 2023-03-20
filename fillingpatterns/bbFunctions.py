import numpy as np 
from . import dotdict
import pandas as pd
dotdict=dotdict.dotdict


def computeBBMatrix(numberOfLRToConsider):
        """
        It returns a beam-beam matrix, that is a representation of the beam-beam encounters scheduled for a filled machine (3564 bunches).
        To obtain the BB encounters scheduled of the bunch N of B1 (for the machine totally full) you have to consider the N-row (e.g., BBMatrix[N,:]).
        To obtain the BB encounters scheduled of the bunch N of B2 (for the machine totally full) you have to consider the N-column (e.g., BBMatrix[:,N]).
                
        The numberOfLRToConsider represents the long-range number at the right and left of each IR (the total number of LR per IR is therefore 2 x numberOfLRToConsider).
        numberOfLRToConsider can be a scalar of an array of 3 scalars (longer arrays will be truncated to the length of 3).
        In case numberOfLRToConsider is as array the first scalar represent the long-range number at the right and left of ATLAS/CMS. 
        The second and third element of the array will be the number of long range encounters (at the right and left) in IR2 and IR8, respectively.
        
        The matrix element has a value 1, 2, 5 and 8 when there is a HO  in IP1, 2, 5 and 8, respectively.
        The matrix element has a value 10, 20, 50 and  80 when there is a LR respectively in IR1, 2, 5 and 8, respectively.
        
        It assumes that the positions of the IPs and the convention of the B1/B2 bunch numbering are such that:
        
        1. B1 Bunch 0 meets B2 Bunch 0 in IP1 and 5.
        2. B1 Bunch 0 meets B2 Bunch 891 in IP2.
        3. B1 Bunch 0 meets B2 Bunch 2670 in IP8.
        4. B2 Bunch 0 meets B1 Bunch 0 in IP1 and 5.
        5. B2 Bunch 0 meets B1 Bunch 2673 in IP2.
        6. B2 Bunch 0 meets B1 Bunch 894 in IP8.
        === EXAMPLE 1 ===
        myMatrix=computeBBMatrix(numberOfLRToConsider=20)
        #in this case the total number of LR will be 160: 40 in IR1/5, 40 in IR2, and 40 in IR8.
        === EXAMPLE 2 ===
        myMatrix=computeBBMatrix(numberOfLRToConsider=[20,15,16])
        #in this case the total number of LR will be 142: 40 in IR1/5, 30 in IR2, and 32 in IR8.
        """ 
        availableBunchSlot=3564
        BBMatrixLHC =np.zeros([availableBunchSlot,availableBunchSlot]);

        if isinstance(numberOfLRToConsider,int):
            numberOfLRToConsiderATLASCMS=numberOfLRToConsider
            numberOfLRToConsiderALICE=numberOfLRToConsider
            numberOfLRToConsiderLHCB=numberOfLRToConsider
        else:
            numberOfLRToConsiderATLASCMS=numberOfLRToConsider[0]
            numberOfLRToConsiderALICE=numberOfLRToConsider[1]
            numberOfLRToConsiderLHCB=numberOfLRToConsider[2]

        # HO in IP1 and IP5
        index=np.arange(availableBunchSlot)
        BBMatrixLHC[index,index]=1

        # BBLR in IP1 and IP5
        for i in range(1,numberOfLRToConsiderATLASCMS+1):
            index=np.arange(availableBunchSlot-i)
            BBMatrixLHC[index,index+i]=10
            BBMatrixLHC[index,index-i]=10
            BBMatrixLHC[index+i,index]=10
            BBMatrixLHC[index-i,index]=10

        # HO in IP2
        IP2slot=int(availableBunchSlot/4)
        index=np.arange(availableBunchSlot-IP2slot)
        BBMatrixLHC[index,index+IP2slot]=2
        index=np.arange(availableBunchSlot-IP2slot,availableBunchSlot)
        BBMatrixLHC[index,index-(availableBunchSlot-IP2slot)]=2

        # BBLR in IP2
        for i in range(1,numberOfLRToConsiderALICE+1):
            index=np.arange(availableBunchSlot-IP2slot-i)
            BBMatrixLHC[index,index+IP2slot+i]=20
            BBMatrixLHC[index+i,index+IP2slot]=20
            BBMatrixLHC[index,index+IP2slot-i]=20
            BBMatrixLHC[index-i,index+IP2slot]=20
            index=np.arange(availableBunchSlot-IP2slot,availableBunchSlot-i)
            BBMatrixLHC[index,index-(availableBunchSlot-IP2slot)+i]=20
            BBMatrixLHC[index+i,index-(availableBunchSlot-IP2slot)]=20
            BBMatrixLHC[index,index-(availableBunchSlot-IP2slot)-i]=20
            BBMatrixLHC[index-i,index-(availableBunchSlot-IP2slot)]=20

        # HO in IP8
        IP8slot=int(availableBunchSlot/4*3-3)
        index=np.arange(availableBunchSlot-IP8slot)
        BBMatrixLHC[index,index+IP8slot]=8
        index=np.arange(availableBunchSlot-IP8slot,availableBunchSlot)
        BBMatrixLHC[index,index-(availableBunchSlot-IP8slot)]=8

        # BBLR in IP8
        for i in range(1,numberOfLRToConsiderLHCB+1):
            index=np.arange(availableBunchSlot-IP8slot-i)
            BBMatrixLHC[index,index+IP8slot+i]=80
            BBMatrixLHC[index+i,index+IP8slot]=80
            BBMatrixLHC[index,index+IP8slot-i]=80
            BBMatrixLHC[index-i,index+IP8slot]=80
            index=np.arange(availableBunchSlot-IP8slot,availableBunchSlot-i)
            BBMatrixLHC[index,index-(availableBunchSlot-IP8slot)+i]=80
            BBMatrixLHC[index+i,index-(availableBunchSlot-IP8slot)]=80
            BBMatrixLHC[index,index-(availableBunchSlot-IP8slot)-i]=80
            BBMatrixLHC[index-i,index-(availableBunchSlot-IP8slot)]=80

        return BBMatrixLHC;

def _bunch_BB_pattern(Bunch,BBMatrixLHC):
    '''
    It returns the beam-beam pattern of a bunch of B1 and B2 [adimensional array of integer] in ALICE, ATLAS, CMS, LHCB.
    
    - Bunch [adimensional integer]: the bunch number in B1 and B2 to consider.
    - BBMatrixLHC [adimensional integer matrix]: the beam-beam matrix to consider (see bbFunctions.computeBBMatrix?).
    The returned array is ordered with respect to the positive direction of B1 (clockwise in LHC). 
    WARNING: the bunch number is defined wrt the negative direction of each beam.
    
    Conventions :
    
    ...=>B1 bunch 1 => B1 bunch 0 => |[IP]| <= B2 bunch 0 <= B2 bunch 1 <=...
    
    This means that B1 bunch 0 will meet B2 bunch 1 on the positive side of the IP and B2 bunch 0 will meet B1 bunch 0 on the negative side of IP1 (positive/negative wrt the B1).
    
    '''
    
    #### BEAM 1 ####
    BBVector=BBMatrixLHC[Bunch,:]
    
    numberOfLRToConsider=len(np.where(BBMatrixLHC[Bunch,:]==10)[0])/2
    HO_in_IP=BBVector==1
    LR_in_IP=BBVector==10
    aux=np.where((LR_in_IP) | (HO_in_IP))[0]
    np.where(aux==Bunch)[0]
    B1=np.roll(aux,  numberOfLRToConsider-np.where(aux==np.where(HO_in_IP)[0][0])[0])
    resultsB1=dotdict({'atATLAS':B1,
         'atCMS':B1})
    
    HO_in_IP=BBVector==2
    LR_in_IP=BBVector==20
    numberOfLRToConsider=len(np.where(BBMatrixLHC[Bunch,:]==20)[0])/2
    aux=np.where((LR_in_IP) | (HO_in_IP))[0]
    np.where(aux==Bunch)[0]
    B1=np.roll(aux, numberOfLRToConsider-np.where(aux==np.where(HO_in_IP)[0][0])[0])
    resultsB1.update({'atALICE':B1})
    
    HO_in_IP=BBVector==8
    LR_in_IP=BBVector==80
    numberOfLRToConsider=len(np.where(BBMatrixLHC[Bunch,:]==80)[0])/2
    aux=np.where((LR_in_IP) | (HO_in_IP))[0]
    np.where(aux==Bunch)[0]
    B1=np.roll(aux, numberOfLRToConsider-np.where(aux==np.where(HO_in_IP)[0][0])[0])
    resultsB1.update({'atLHCB':B1})
    
    #### BEAM 2 ####
    BBVector=BBMatrixLHC[:,Bunch]
    numberOfLRToConsider=len(np.where(BBMatrixLHC[Bunch,:]==10)[0])/2
    HO_in_IP=BBVector==1
    LR_in_IP=BBVector==10
    aux=np.where((LR_in_IP) | (HO_in_IP))[0]
    np.where(aux==Bunch)[0]
    B2=np.roll(aux,  numberOfLRToConsider-np.where(aux==np.where(HO_in_IP)[0][0])[0])
    resultsB2=dotdict({'atATLAS':B2[::-1], 
         'atCMS':B2[::-1]}) # To note the inverstion of the direction
   
    HO_in_IP=BBVector==2
    LR_in_IP=BBVector==20
    numberOfLRToConsider=len(np.where(BBMatrixLHC[Bunch,:]==20)[0])/2
    aux=np.where((LR_in_IP) | (HO_in_IP))[0]
    np.where(aux==Bunch)[0]
    B2=np.roll(aux, numberOfLRToConsider-np.where(aux==np.where(HO_in_IP)[0][0])[0])
    resultsB2.update({'atALICE':B2[::-1]})
    
    HO_in_IP=BBVector==8
    LR_in_IP=BBVector==80
    numberOfLRToConsider=len(np.where(BBMatrixLHC[Bunch,:]==80)[0])/2
    aux=np.where((LR_in_IP) | (HO_in_IP))[0]
    np.where(aux==Bunch)[0]
    B2=np.roll(aux, numberOfLRToConsider-np.where(aux==np.where(HO_in_IP)[0][0])[0])
    resultsB2.update({'atLHCB':B2[::-1]})
    # The encounters seen by  B1/2 are in increasing/decreasing order
    return dotdict({'atB1':resultsB1, 'atB2':resultsB2})

def BBEncounterSchedule(B1_fillingScheme,B2_fillingScheme,BBMatrixLHC):
    """
    It returns a dictionary structure with the BB encounters of B1 and B2 taking into account the filling schemes.
    - B1_fillingScheme [adimensional integer array]: the B1 filling scheme.
    - B2_fillingScheme [adimensional integer array]: the B2 filling scheme.
    The dictionary structure has the following hierarchy:
    - BEAM >> BUNCH >> EXPERIMENT >> ENCOUNTERS
    - BEAM >> BUNCH >> EXPERIMENT >> POSITIONS
    All the positions are referred to the positive direction of B1 (clockwise in LHC).
    WARNING: the bunch number is defined wrt the negative direction of each beam.
    
    === EXAMPLE 1 ===
    from cl2pd import bbFunctions
    from cl2pd import importData
   
    np=importData.np  
    BBMatrixLHC=bbFunctions.computeBBMatrix(numberOfLRToConsider=25)
    B1_bunches=np.array([0,1,2])
    B2_bunches=np.array([0,1,2])
    results=beam_BB_pattern(B1_bunches, B2_bunches, BBMatrixLHC)
    #or more realistically
    fillingSchemeDF=importData.LHCCals2pd(['LHC.BCTFR.A6R4.B%:BUNCH_FILL_PATTERN'],6972, ['FLATTOP'],flag='next')
    B1_bunches = fillingSchemeDF['LHC.BCTFR.A6R4.B1:BUNCH_FILL_PATTERN'].iloc[0]
    B2_bunches = fillingSchemeDF['LHC.BCTFR.A6R4.B2:BUNCH_FILL_PATTERN'].iloc[0]
    results=bbFunctions.BBEncounterSchedule(np.where(B1_bunches)[0], np.where(B2_bunches)[0], BBMatrixLHC)
    """
    experiments=['atALICE','atATLAS','atCMS','atLHCB']

    #B1
    B1_BB_pattern=dotdict()
    for i in B1_fillingScheme:
        bunch_aux=dotdict()
        for j in experiments:
            results=_bunch_BB_pattern(i,BBMatrixLHC)
            B2=results['atB1'][j]
            aux=B2[np.in1d(B2,B2_fillingScheme)]
            myPosition=np.arange(-(len(B2)-1)/2,(len(B2)-1)/2+1)
            bunch_aux.update({j: {'atEncounters' : aux,'atPositions':myPosition[np.in1d(B2,B2_fillingScheme)]}})
            B1_BB_pattern.update({'at'+format(i,'04d'):bunch_aux})

    #B2
    B2_BB_pattern=dotdict()
    for i in B2_fillingScheme:
        bunch_aux=dotdict()
        for j in experiments:
            results=_bunch_BB_pattern(i,BBMatrixLHC)
            B1=results['atB2'][j]
            aux=B1[np.in1d(B1,B1_fillingScheme)]
            myPosition=np.arange(-(len(B1)-1)/2,(len(B1)-1)/2+1)
            bunch_aux.update({j: {'atEncounters' : aux,'atPositions':myPosition[np.in1d(B1,B1_fillingScheme)]}})
            B2_BB_pattern.update({'at'+format(i,'04d'):bunch_aux})

    beam_BB_pattern=dotdict({'atB1':B1_BB_pattern,'atB2':B2_BB_pattern})
    return beam_BB_pattern

def B1CollisionScheduleDF (B1_bunches, B2_bunches, numberOfLRToConsider):
    
    '''
    For given two series of booleans which represent bunches and a array that represent long range collisions, 
    this function returns dataframe related to their collisions from perspective of beam 1.
    
    ===EXAMPLE===
    fillingSchemeDF=importData.LHCFillsAggregation(['LHC.BCTFR.A6R4.B%:BUNCH_FILL_PATTERN'],6666, ['FLATTOP'],flag='next')
    B1_bunches = fillingSchemeDF['LHC.BCTFR.A6R4.B1:BUNCH_FILL_PATTERN'].iloc[0]
    B2_bunches = fillingSchemeDF['LHC.BCTFR.A6R4.B2:BUNCH_FILL_PATTERN'].iloc[0]
    B1CollisionScheduleDF(B1_bunches, B2_bunches, 25)
    
    '''
    bunch_value = 1.0
    # Transforming bunches in to boolean array
    B1_bunches = np.array(B1_bunches) == bunch_value
    B2_bunches = np.array(B2_bunches) == bunch_value
    
    # For debugging
    # pdb.set_trace()
    
    # Get indexes of Bunches
    B1_bunches_index = np.where(B1_bunches)[0]
    B2_bunches_index = np.where(B2_bunches)[0]
    
    if isinstance(numberOfLRToConsider, int):
        numberOfLRToConsider = [numberOfLRToConsider, numberOfLRToConsider, numberOfLRToConsider]
    
    B1df = pd.DataFrame() 
       
    for n in B1_bunches_index:
    
    
        # First check for collisions in ALICE

        # Formula for head on collision in ALICE is 
        # (n + 891) mod 3564 = m
        # where n is number of bunch in B1, and m is number of bunch in B2
        
        # Formula for head on collision in ATLAS/CMS is 
        # n = m
        # where n is number of bunch in B1, and m is number of bunch in B2
        
        # Formula for head on collision in LHCb is 
        # (n + 2670) mod 3564 = m
        # where n is number of bunch in B1, and m is number of bunch in B2

        head_on_names = ["HO partner in ALICE", "HO partner in ATLAS/CMS", "HO partner in LHCB"]
        secondary_names = ["# of LR in ALICE", "# of LR in ATLAS/CMS", "# of LR in LHCB"]
        encounters_names = ["BB partners in ALICE", "BB partners in ATLAS/CMS", "BB partners in LHCB"]
        positions_names = ["Positions in ALICE", "Positions in ATLAS/CMS", "Positions in LHCB"]
        
        colide_factor_list = [891, 0, 2670]
        number_of_bunches = 3564
        
        # i == 0 for ALICE
        # i == 1 for ATLAS and CMS
        # i == 2 for LHCB
        
        dictonary = {}
        
        for i in range(0,3):
        
            collide_factor = colide_factor_list[i]
            m = (n + collide_factor) % number_of_bunches
            
            #pdb.set_trace()
            # if this Bunch is true, than there is head on collision
            if B2_bunches[m]:
                head_on = m
            else:
                head_on = np.nan
            
            ## Check if beam 2 has bunches in range  m - numberOfLRToConsider to m + numberOfLRToConsider 
            ## Also have to check if bunches wrap around from 3563 to 0 or vice versa
            
            
            bunches_ineraction_temp = np.array([])
            encounters = np.array([])
            positions = np.array([])
            
            first_to_consider = m - numberOfLRToConsider[i]
            last_to_consider = m + numberOfLRToConsider[i] + 1
            
            numb_of_long_range = 0
                
            if first_to_consider < 0:
                bunches_ineraction_partial = np.where(B2_bunches[(number_of_bunches + first_to_consider):(number_of_bunches)])[0]
                                
                # This represents the absolute position of the bunches 
                encounters = np.append(encounters, number_of_bunches + first_to_consider + bunches_ineraction_partial )
                
                #This represents the relative position to the head-on bunch
                positions = np.append(positions, first_to_consider + bunches_ineraction_partial)
                
                #Set this varibale so later the normal syntax wihtout the wrap around checking can be used
                first_to_consider = 0
                
            if last_to_consider > number_of_bunches:
                bunches_ineraction_partial = np.where(B2_bunches[0:last_to_consider - number_of_bunches])[0]
                                
                # This represents the absolute position of the bunches 
                encounters = np.append(encounters, bunches_ineraction_partial )
                
                #This represents the relative position to the head-on bunch
                positions = np.append(positions, number_of_bunches - m + bunches_ineraction_partial)
                
                
                last_to_consider = number_of_bunches
            
            bunches_ineraction_partial = np.append(bunches_ineraction_temp, np.where(B2_bunches[first_to_consider:last_to_consider])[0])
                        
            # This represents the absolute position of the bunches 
            encounters = np.append(encounters, first_to_consider + bunches_ineraction_partial)
            
            #This represents the relative position to the head-on bunch
            positions = np.append(positions, bunches_ineraction_partial - (m - first_to_consider))
            
            
            # Substract head on collision from number of secondary collisions
            numb_of_long_range = len(positions) - int(B2_bunches[m])

            
            dictonary.update({head_on_names[i] : { n : head_on }, secondary_names[i] : { n : numb_of_long_range }, \
                             encounters_names[i] : { n : encounters },  positions_names[i] : { n : positions } })

        B1df = pd.concat([B1df,pd.DataFrame(dictonary)])


    return B1df

def B2CollisionScheduleDF (B1_bunches, B2_bunches, numberOfLRToConsider):
    
    '''
    For given two series of booleans which represent bunches and a array that represent long range collisions, 
    this function returns dataframe related to their collisions from perspective of beam 1.
    
    ===EXAMPLE===
    fillingSchemeDF=importData.LHCFillsAggregation(['LHC.BCTFR.A6R4.B%:BUNCH_FILL_PATTERN'],6666, ['FLATTOP'],flag='next')
    B1_bunches = fillingSchemeDF['LHC.BCTFR.A6R4.B1:BUNCH_FILL_PATTERN'].iloc[0]
    B2_bunches = fillingSchemeDF['LHC.BCTFR.A6R4.B2:BUNCH_FILL_PATTERN'].iloc[0]
    B2CollisionScheduleDF(B1_bunches, B2_bunches, 25)
    
    '''
    bunch_value = 1.0
    # Transforming bunches in to boolean array
    B1_bunches = np.array(B1_bunches) == bunch_value
    B2_bunches = np.array(B2_bunches) == bunch_value
    
    # For debugging
    # pdb.set_trace()
    
    # Get indexes of Bunches
    B1_bunches_index = np.where(B1_bunches)[0]
    B2_bunches_index = np.where(B2_bunches)[0]
    
    if isinstance(numberOfLRToConsider, int):
        numberOfLRToConsider = [numberOfLRToConsider, numberOfLRToConsider, numberOfLRToConsider]
    
    B2df = pd.DataFrame() 
       
    for n in B2_bunches_index:
    
    
        # First check for collisions in ALICE

        # Formula for head on collision in ALICE is 
        # (n + 891) mod 3564 = m
        # where n is number of bunch in B1, and m is number of bunch in B2
        
        # Formula for head on collision in ATLAS/CMS is 
        # n = m
        # where n is number of bunch in B1, and m is number of bunch in B2
        
        # Formula for head on collision in LHCb is 
        # (n + 2670) mod 3564 = m
        # where n is number of bunch in B1, and m is number of bunch in B2

        head_on_names = ["HO partner in ALICE", "HO partner in ATLAS/CMS", "HO partner in LHCB"]
        secondary_names = ["# of LR in ALICE", "# of LR in ATLAS/CMS", "# of LR in LHCB"]
        encounters_names = ["BB partners in ALICE", "BB partners in ATLAS/CMS", "BB partners in LHCB"]
        positions_names = ["Positions in ALICE", "Positions in ATLAS/CMS", "Positions in LHCB"]
        
        colide_factor_list = [891, 0, 2670]
        number_of_bunches = 3564
        
        # i == 0 for ALICE
        # i == 1 for ATLAS and CMS
        # i == 2 for LHCB
        
        dictonary = {}
        
        for i in range(0,3):
        
            collide_factor = colide_factor_list[i]
            m = (n - collide_factor) % number_of_bunches
            
            #pdb.set_trace()
            # if this Bunch is true, than there is head on collision
            if B1_bunches[m]:
                head_on = m
            else:
                head_on = np.nan
            
            ## Check if beam 2 has bunches in range  m - numberOfLRToConsider to m + numberOfLRToConsider 
            ## Also have to check if bunches wrap around from 3563 to 0 or vice versa
            
            
            bunches_ineraction_temp = np.array([])
            encounters = np.array([])
            positions = np.array([])
            
            first_to_consider = m - numberOfLRToConsider[i]
            last_to_consider = m + numberOfLRToConsider[i] + 1
            
            numb_of_long_range = 0
                
            if first_to_consider < 0:
                bunches_ineraction_partial = np.where(B1_bunches[(number_of_bunches + first_to_consider):(number_of_bunches)])[0]
                                
                # This represents the absolute position of the bunches 
                encounters = np.append(encounters, number_of_bunches + first_to_consider + bunches_ineraction_partial )
                
                #This represents the relative position to the head-on bunch
                positions = np.append(positions, first_to_consider + bunches_ineraction_partial)
                
                #Set this varibale so later the normal syntax wihtout the wrap around checking can be used
                first_to_consider = 0
                
            if last_to_consider > number_of_bunches:
                bunches_ineraction_partial = np.where(B1_bunches[0:last_to_consider - number_of_bunches])[0]
                                
                # This represents the absolute position of the bunches 
                encounters = np.append(encounters, bunches_ineraction_partial )
                
                #This represents the relative position to the head-on bunch
                positions = np.append(positions, number_of_bunches - m + bunches_ineraction_partial)
                
                
                last_to_consider = number_of_bunches
            
            bunches_ineraction_partial = np.append(bunches_ineraction_temp, np.where(B1_bunches[first_to_consider:last_to_consider])[0])
                        
            # This represents the absolute position of the bunches 
            encounters = np.append(encounters, first_to_consider + bunches_ineraction_partial)
            
            #This represents the relative position to the head-on bunch
            positions = np.append(positions, bunches_ineraction_partial - (m - first_to_consider))            
            
            # Substract head on collision from number of secondary collisions
            numb_of_long_range = len(positions) - int(B1_bunches[m])

            
            dictonary.update({head_on_names[i] : { n : head_on }, secondary_names[i] : { n : numb_of_long_range }, \
                             encounters_names[i] : { n : encounters },  positions_names[i] : { n : positions } })

        B2df = pd.concat([B2df,pd.DataFrame(dictonary)])


    return B2df

def _BeamFilling2str(booleanList, flagChars = ['e', 'b']):
    '''
    
    For a given boolean list, this function returns a string which represents the counting operation of boolean values. For example
    booleanList = [False, False, True, False, False, False] would return "2e1b3e"
    Character flags that represent true or false can be edited with flagChars variable
    
    b - bunch ~ True
    e - empty ~ False
    
    ===EXAMPLE===
    genericFillingPatern([True, True, True, False, False, False, True])
        
    '''

    
    #pdb.set_trace()
    
    counter = [0, 0]
    result = ''
    
    lengthList = len(booleanList)
    if lengthList == 0:
        raise Error()
    currentStatus = booleanList[0]
    
    for i in range(0, lengthList):
        
        if currentStatus == booleanList[i]:
            counter[int(currentStatus)] = counter[int(currentStatus)] + 1
            if i == lengthList - 1:
                result = result + str(counter[int(currentStatus)]) + flagChars[int(currentStatus)]    
        else:
            result = result + str(counter[int(currentStatus)]) + flagChars[int(currentStatus)]
            counter[int(currentStatus)] = 0
            
            currentStatus = not currentStatus
            counter[int(currentStatus)] = 1
            
            if i == lengthList - 1:
                result = result + str(counter[int(currentStatus)]) + flagChars[int(currentStatus)]
            
    return result

def BeamFilling2str (bunchFilling, flagChars = ['e', 'b']):
    
    '''
    
    Simple function that addapts the _BeamFilling2str for use in the cals2pd package
    
    Text from _BeamFilling2str:    
    For a given boolean list, this function returns a string which represents the counting operation of boolean values. For example
    booleanList = [False, False, True, False, False, False] would return "2e1b3e"
    Character flags that represent true or false can be edited with flagChars variable
    
    b - bunch ~ True
    e - empty ~ False
    
    ===EXAMPLE===
    fillingSchemeDF=importData.LHCFillsAggregation(['LHC.BCTFR.A6R4.B%:BUNCH_FILL_PATTERN'],6666, ['FLATTOP'],flag='next')
    BeamFilling2str(fillingSchemeDF['LHC.BCTFR.A6R4.B1:BUNCH_FILL_PATTERN'].iloc[0])
        
    '''
    
    bunchesBoolean = np.array(bunchFilling) == 1.0
    return _BeamFilling2str(bunchesBoolean, flagChars)



def CollisionsInTheGivenSlot_vec(B1_filling_scheme, B2_filling_scheme, slot_number, numberOfLRToConsider):
    '''
    This function returns three vectors with the information, respectively, about 
    the partner bunch of filling_scheme in the collision with filling_scheme_to_be_rolled, 
    about the position of the collision saved in the previous vector, and finally about the
    total number of collision for that particle of the filling_scheme_to_be_rolled
    Args:
        filling_scheme_to_be_rolled: a boolean np.array, that represent the fillig scheme
        to be rolled in order to compute the collision
        filling_scheme: a boolean np.array, that represent the other filling scheme
        IPN_pos: int that give the position of the detector, around which we want to compute the LR
        n_LR: is a int that represent how many LR the user want to consider around the detector
    Return:
        np.array(v): information about the partner of the collision (filling_scheme)
        np.array(v_pos): information about the position around the detector where happen the collision
        np.array(tot_LR): how many collision for that particle in filling_scheme_to_be_rolled
        or 
        pd.DataFrame that give the information about the partner beam collimator "BB in LR"
        and the position of the LR 
        np.array that give the information about the posi( to verify)
    '''
    if numberOfLRToConsider== 0:
        n_bunches = len(B1_filling_scheme)
        filling_scheme_rolled = np.concatenate([B1_filling_scheme[-(slot_number%n_bunches):],\
            B1_filling_scheme[:-(slot_number%n_bunches)]])
        index_events = (filling_scheme_rolled*B2_filling_scheme)
        return  pd.DataFrame({'BB_collision' : np.where(B2_filling_scheme*index_events)[0], 'pos_collision' : slot_number},\
            index = (np.where(filling_scheme_rolled*index_events)[0]+(-slot_number))%n_bunches), \
            np.concatenate([index_events[-((-slot_number)%n_bunches):],index_events[:-((-slot_number)%n_bunches)]])
    else:
        n_bunches = len(B1_filling_scheme)
        v1 = np.empty([n_bunches,2*numberOfLRToConsider+1])
        v1[:] = np.nan
        v1_pos = np.empty([n_bunches,2*numberOfLRToConsider+1])
        v1_pos[:] = np.nan
        tot_LR = np.zeros(n_bunches)
        count = 0
        for i in (np.arange(numberOfLRToConsider*2+1) - numberOfLRToConsider + slot_number):
            if i !=slot_number:
                s = np.empty([n_bunches,1])
                s[:] = np.nan
                s_pos = np.empty([n_bunches,1])
                s_pos[:] = np.nan
                filling_scheme_rolled = np.concatenate([B1_filling_scheme[-(i%n_bunches):],\
                B1_filling_scheme[:-(i%n_bunches)]])
                index_events = (filling_scheme_rolled*B2_filling_scheme)   
                v1[(np.where(index_events)[0]+(-i))%n_bunches,count] = np.concatenate([s[(np.where(index_events)[0]+(-numberOfLRToConsider))%n_bunches],\
                    np.array(list(np.where(index_events))).T],axis = 1)[:,1]
                
                
                v1_pos[(np.where(index_events)[0]+(-i))%n_bunches,count] = np.concatenate([s_pos[(np.where(index_events)[0]+(-numberOfLRToConsider))%n_bunches],\
                    i*np.ones([len(np.where(index_events)[0]),1])],axis = 1)[:,1]
                count +=1
                tot_LR +=np.concatenate([index_events[-((-i)%n_bunches):],index_events[:-((-i)%n_bunches)]])
        t = ~np.isnan(v1)
        v = [v1[ii][t[ii]] for ii in np.arange(n_bunches)]

        t_pos = ~np.isnan(v1)
        v_pos = [v1_pos[ii][t_pos[ii]] for ii in np.arange(n_bunches)]
        return v1 ,v1_pos ,tot_LR

def events_in_slots(filling_scheme_to_be_rolled, filling_scheme,LR):
    '''
    Create a dataframe that represent the particle of the filling_scheme that 
    collide with the particles of the filling_scheme_to_be_rolled in the position LR,
    moreover it returns a vector that give the information of how many collision 
    Args:
        filling_scheme_to_be_rolled: a boolean np.array, that represent the fillig scheme
        to be rolled in order to compute the collision
        filling_scheme: a boolean np.array, that represent the other filling scheme
        LR: a int that give the position where to compute the collision of the two beams
    Return:
        pd.DataFrame that give the information about the partner beam collimator "BB in LR"
        and the position of the LR 
        np.array that give the information about the posi( to verify)
    '''
    n_bunches = len(filling_scheme_to_be_rolled)
    filling_scheme_rolled = np.concatenate([filling_scheme_to_be_rolled[-(LR%n_bunches):],\
        filling_scheme_to_be_rolled[:-(LR%n_bunches)]])
    index_events = (filling_scheme_rolled*filling_scheme)
    return  pd.DataFrame({'BB_LR' : np.where(filling_scheme*index_events)[0], 'pos_LR' : LR},\
    index = (np.where(filling_scheme_rolled*index_events)[0]+(-LR))%n_bunches), \
        np.concatenate([index_events[-((-LR)%n_bunches):],index_events[:-((-LR)%n_bunches)]])
 

def events_in_slots_vec(filling_scheme_to_be_rolled, filling_scheme, IPN_pos, n_LR):
    '''
    This function returns three vectors with the information, respectively, about 
    the partner bunch of filling_scheme in the collision with filling_scheme_to_be_rolled, 
    about the position of the collision saved in the previous vector, and finally about the
    total number of collision for that particle of the filling_scheme_to_be_rolled
    Args:
        filling_scheme_to_be_rolled: a boolean np.array, that represent the fillig scheme
        to be rolled in order to compute the collision
        filling_scheme: a boolean np.array, that represent the other filling scheme
        IPN_pos: int that give the position of the detector, around which we want to compute the LR
        n_LR: is a int that represent how many LR the user want to consider around the detector
    Return:
        np.array(v): information about the partner of the collision (filling_scheme)
        np.array(v_pos): information about the position around the detector where happen the collision
        np.array(tot_LR): how many collision for that particle in filling_scheme_to_be_rolled
    '''
    n_bunches = len(filling_scheme_to_be_rolled)
    v1 = np.empty([n_bunches,2*n_LR+1])
    v1[:] = np.nan
    v1_pos = np.empty([n_bunches,2*n_LR+1])
    v1_pos[:] = np.nan
    tot_LR = np.zeros(n_bunches)
    count = 0
    for i in (np.arange(n_LR*2+1) - n_LR + IPN_pos):
        if i !=IPN_pos:
            s = np.empty([n_bunches,1])
            s[:] = np.nan
            s_pos = np.empty([n_bunches,1])
            s_pos[:] = np.nan
            filling_scheme_rolled = np.concatenate([filling_scheme_to_be_rolled[-(i%n_bunches):],\
            filling_scheme_to_be_rolled[:-(i%n_bunches)]])
            index_events = (filling_scheme_rolled*filling_scheme)   
            v1[(np.where(index_events)[0]+(-i))%n_bunches,count] = np.concatenate([s[(np.where(index_events)[0]+(-n_LR))%n_bunches],\
                np.array(list(np.where(index_events))).T],axis = 1)[:,1]
            
            
            v1_pos[(np.where(index_events)[0]+(-i))%n_bunches,count] = np.concatenate([s_pos[(np.where(index_events)[0]+(-n_LR))%n_bunches],\
                i*np.ones([len(np.where(index_events)[0]),1])],axis = 1)[:,1]
            count +=1
            tot_LR +=np.concatenate([index_events[-((-i)%n_bunches):],index_events[:-((-i)%n_bunches)]])
    t = ~np.isnan(v1)
    v = [v1[ii][t[ii]] for ii in np.arange(n_bunches)]

    t_pos = ~np.isnan(v1)
    v_pos = [v1_pos[ii][t_pos[ii]] for ii in np.arange(n_bunches)]
    return v ,v_pos ,tot_LR



def bbschedule(bool_slotsB1,bool_slotsB2, numberOfLRToConsider, Dict1 = {'ATLAS/CMS':0,'LHCB':2670, 'ALICE': 891}, n_slots = 3564, Dict_nLRs = {'Nan' : np.NaN}):
    ''' 
    This function return two pd.DataFrame, associated to the two beams, with the information 
    about the purtner bunch for HO and LR, about the position of collision with respect 
    to the detector and also how many LR there are around that detector. The two filling schemes
    are loaded directly from the two boolean vectors as input.
    Args:
        bool_slotsB1: a boolean np.array that give the information about the position of the
        filling scheme of B1 
        bool_slotsB2: a boolean np.array that give the information about the position of the 
        filling scheme of B2
        numberofLRToConsider : a number of LR to consider in every detector of the accellerator
        Dic1 : a dictionary that give the position of the detectors, seen from the point of view of B1, 
        and also how many LR are considered in each detector, set with the ones of LHC
        n_slots : number of length of the accellerator, in this case is set as the number of harmonic RF of LHC divided by 10
        Dict_nLRs : a dictionary that give the number of LR to consider in every detector of the dictionary Dict1, 
        they must have the same length, is set as np.nan
    Returns:
        pd.DataFrame (df_B1): that give all the information described above about LR and HO, seen from
        a poinf of view of B1
        pd.DataFrame (df_B2): that give all the information described above about LR and HO, seen from
        a point of view of B2

    '''

    assert np.isnan([numberOfLRToConsider,Dict_nLRs[random.choice(list(Dict_nLRs.keys()))]]).any(), "one between numberOfLRToConsider and Dict_nLRs should be NaN"
    assert ~np.isnan([numberOfLRToConsider,Dict_nLRs[random.choice(list(Dict_nLRs.keys()))]]).all(), "one between numberOfLRToConsider and Dict_nLRs should be NaN"
    keys = list(Dict1.keys())#
    if np.isnan(Dict_nLRs[random.choice(list(Dict_nLRs.keys()))]):
        Dic1 = {}
        Dic2 = {}
        for i in range(len(keys)):
            Dic1[f'{keys[i]}'] = [Dict1[f'{keys[i]}'],numberOfLRToConsider]
            Dic2[f'{keys[i]}'] = [np.mod(n_slots-Dict1[f'{keys[i]}'],n_slots),numberOfLRToConsider]
    else:
        keys_LR = list(Dict_nLRs.keys())
        assert len(keys) == len(keys_LR), "Dict1 and Dict_nLRs must have the same lengths of keys"
        Dic1 = {}
        Dic2 = {}
        for i in range(len(keys)):
            Dic1[f'{keys[i]}'] = [Dict1[f'{keys[i]}'],Dict_nLRs[f'{keys_LR[i]}']]
            Dic2[f'{keys[i]}'] = [np.mod(n_slots-Dict1[f'{keys[i]}'],n_slots),Dict_nLRs[f'{keys_LR[i]}']]
    
    
    
    ones_B1 = np.where(bool_slotsB1)
    ones_B2 = np.where(bool_slotsB2)


    df_B1 = pd.DataFrame(index = ones_B1[0])

    df_B2 = pd.DataFrame(index = ones_B2[0])
    
    for j in np.arange(len(Dic1.keys())):
        
        IPN_B1 = list(Dic1.keys())[j]
        [v,v_pos,sum_v]\
         = events_in_slots_vec(bool_slotsB1,bool_slotsB2,Dic1[IPN_B1][0],Dic1[IPN_B1][1])
        df_B1[f'# of LR in {IPN_B1}'] = sum_v[ones_B1[0]]
        df_B1[f'HO partner in {IPN_B1}'] = events_in_slots(bool_slotsB1,bool_slotsB2,Dic1[IPN_B1][0])[0].iloc[:,0]
        
        df_B1[f'BB partners in {IPN_B1}'] = [list(v[ii]) for ii in ones_B1[0]]
        df_B1[f'Positions in {IPN_B1}'] = [list(v_pos[ii]-Dic1[IPN_B1][0]) for ii in ones_B1[0]]
    
    for j in np.arange(len(Dic1.keys())):
        
        IPN_B2 = list(Dic2.keys())[j]
        [v,v_pos,sum_v]\
         = events_in_slots_vec(bool_slotsB2,bool_slotsB1,Dic2[IPN_B2][0],Dic2[IPN_B2][1])
        df_B2[f'# of LR in {IPN_B2}'] = sum_v[ones_B2[0]]
        df_B2[f'HO partner in {IPN_B2}'] = events_in_slots(bool_slotsB2,bool_slotsB1,Dic2[IPN_B2][0])[0].iloc[:,0]
        df_B2[f'BB partners in {IPN_B2}'] = [list(v[ii]) for ii in ones_B2[0]]
        df_B2[f'Positions in {IPN_B2}'] = [list(v_pos[ii]-Dic1[IPN_B2][0]) for ii in ones_B2[0]]
    return df_B1,df_B2