import numpy as np
from numpy import linalg as LA
import time


### EXAMPLE ON HOW TO USE THIS:
### IN MAIN SCRIPT FILE,

#import numpy as np
#from numpy import linalg as LA
#import matplotlib.pyplot as plt
#import time
#from random import seed
#from random import random
#
#from Cq_generate_concatenated_utilities import Cq_generate_concatenated
#
## seed random number generator
#seed(1)
#
#################################################
### GENERATE FICTICIOUS BITSTREAM FOR TESTING ###
#################################################
#
### RAND
#bitlength = 16000
#bits1 = np.array(np.zeros(bitlength)) # All zeros
#for i in range(0,bitlength):
#    x = round(random())
#    bits1[i] = x
#
### 3-2 GMP
#file = open('3-2GMP.txt','r')
#data = file.read().split(' ')
#bits2 = np.array(np.zeros(len(data))) # All zeros
#for i in range(0,len(data)):
#    bits2[i] = data[i]
#
### ASSIGN BITS TO BE WHICHEVER PROCESS:
#bits = bits2


### CODE WAS TRANSLATED FROM MATLAB TO PYTHON3 ON 2 MAY 2020.
### DISCLAIMER: THE CODE MAY NOT BE OPTIMISED SINCE IT WAS TRANSLATED
###             DIRECTLY FROM MATLAB. MATTHEW HAS BEEN A MATLAB USER
###             FOR MANY YEARS AND HAS ONLY JUST BEGUN THE SWITCH
###             TO PYTHON.


def Cq_generate_concatenated(L,bits):

    ###############

    t = time.time()
    
    #################
    ### BITSTREAM ###
    #################
    
    lengthOfBits = (bits.shape[0])
    
    ######################################################
    ### USER INPUT LENGTH OF WORD SEQUENCES TO ANALYSE ###
    ######################################################
    
    L = L
    # print('Observable sequence length L =',L)
    
    #####################################################
    ### TO CONSTRUCT MATRIX WITH INCREASING TIMESTEPS ###
    #####################################################
    
    bitsMatrix = np.zeros((lengthOfBits-2*L+L,L+1))

    [rows, cols] = bitsMatrix.shape
    # print(rows,cols)
    # print(bitsMatrix.shape[0],bitsMatrix.shape[1])
    
    for i in range(0,lengthOfBits-2*L+L):
        bitsMatrix[i,:] = bits[i:(L+1)+i]
    
    # print(bitsMatrix)

    #################################
    ### CONDITIONAL PROBABILITIES ###
    #################################

    condCounts = np.zeros((2**L,2))

    for i in range(0,bitsMatrix.shape[0]):

        ### ROW NUMBER
        row_toconvert = list(map(int, bitsMatrix[i,0:L]))
        binary_string = ""
        for x in row_toconvert:
            binary_string += str(x)

        which_row = int(binary_string,2)
        # print('ROW: bin =', binary_string,', dec = ',which_row)

        ### COL NUMBER
        col_toconvert = bitsMatrix[i,L]
        which_col = int(col_toconvert)
        # print('COL: bin = ',col_toconvert,', dec = ',which_col)

        ### ADD TO CONDCOUNTS
        condCounts[which_row,which_col] += 1

    # print('condCount: \n',condCounts)

    condProb = np.zeros((2**L,2))
    for i in range(0,2**L):
        for j in range(0,2):
            condProb[i,j] = condCounts[i,j]/sum(condCounts[i,:])
            if np.isnan(condProb[i,j]) == True:
                condProb[i,j] = 0

    # print('condProbs: \n',np.round(condProb,4))

    ##########################################
    ### TO GENERATE MATRIX OF SIZE (2^L,L) ###
    ### EACH ROW CORRESPONDS TO THE BINARY ###
    ### REPRESENTATION OF EACH ROW NUMBER  ###
    ### I.E. GENERATE PASTS/FUTURES        ###
    ##########################################

    base = 2
    maxdecimal = base**L
    nums = np.arange(0,maxdecimal)

    past = np.zeros((2**L,L))
    for i in range(0,maxdecimal):
        bin_string = np.binary_repr(i,width=L)
        bin_int = list(map(int, bin_string))
        past[i,:] = (bin_int) 

    future = past
    # print('Pasts: \n',past,'\nFutures: \n',future)


    #######################################################
    ### EVALUATING 2^L BY 2^L CONDITIONAL PROBABILITIES ###
    ### I.E. CONCATENATING THE SINGLE FUTURE COND PROBS ###
    #######################################################

    condfut_no_need_to_normalise = np.zeros((2**L,2**L))
    # print(condfut_no_need_to_normalise)

    for i in range(0,2**L):
        for j in range(0,2**L):

            # STEP (1) -- Past-Future combination
            PF = np.hstack((past[i,:],future[j,:]))

            # STEP(2) -- Concatenated-combinations
            mat = np.zeros((L,L+1))
            for k in range(0,L):
                mat[k,:] = PF[k:k+L+1]

            # STEP (3) -- Getting indiv. single future condprobs
            prob_from_mat = np.zeros((1,L))
    #         print(prob_from_mat)
            for k in range(0,L):

                ### ROW NUMBER
                row_string = list(map(int, mat[k,0:L]))
                binary_string = ""
                for x in row_string:
                    binary_string += str(x)
                row_in_dec = int(binary_string,2)

                ### COL NUMBER
                col_string = mat[k,L]
                col_in_dec = int(col_string)

    #             print(mat[k,0:L], row_in_dec, '......[',mat[k,L],']', col_in_dec)
                prob_from_mat[0,k] = condProb[row_in_dec,col_in_dec]

    #         print(prob_from_mat)
            condfut_no_need_to_normalise[i,j] = np.prod(prob_from_mat)

    # print('condProbs_Lfutures: \n',np.round(condfut_no_need_to_normalise,4)) 
    # Sum of each row = 1

    ######################################################
    ### SQRT THE COND PROBABILITIES FOR QUANTUM STATES ###
    ######################################################

    sqrt_condFuture_probs = np.sqrt(condfut_no_need_to_normalise)
    # print('Sqrt_condfutures: \n',np.round(sqrt_condFuture_probs,4))
    # Each row contains probability amplitudes of future combinations

    #########################################
    ### FINDING PROBABILITY OF EACH STATE ###
    #########################################

    bitsMatrix2 = np.zeros((lengthOfBits-L+1,L))
    for i in range(0,lengthOfBits-L+1):
        bitsMatrix2[i,:] = bits[i:L+i]

    stateProbCount = np.zeros((2**L,1))
    for i in range(0,lengthOfBits-L+1):

        row_string = list(map(int, bitsMatrix2[i,:]))
        binary_string = ""
        for x in row_string:
            binary_string += str(x)
        row_in_dec = int(binary_string,2)

        stateProbCount[row_in_dec,0] += 1

    stateProbVec = np.zeros((2**L,1))
    stateProbVec = stateProbCount/np.sum(stateProbCount)

    # print('Prob of each state, stateProbVec: \n',stateProbVec)


    # ### FINDING \rho ###
    ####################

    rho = np.zeros((2**L,2**L))
    # print('rho: \n',np.round(rho,4))
    for i in range(0,2**L):
        state = np.matrix(sqrt_condFuture_probs[i,:])
        state_prob = stateProbVec[i]
        rho_temp = state_prob[0] * (np.transpose(state) * state)
        rho += rho_temp

    # print('rho: \n',np.round(rho,4))



    ###################################
    ### FINDING EIGENVALUES OF \rho ###
    ###################################

    #The normalized (unit “length”) eigenvectors, such that the column v[:,i]
    #is the eigenvector corresponding to the eigenvalue w[i].
    eigvals, eigvecs = LA.eig(rho)
    # print(np.round(eigvals,15))
    for i in range(0,2**L):
        if eigvals[i] < 0:
            eigvals[0] = 0

    log2eigvals = np.real(np.log2(eigvals))

    # print(np.round(log2eigvals,2))
    for i in range(0,2**L):
        if np.isinf(-log2eigvals[i]) == True:
            log2eigvals[i] = 0

        if np.isnan(log2eigvals[i]) == True:
            log2eigvals[i] = 0

    # print(np.round(log2eigvals,2))

    ######################
    ### CALCULATING CQ ###
    ######################

    eigvals_log2eigvals = eigvals*log2eigvals

    # print(np.round(eigvals_log2eigvals,2))

    Cq_concatenated = (-1)*sum(eigvals_log2eigvals)
    # print('Cq = ',np.round(Cq_concatenated,4))




    ##############################################################
    elapsed = time.time() - t
    # print('------------------------------------------------------')
    print('Code completed in ',elapsed,' seconds')
    
    return Cq_concatenated, sqrt_condFuture_probs, stateProbVec, rho 