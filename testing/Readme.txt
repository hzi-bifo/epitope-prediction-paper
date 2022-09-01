This folder contains all the datasets that can be tested directly with the trained EpitopeVec models.
Please use the test.py (in the root directory) file to test.
usage: python test.py modelname inputfilename


abcpred16.txt: ABCPred peptides with fixed length of 16

bcpred_20.txt: BCPreds peptides with fixed length of 20

blind_full.txt: Blind dataset published with ABCPred method. The length of the peptides is not fixed.

chen_full.txt: Peptides used in training of AAP method by Chen et al. Fixed length of 20

ibce-ind.txt: peptides used for the training of iBCE-EL method. variable length

ibce-ind.txt: peptides used for independent testing of the iBCE-EL method. variable length

lbtope_fixed_nr-combined.txt: Non-renundant peptides with fixed length of 20 published with the LBTope method.
