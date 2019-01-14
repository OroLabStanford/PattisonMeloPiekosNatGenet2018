# PattisonMeloPiekosNatGen2018
These are the latest versions of the customized code utilized in the Pattison et al. 2018 Nature Genetics paper. Author of the code is Samantha Piekos.

AnchorLoops.py is a custom python3 script that identifies chromatin contacts that contains an element of interest in either of the two bins. Prints out the chromatin contact and feature of interest pairs to a new file. This output file is in the correct format for the HiChIP input file for the Deg1LoopChecker.py script.

Deg1LoopChecker.py is a custom python3 script that identifies two distal genomic elements that are connected to each other via chromatin looping ensuring that each element is in it's own individual contact bin. It writes information on the coordinates of the two elements and the associated chromatin loop to an output file.
