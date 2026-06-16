import numpy as np
import os

#tot_masses = np.array([1.00])
#shell_masses = np.array([0.050])
tot_masses = np.array([0.90, 1.00, 1.10, 1.20])
shell_masses = np.array([0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100])

#os.chdir('/"path2Cattails"/Cattails/Typha/WD_maker/')
os.chdir('/home/dhasenou/Research/Cattails/Typha/WD_maker/')

# clear output file and make a new one
if os.path.isfile("./output.txt") == True:
    os.system("rm output.txt")
    os.system("touch output.txt")
# make directory for WD model files
if os.path.isdir("./initial_WD_models") == False:
    os.system("mkdir ./initial_WD_models")

for mtot in tot_masses:
    for mshell in shell_masses:
        mcore = mtot - mshell
        if mcore > 0.2: # WD mass should be at least .2 Msun
            # execute the WD maker script with the mtot and mshell value
            print("executing WD_maker.exe for tot mass = %.3f and shell mass = %.3f... " %(mtot, mshell))
            os.system('./WD_maker.exe %.3f %.3f >> output.txt' % (mtot, mshell))
            print("finished.")
            # rename the initial.dat file 
            print("renaming initial.dat ... ")
            totmass_str = str(mtot).replace(".","_")
            shellmass_str = str(mshell).replace(".","_")
            filename = 'initial-' + totmass_str + '-' + shellmass_str + '.dat'
            os.system("mv initial.dat ./initial_WD_models/%s" % filename)
