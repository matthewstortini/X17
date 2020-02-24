import numpy as np
import h5py
import pandas as pd
import sys
import matplotlib.pyplot as plt
plt.style.use('style.mplstyle')

def main():

    spectrum()

def spectrum():
    """
    This function simply plots an energy histogram. Run postprocesshdf5.py (function process()) on g4simple output file to get desired files 
    for the dataframes defined below.
    """ 

    if(len(sys.argv) != 2):
        print('Usage: spectrum.py [input .hdf5 file (with extension)]')
        sys.exit()

    # read in pandas dataframe, and pull out data of interest
    df =  pd.read_hdf("{}".format(sys.argv[1]), key="procdf")

    # set up data to plot, and perform any necessary unit conversion
    m = list(df['theta'])

    print(len(m))

    # plot simulation data
    plt.hist(m, np.arange(0,190,0.1), histtype='step', color = 'black')
    
    # set plot aesthetics
    plt.xlim(0,190)
    plt.ylim(0,plt.ylim()[1])
    plt.xlabel('theta', ha='right', x=1.0)
    plt.ylabel('counts', ha='right', y=1.0)
    plt.tight_layout()
    #plt.semilogy()
    #plt.semilogx()
    plt.show()


if __name__ == '__main__':
	main()

