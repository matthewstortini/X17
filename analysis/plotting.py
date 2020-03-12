import numpy as np
import h5py
import pandas as pd
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import json, os
plt.style.use('style.mplstyle')

def main():

    angles()
    #depths()
    #resonance_positions()
    #spectrum()

def angles():

    if(len(sys.argv) != 2):
        print('Usage: plotting.py [input .hdf5 file (with extension)]')
        sys.exit()

    with open("data.json") as f:
        data = json.load(f)
    data_dir = os.path.expandvars(data["data_dir"])

    # read in pandas dataframe, and pull out data of interest
    df =  pd.read_hdf("{}/{}".format(data_dir,sys.argv[1]), key="procdf")    

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


def depths():

    if(len(sys.argv) != 2):
        print('Usage: plotting.py [input .hdf5 file (with extension)]')
        sys.exit()

    with open("data.json") as f:
        data = json.load(f)
    data_dir = os.path.expandvars(data["data_dir"])

    # read in pandas dataframe, and pull out data of interest
    df =  pd.read_hdf("{}/{}".format(data_dir,sys.argv[1]), key="procdf")
    df['z'] = df['z']*1000+7.5
    df['KE'] = df['KE']*1000

    plt.hist2d(df['z'], df['KE'], np.arange(-1,1001,1), norm=LogNorm())
    plt.xlabel('depth (micron)', ha='right', x=1.0)
    plt.ylabel('energy (keV)', ha='right', y=1.0)
    plt.ylim(300,1000)
    plt.xlim(0,60)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('counts')
    plt.tight_layout()
    plt.show()


def resonance_positions():

    if(len(sys.argv) != 2):
        print('Usage: plotting.py [input .hdf5 file (with extension)]')
        sys.exit()

    with open("data.json") as f:
        data = json.load(f)
    data_dir = os.path.expandvars(data["data_dir"])

    # read in pandas dataframe, and pull out data of interest
    df =  pd.read_hdf("{}/{}".format(data_dir,sys.argv[1]), key="procdf")
    df['x'] = df['x']*1000
    df['y'] = df['y']*1000
    df['z'] = df['z']*1000+7.5

    plt.hist(df['z'], np.arange(-1,15,0.1), histtype='step', color = 'black')
    plt.xlabel('z depth (micron)', ha='right', x=1.0)
    plt.ylabel('counts', ha='right', y=1.0)
    #plt.ylim(300,1000)
    #plt.xlim(0,60)
    plt.tight_layout()
    plt.show()


def spectrum():

    if(len(sys.argv) != 2):
        print('Usage: plotting.py [input .hdf5 file (with extension)]')
        sys.exit()

    with open("data.json") as f:
        data = json.load(f)
    data_dir = os.path.expandvars(data["data_dir"])

    # read in pandas dataframe, and pull out data of interest
    df =  pd.read_hdf("{}/{}".format(data_dir,sys.argv[1]), key="procdf")
    df['KE'] = df['KE']*1000

    plt.hist(df['KE'], np.arange(0,1000,1), histtype='step', color = 'black')
    plt.xlabel('capture KE (keV)', ha='right', x=1.0)
    plt.ylabel('counts', ha='right', y=1.0)
    #plt.ylim(300,1000)
    #plt.xlim(0,60)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
	main()

