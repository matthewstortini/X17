import numpy as np
import h5py
import pandas as pd
import sys
import time
from numba import jit
from tqdm import tqdm
import json, os
import random
#np.set_printoptions(threshold=sys.maxsize)

def main():

    #true_angles()
    #proton_depths()
    capture_positions()
    #pair_energy_loss()

def true_angles():
    ''' created to pull out angles between e+e- when running in x17 mode. this computes the actual angle between the particles, not
        necessarily the angle between them that may be measured (for instance if one scatters before being measured). '''

    if(len(sys.argv) != 3):
        print('Usage: postprocesshdf5.py [input filename.hdf5 (with extension)] [output filename.hdf5 (with extension)]')
        sys.exit()

    start = time.time()
    print('In Progress...')

    pd.options.mode.chained_assignment = None

    # set data directory where data is stored
    with open("data.json") as f:
        data = json.load(f)
    data_dir = os.path.expandvars(data["data_dir"])

    # import data to be read
    g4sfile = h5py.File('{}/{}'.format(data_dir,sys.argv[1]), 'r')
    g4sntuple = g4sfile['default_ntuples']['g4sntuple']

    # Taking data from g4sntuple and organizing it into a pandas dataframe.
    g4sdf = pd.DataFrame(np.array(g4sntuple['event']['pages']), columns=['event'])
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['pdx']['pages']),
                       columns=['pdx']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['pdy']['pages']),
                       columns=['pdy']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['pdz']['pages']),
                       columns=['pdz']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['step']['pages']),
                       columns=['step']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['pid']['pages']),
                       columns=['pid']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['parentID']['pages']),
                       columns=['parentID']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['KE']['pages']),
                       columns=['KE']), lsuffix = '_caller', rsuffix = '_other')
    
    df = g4sdf.loc[(g4sdf.parentID==0)&(g4sdf.step==0)&(g4sdf.pid!=0)&(g4sdf.pid!=22)]

    procdf= pd.DataFrame(df.groupby(['event'], as_index=False)['pdx','pdy','pdz'].prod())

    pdx_array = procdf['pdx'].values
    pdy_array = procdf['pdy'].values
    pdz_array = procdf['pdz'].values

    theta_array = np.arccos(pdx_array+pdy_array+pdz_array)*180/np.pi

    procdf['theta'] = theta_array
    procdf.drop(columns=['pdx', 'pdy','pdz'])

    # Save the pandas dataframe.
    procdf.to_hdf('{}/{}'.format(data_dir,sys.argv[2]), key='procdf', mode='w')


def proton_depths():
    ''' this function was created to track the proton when running in gun mode to see how far it penetrates into the
        geometries one chooses to track in the macro file when running the simulation. '''

    if(len(sys.argv) != 3):
        print('Usage: postprocesshdf5.py [input filename.hdf5 (with extension)] [output filename.hdf5 (with extension)]')
        sys.exit()

    start = time.time()
    print('In Progress...')

    pd.options.mode.chained_assignment = None

    # set data directory where data is stored
    with open("data.json") as f:
        data = json.load(f)
    data_dir = os.path.expandvars(data["data_dir"])

    # import data to be read
    g4sfile = h5py.File('{}/{}'.format(data_dir,sys.argv[1]), 'r')
    g4sntuple = g4sfile['default_ntuples']['g4sntuple']

    # Taking data from g4sntuple and organizing it into a pandas dataframe.
    g4sdf = pd.DataFrame(np.array(g4sntuple['event']['pages']), columns=['event'])
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['step']['pages']),
                       columns=['step']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['pid']['pages']),
                       columns=['pid']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['z']['pages']),
                       columns=['z']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['KE']['pages']),
                       columns=['KE']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['volID']['pages']),
                       columns=['volID']), lsuffix = '_caller', rsuffix = '_other')
    
    procdf = g4sdf.loc[(g4sdf.pid==2212)&(g4sdf.step>0)&(g4sdf.KE>0)&(g4sdf.volID!=0)]

    # Save the pandas dataframe.
    procdf.to_hdf('{}/{}'.format(data_dir,sys.argv[2]), key='procdf', mode='w')


def capture_positions():
    ''' for postprocessing data after running in capture mode, with /tracking/recordAllSteps in the macro. the
        purpose of this script is to pull out the position that protons are captured at. this data can then be
        used to set the position of X17s when running in X17 mode. to convert the hdf5 file output here to .txt
        files that the X17 PrimaryGenerator uses, use the script write_to_text.py -- this writes the positions
        to a .txt file and saves them to the correct directory. '''

    if(len(sys.argv) != 3):
        print('Usage: postprocesshdf5.py [input filename.hdf5 (with extension)] [output filename.hdf5 (with extension)]')
        sys.exit()

    start = time.time()
    print('In Progress...')

    pd.options.mode.chained_assignment = None

    # set data directory where data is stored
    with open("data.json") as f:
        data = json.load(f)
    data_dir = os.path.expandvars(data["data_dir"])

    # import data to be read
    g4sfile = h5py.File('{}/{}'.format(data_dir,sys.argv[1]), 'r')
    g4sntuple = g4sfile['default_ntuples']['g4sntuple']

    # Taking data from g4sntuple and organizing it into a pandas dataframe.
    g4sdf = pd.DataFrame(np.array(g4sntuple['event']['pages']), columns=['event'])
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['step']['pages']),
                       columns=['step']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['pid']['pages']),
                       columns=['pid']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['x']['pages']),
                       columns=['x']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['y']['pages']),
                       columns=['y']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['z']['pages']),
                       columns=['z']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['pdx']['pages']),
                       columns=['pdx']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['pdy']['pages']),
                       columns=['pdy']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['pdz']['pages']),
                       columns=['pdz']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['KE']['pages']),
                       columns=['KE']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['mass']['pages']),
                       columns=['mass']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['beta']['pages']),
                       columns=['beta']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['volID']['pages']),
                       columns=['volID']), lsuffix = '_caller', rsuffix = '_other')
    
    # here i only want to look at captures in volume 1, but one can change this if they please
    procdf = g4sdf.loc[(g4sdf.step!=0)&(g4sdf.pid==2212)&(g4sdf.volID==1)]

    # calculate momenta and add to dataframe
    procdf['px'] = procdf['mass']*procdf['beta']*procdf['pdx']/np.sqrt(1-procdf['beta']*procdf['beta'])
    procdf['py'] = procdf['mass']*procdf['beta']*procdf['pdy']/np.sqrt(1-procdf['beta']*procdf['beta'])
    procdf['pz'] = procdf['mass']*procdf['beta']*procdf['pdz']/np.sqrt(1-procdf['beta']*procdf['beta'])

    # drop unnecessary columns
    procdf = procdf.drop(columns=['volID','step','event','pid','pdx','pdy','pdz','mass','beta'])

    print(len(procdf))

    # Save the pandas dataframe.
    procdf.to_hdf('{}/{}'.format(data_dir,sys.argv[2]), key='procdf', mode='w')


def pair_energy_loss():
    ''' this function was created to see how much energy e+e- pair loses in foil. We only keep events for which less than
        100 keV is lost. '''

    if(len(sys.argv) != 3):
        print('Usage: postprocesshdf5.py [input filename.hdf5 (with extension)] [output filename.hdf5 (with extension)]')
        sys.exit()

    start = time.time()
    print('In Progress...')

    pd.options.mode.chained_assignment = None

    # set data directory where data is stored
    with open("data.json") as f:
        data = json.load(f)
    data_dir = os.path.expandvars(data["data_dir"])

    # import data to be read
    g4sfile = h5py.File('{}/{}'.format(data_dir,sys.argv[1]), 'r')
    g4sntuple = g4sfile['default_ntuples']['g4sntuple']

    # Taking data from g4sntuple and organizing it into a pandas dataframe.
    g4sdf = pd.DataFrame(np.array(g4sntuple['event']['pages']), columns=['event'])
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['step']['pages']),
                       columns=['step']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['pid']['pages']),
                       columns=['pid']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['Edep']['pages']),
                       columns=['Edep']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['volID']['pages']),
                       columns=['volID']), lsuffix = '_caller', rsuffix = '_other')

    df = g4sdf.loc[(g4sdf.volID<=1)]

    procdf = pd.DataFrame(df.groupby(['event'], as_index=False)['Edep'].sum())

    procdf = procdf.loc[(procdf.Edep<0.100)]

    print(len(procdf))

    # Save the pandas dataframe.
    procdf.to_hdf('{}/{}'.format(data_dir,sys.argv[2]), key='procdf', mode='w')


if __name__ == '__main__':
        main()
