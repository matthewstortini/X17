import numpy as np
import h5py
import pandas as pd
import sys
import time
from numba import jit
from tqdm import tqdm
import json, os
#np.set_printoptions(threshold=sys.maxsize)

def main():

    #process_angles()
    #process_depth()
    resonance_positions()

def process_angles():

    if(len(sys.argv) != 3):
        print('Usage: postprocesshdf5.py [input filename.hdf5 (with extension)] [output filename.hdf5 (with extension)]')
        sys.exit()

    start = time.time()
    print('In Progress...')

    pd.options.mode.chained_assignment = None

    # Import g4simple data
    with open("data.json") as f:
        data = json.load(f)
    data_dir = os.path.expandvars(data["data_dir"])

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
    procdf.to_hdf('{}'.format(data_dir,sys.argv[2]), key='procdf', mode='w')


def process_depth():

    if(len(sys.argv) != 3):
        print('Usage: postprocesshdf5.py [input filename.hdf5 (with extension)] [output filename.hdf5 (with extension)]')
        sys.exit()

    start = time.time()
    print('In Progress...')

    pd.options.mode.chained_assignment = None

    # Import g4simple data
    with open("data.json") as f:
        data = json.load(f)
    data_dir = os.path.expandvars(data["data_dir"])

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
    procdf.to_hdf('{}'.format(data_dir,sys.argv[2]), key='procdf', mode='w')


def resonance_positions():

    if(len(sys.argv) != 3):
        print('Usage: postprocesshdf5.py [input filename.hdf5 (with extension)] [output filename.hdf5 (with extension)]')
        sys.exit()

    start = time.time()
    print('In Progress...')

    pd.options.mode.chained_assignment = None

    # Import g4simple data
    with open("data.json") as f:
        data = json.load(f)
    data_dir = os.path.expandvars(data["data_dir"])

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
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['KE']['pages']),
                       columns=['KE']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['volID']['pages']),
                       columns=['volID']), lsuffix = '_caller', rsuffix = '_other')
   
    df = g4sdf.loc[(g4sdf.pid==2212)&(g4sdf.volID==2)&(g4sdf.KE>0.4409)&(g4sdf.KE<0.4411)]
    procdf = pd.DataFrame(df.groupby(['event'], as_index=False)['x','y','z'].mean())

    # Save the pandas dataframe.
    procdf.to_hdf('{}/{}'.format(data_dir,sys.argv[2]), key='procdf', mode='w')


if __name__ == '__main__':
        main()
