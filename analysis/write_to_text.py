import h5py
import sys
import pandas as pd
import json, os

if(len(sys.argv) != 2):
    print('Usage: write_to_text.py [input .hdf5 file (with extension)]')
    sys.exit()

with open("data.json") as f:
    data = json.load(f)
data_dir = os.path.expandvars(data["data_dir"])

# read in pandas dataframe
df =  pd.read_hdf("{}/{}".format(data_dir,sys.argv[1]), key="procdf")

print(df)

x_positions_file = open("X17_xpositions.txt", "w")
for i in range(len(df)):
    x_positions_file.write("{}\n".format(df['x'][i]))
x_positions_file.close()

y_positions_file = open("X17_ypositions.txt", "w")
for i in range(len(df)):
    y_positions_file.write("{}\n".format(df['y'][i]))
y_positions_file.close()

z_positions_file = open("X17_zpositions.txt", "w")
for i in range(len(df)):
    z_positions_file.write("{}\n".format(df['z'][i]))
z_positions_file.close()
