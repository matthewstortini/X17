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
x17_dir = os.path.expandvars(data["x17_dir"])

# read in pandas dataframe
df =  pd.read_hdf("{}/{}".format(data_dir,sys.argv[1]))

df = df.reset_index(drop=True)

print("writing x positions")
x_positions_file = open("{}/capture_xpositions.txt".format(x17_dir), "w")
for i in range(len(df)):
    x_positions_file.write("{}\n".format(df['x'][i]))
x_positions_file.close()

print("writing y positions")
y_positions_file = open("{}/capture_ypositions.txt".format(x17_dir), "w")
for i in range(len(df)):
    y_positions_file.write("{}\n".format(df['y'][i]))
y_positions_file.close()

print("writing z positions")
z_positions_file = open("{}/capture_zpositions.txt".format(x17_dir), "w")
for i in range(len(df)):
    z_positions_file.write("{}\n".format(df['z'][i]))
z_positions_file.close()

print("writing x momenta")
x_momenta_file = open("{}/capture_xmomenta.txt".format(x17_dir), "w")
for i in range(len(df)):
    x_momenta_file.write("{}\n".format(df['px'][i]))
x_momenta_file.close()

print("writing y momenta")
y_momenta_file = open("{}/capture_ymomenta.txt".format(x17_dir), "w")
for i in range(len(df)):
    y_momenta_file.write("{}\n".format(df['py'][i]))
y_momenta_file.close()

print("writing z momenta")
z_momenta_file = open("{}/capture_zmomenta.txt".format(x17_dir), "w")
for i in range(len(df)):
    z_momenta_file.write("{}\n".format(df['pz'][i]))
z_momenta_file.close()
