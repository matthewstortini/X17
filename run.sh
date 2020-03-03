sed 's/GenericName/SimulationName/' run.mac > Run1.mac
./X17 Run1.mac
mv SimulationName.hdf5 $X17DATADIR
rm Run1.mac
