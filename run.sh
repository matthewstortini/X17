# make SIMULATIONNAME equal to what you want to name your simulation file
SIMULATIONNAME="mysimulation"

sed "s/GenericName/$SIMULATIONNAME/" run.mac > Run1.mac
./X17 Run1.mac
mv $SIMULATIONNAME.hdf5 $X17DATADIR
rm Run1.mac
