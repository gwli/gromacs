TOP=`pwd`/..
exec_proc=$TOP/build/bin/gmx
data=$TOP/examples/protein

cd $data
cmd="$exec_proc mdrun -v -deffnm em-vac"

cuda-gdb  --args $cmd
