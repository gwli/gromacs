TOP=`pwd`/..
exec_proc=$TOP/build/bin/gmx
data=$TOP/examples/protein

cd $data
cmd="$exec_proc mdrun -v -deffnm em-vac"
sp="/raid/tools/SP/SystemProfiler-linux-internal-3.9.3364-05e1bbd/Target-x86_64/x86_64/sp"

rm -fr ./*.QDSTRM
$sp profile --delay=2 --duration=10 -o mgemm.QDSTRM -t cuda,cublas,curand,pthread $cmd
