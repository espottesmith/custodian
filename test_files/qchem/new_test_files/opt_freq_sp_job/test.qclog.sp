You are running Q-Chem version: 5.0.2

#
# job setting
#
local host:  nid04222
current dir: /scratch1/scratchdirs/ewcss/with_water
input file:  /scratch1/scratchdirs/ewcss/with_water/49_18490/pro_84283.in
output file: /scratch1/scratchdirs/ewcss/with_water/49_18490/pro_84283.out
nprocs     : 1
nthreads   : 24
#
# qchem installation setting
#
QC:          /global/common/edison/software/qchem/5.1beta
QCAUX:       /global/common/edison/software/qchem/5.1beta/qcaux
QCPROG:      /global/common/edison/software/qchem/5.1beta/exe/qcprog.exe
QCPROG_S:    /global/common/edison/software/qchem/5.1beta/exe/qcprog.exe
PARALLEL:    -DPARALLEL
QCMPI:       slurm
#
# qchem directory setting
#
qcrun:       qchem3444
QCSCRATCH:   /dev/shm/qcscratch/
QCLOCALSCR:  /tmp
QCTMPDIR:    /tmp
QCFILEPREF:  /tmp/qchem3444
QCSAVEDIR:   
workdirs:    /tmp/qchem3444
workdir0:    /tmp/qchem3444
partmpdirs =  
#
# parallel setting
#
QCRSH:           ssh
QCMPI:           slurm
QCMPIRUN:        srun
QCMACHINEFILE:   

#
# env setting
#
exported envs:   QC QCAUX QCSCRATCH QCRUNNAME QCFILEPREF QCPROG QCPROG_S GUIFILE
remove work dirs /tmp/qchem3444.0 -- /tmp/qchem3444.0
rm -rf /tmp/qchem3444.0
rm -rf /tmp/qchem3444
