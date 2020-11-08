#!/bin/bash
#PJM -L "rscunit=ito-a" 
#PJM -L "rscgrp=ito-ss"
#PJM -L "vnode=1"
#PJM -L "vnode-core=36"
#PJM -L "elapse=3:00:00"
#PJM -j
#PJM -X

ulimit -s unlimited
module load intel/2019.4
#module load openmpi/3.1.3-cuda9.1-intel18.3

# Open MPI
#export PATH="${HOME}/local/openmpi/4.0.4/bin:$PATH"
#export LD_LIBRARY_PATH="${HOME}/local/openmpi/4.0.4/lib:$LD_LIBRARY_PATH"

NUM_NODES=${PJM_VNODES}
NUM_CORES=36
NUM_PROCS=36

export I_MPI_PERHOST=$NUM_CORES
export I_MPI_FABRICS=shm:ofi

export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/bin/pjrsh
export I_MPI_HYDRA_HOST_FILE=${PJM_O_NODEINF}

mpiexec.hydra -n $NUM_PROCS ./real.exe
mpiexec.hydra -n $NUM_PROCS ./wrf.exe namelist.input
