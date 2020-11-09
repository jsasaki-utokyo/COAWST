#!/bin/bash
#PJM -L "rscunit=ito-a" 
#PJM -L "rscgrp=ito-ss-dbg"
#PJM -L "vnode=1"
#PJM -L "vnode-core=36"
#PJM -L "elapse=1:00:00"
#PJM -j
#PJM -X

module load intel/2019.4

NUM_NODES=${PJM_VNODES}
NUM_CORES=36
NUM_PROCS=36

export I_MPI_PERHOST=$NUM_CORES
export I_MPI_FABRICS=shm:ofi

export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/bin/pjrsh
export I_MPI_HYDRA_HOST_FILE=${PJM_O_NODEINF}

mpiexec.hydra -n $NUM_PROCS ./coawstM Projects/Sandy/ocean_sandy.in

