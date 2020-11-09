#!/bin/bash
#PJM -L "rscunit=ito-a" 
#PJM -L "rscgrp=ito-qq"
#PJM -L "vnode=1"
#PJM -L "vnode-core=3"
#PJM -L "elapse=3:00:00"
#PJM -j
#PJM -X

module load intel/2019.4
export WRF_EM_CORE=1
. ./coawst.bash -j 3 -nocleanwrf
