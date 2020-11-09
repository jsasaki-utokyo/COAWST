#!/bin/bash
#PJM -L "rscunit=ito-a" 
#PJM -L "rscgrp=ito-single"
#PJM -L "vnode=1"
#PJM -L "vnode-core=1"
#PJM -L "elapse=3:00:00"
#PJM -j
#PJM -X

module load intel/2019.4
./coawstS < Projects/Sandy/ocean_sandy.in

