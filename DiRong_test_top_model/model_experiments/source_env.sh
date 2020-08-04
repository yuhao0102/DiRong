#!/bin/sh


MYPATH=$(readlink -f "${BASH_SOURCE[0]}")
MYPATH=$(dirname "$MYPATH")

source $MYPATH/../model_platform/scripts/register_platform.sh

MYPATH=$(readlink -f "${BASH_SOURCE[0]}")
MYPATH=$(dirname "$MYPATH")

source $MYPATH/../inputdata/register_inputdata.sh

module load netcdf/4.4.1-parallel-icc18-wzm
module load lapack/3.8.0-fenggl
module load mpi/intel/18.0.2-thc
module load intel/18.0.2-thc

