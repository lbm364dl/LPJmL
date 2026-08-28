#!/bin/bash
# Build and run the divide_cells() unit test against the current libraries.
set -e
cd "$(dirname "$0")/.."
mpicc -Wall -m64 -DUSE_RAND48 -DUSE_MPI -DSAFE -DUSE_NETCDF -DUSE_UDUNITS \
      -DPERMUTE -DSTRICT_JSON -O2 -Iinclude -o /tmp/test_divide bench/test_divide.c \
      lib/liblpj.a lib/libtools.a lib/libnum.a lib/libbstruct.a lib/libtools.a \
      -lm -lnetcdf -ludunits2 -ljson-c
exec /tmp/test_divide
