#!/bin/sh

# At this time this script is not able to read $(MAKE) from the test harness;
# we are assuming 'make' will do the job

make serialtests
cd ../tests-solvers
make serialtests
cd ../tests-std-framework
make serialtests
cd ../tests-std-mesh
make serialtests
cd ../tests-utils
make serialtests
cd ../tests-core
