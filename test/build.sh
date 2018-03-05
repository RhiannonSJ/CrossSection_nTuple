#!/bin/bash

export FHICL_DIR=/sbnd/app/users/rsjones/LArSoft_v06_69_00/LArSoft-v06_70_00/srcs/analysistree/analysistree
cd $MRB_BUILDDIR
mrbsetenv
mrb i -j4
mrbslp
cd -
export FHICL_FILE_PATH=${FHICL_DIR}:${FHICL_FILE_PATH}
