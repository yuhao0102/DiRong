#/bin/bash

cd model_platform/models/atm/atm_demo/
make
cp atm_demo ../../../../
cd -

cd model_platform/models/ocn/ocn_demo
make
cp ocn_demo ../../../../
cd -
