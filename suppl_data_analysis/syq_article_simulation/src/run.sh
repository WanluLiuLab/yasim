#!/usr/bin/env bash
set -uex

git clone http://github.com/yukiteruono/pbsim3
cd pbsim3 || exit 1
git checkout f3e5f1cf3d0e8346b5e4598ac238b2b570b223e8
autoreconf -i
./configure
make -j20
ln -sf "$(pwd)"/src/pbsim ../../bin/pbsim3
cd .. || exit 1

git clone http://github.com/yukiteruono/pbsim2
cd pbsim2 || exit 1
git checkout eeb5a19420534a0f672c81db2670117e62a9ee38
autoreconf -i
./configure
make -j20
ln -sf "$(pwd)"/src/pbsim ../../bin/pbsim2
cd .. || exit 1

git clone https://salsa.debian.org/med-team/pbsim # Network issues
cd pbsim || exit 1
git checkout a7a3a39c774b240d482e1f444bc40cefdd93227f
autoreconf -i
./configure
make -j20
ln -sf "$(pwd)"/src/pbsim ../../bin/pbsim
cd .. || exit 1
