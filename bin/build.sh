#!/usr/bin/env bash

# Find the absolute path to the directory containing this script
BIN_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Change into that directory so relative paths work perfectly
cd "${BIN_DIR}"

if [ -d HMM/build ]; then
  rm -rf HMM/build
fi

cd HMM
mkdir build
cd build
cmake ..
make -j 4