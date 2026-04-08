#!/bin/bash
set -e

# Make sure we are in the package root
if [ ! -d "tests/testthat/testdata" ]; then
  echo "Error: Please run this script from the root directory of the package."
  exit 1
fi

TMP_DIR=$(mktemp -d)
echo "Created temporary directory: $TMP_DIR"

echo "Starting R build script..."
Rscript inst/exec/build_testdata.R "$TMP_DIR"

echo "Copying files to tests/testthat/testdata/STAT3"
rm -rf tests/testthat/testdata/STAT3/*
cp -r "$TMP_DIR/STAT3/"* tests/testthat/testdata/STAT3/

echo "Copying files to tests/testthat/testdata/CTNNB1"
rm -rf tests/testthat/testdata/CTNNB1/*
cp -r "$TMP_DIR/CTNNB1/"* tests/testthat/testdata/CTNNB1/

echo "Cleaning up temporary directory..."
rm -rf "$TMP_DIR"

echo "Done rebuilding test cases!"
