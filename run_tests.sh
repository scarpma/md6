#!/bin/bash
make
for dir in test_*/; do
  echo "==========================="
  echo "Doing test ${dir}"
  echo "==========================="
  echo ""
  ./exec/main $dir && (cd $dir && python3 analyze.py)
done
