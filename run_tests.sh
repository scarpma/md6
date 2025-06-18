#!/bin/bash
make
for dir in test_*/; do
  echo "==========================="
  echo "Doing test ${dir}"
  echo "==========================="
  echo ""
  ./exec/main $dir > $dir/log.txt && (cd $dir && python3 analyze.py)
done
