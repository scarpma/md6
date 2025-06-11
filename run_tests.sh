#!/bin/bash
make
for dir in test_*/; do
  echo "==========================="
  echo "Doing test ${dir}"
  echo "==========================="
  echo ""
  if [ -f $dir/"param.in" ]; then
    ./exec/main $dir && (cd $dir && python3 analyze.py)
  else
    ./${dir}/test.sh
  fi
done
