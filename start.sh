#!/bin/bash

set -e

SQUARE=0
TYPE=
INSTANCES=1


usage()
{
cat << EOF
usage: `basename $0` options

This script starts N magics finders

OPTIONS:
  -h  Show this message
  -t  Number of instances (default $INSTANCES)
  -s  Square (default $SQUARE)
  -b  Do for bishop, rook otherwise
EOF
}

anywait()
{
  for pid in "$@"; do
    while kill -0 "$pid"; do
      sleep 1
    done
  done
}

while getopts "ht:s:b" OPTION
do
  case $OPTION in
    h)
      usage
      exit 0
      ;;
    t)
      INSTANCES=$OPTARG
      ;;
    s)
      SQUARE=$OPTARG
      ;;
    b)
      TYPE="-b"
      ;;
    ?)
      usage
      exit 1
      ;;
  esac
done

for (( i=0; i<$INSTANCES; i++ ))
do
  CPUPROFILE=/tmp/magics-$i.prof ./magics -s $SQUARE $TYPE -1 > out-$i.txt &
done;

exit 0
