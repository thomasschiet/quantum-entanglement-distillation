#!/bin/sh

RC=1
while [ $RC -ne 0 ]; do
   julia generateData.jl
   RC=$?
done
