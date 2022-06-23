#!/usr/bin/bash

sed -i 's/COFF/OFF/' $1/*.off
awk -i inplace '{printf $1=="" ? "\n" : "%s %s %s %s", $1, $2, $3, $4 ;  printf $5=="" ? "\n" : " %.3g %.3g %.3g\n", $5 / 255, $6 / 255, $7 / 255 }' outputs/*.off
