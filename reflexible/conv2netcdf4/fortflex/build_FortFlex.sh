#!/bin/bash

rm -f FortFlex.pyf
/opt/anaconda/bin/f2py -m FortFlex -h FortFlex.pyf FortFlex.f
/opt/anaconda/bin/f2py -c --fcompiler=gfortran FortFlex.pyf FortFlex.f
mv FortFlex.so ..

