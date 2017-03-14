#!/bin/bash

f2py -m FortFlex -h FortFlex.pyf FortFlex.f
f2py -c --fcompiler=gfortran FortFlex.pyf FortFlex.f
rm -f FortFlex.pyf
mv *.so ..

