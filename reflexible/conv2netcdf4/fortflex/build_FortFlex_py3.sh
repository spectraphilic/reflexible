#!/bin/bash

f2py3 -m FortFlex -h FortFlex.pyf FortFlex.f
f2py3 -c --fcompiler=gfortran FortFlex.pyf FortFlex.f
rm -f FortFlex.pyf
mv *.so ..

