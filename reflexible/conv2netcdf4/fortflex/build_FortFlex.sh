#!/bin/bash

rm -f FortFlex.pyf
f2py -m FortFlex -h FortFlex.pyf FortFlex.f
f2py -c --fcompiler=libgcc FortFlex.pyf FortFlex.f
mv FortFlex.so ..

