#!/bin/bash
echo "-------------------"
echo "Remove output files"
echo "-------------------"


rm *.odb *.odb_f
rm *~
rm *.pyc
rm *.deps
rm *.log
rm core.*

rm *rpy*

rm *.dat *.prt *.sim *.sta *.txt *.com *.msg *.cae *.rec *.inp *.ipm

rm *.jnl abaqus.guiState* *.pdf
rm Job-1.* *.tmp

rm *.SMABulk