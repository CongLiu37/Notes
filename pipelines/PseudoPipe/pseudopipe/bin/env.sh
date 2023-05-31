#!/bin/sha

if [ ! -z "$PSEUDOPIPE_ENV" ]; then source $PSEUDOPIPE_ENV; return; fi

# Pseudopipe configuration
export PSEUDOPIPE_HOME=/home/c/c-liu/Softwares/PseudoPipe/pseudopipe
export pseudopipe=$PSEUDOPIPE_HOME/core/runScripts.py
export genPgeneResult=$PSEUDOPIPE_HOME/ext/genPgeneResult.sh
export genFullAln=$PSEUDOPIPE_HOME/ext/genFullAln.sh
export fastaSplitter=$PSEUDOPIPE_HOME/ext/splitFasta.py
export sqDummy=$PSEUDOPIPE_HOME/ext/sqDummy.py
export sqDummyR=$PSEUDOPIPE_HOME/ext/sqDummy.R
export blastHandler=$PSEUDOPIPE_HOME/core/processBlastOutput.py
export blastHandlerR=$PSEUDOPIPE_HOME/core/processBlastOutput.R
export extractExLoc=$PSEUDOPIPE_HOME/core/extractKPExonLocations.py

# Python configuration
export pythonExec=/apps/free72/python/2.7.3/bin/python

# R configuration
export RExec=/apps/free81/R/4.0.4/bin/Rscript

# Alignment tools configuration
export formatDB=/home/c/c-liu/miniconda3/bin/formatdb
export blastExec=/home/c/c-liu/miniconda3/bin/blastall
export fastaExec=/apps/free72/fasta/35.4.12/bin/tfasty35
