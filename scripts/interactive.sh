#!/bin/bash
qsub  -lmem=8192mb -lnodes=1:ppn=4 -lwalltime=4:00:00 -qinteractive -I  -M${USER}@nyu.edu
