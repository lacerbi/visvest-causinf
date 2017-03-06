#!/bin/bash
module purge
module load matlab
export MATLABPATH=${MATLABPATH}:${HOME}/MATLAB
matlab -nodisplay
