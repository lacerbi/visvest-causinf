#!/bin/bash
module purge
module load matlab
export MATLABPATH=${MATLABPATH}:/home/la67/MATLAB
matlab -nodisplay
