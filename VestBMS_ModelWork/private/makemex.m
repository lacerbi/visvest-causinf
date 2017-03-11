% COMPILE ALL MEX FILES

files = dir('*.c*');
for iFile = 1:numel(files); mex(files(iFile).name); end