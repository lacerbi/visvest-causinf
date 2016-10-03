function string = VestBMS_getModelName(model)
% VESTBMS_GETMODELNAME return the model string.

% Bimodal-data models
switch model(15)
    case 1; string = 'BP';
    case 2; string = 'GBP';
    case 3; string = 'CX';
    case 4; string = 'SCX';
    case 5; string = 'FF';
    case 6; string = 'BPM';
    case 7; string = 'BPT';
    case 8; string = 'BPS';
    case 9; string = 'PF';
end

if model(16) == 6
    string = [string 'F'];
end

if model(11) == 3
    string = [string 'P'];
end

switch model(9)
    case 1
    case 2; string = [string 'u'];
    case 3; string = [string 'U'];
    case 4; string = [string 'd'];
    case 5; string = [string 'D'];
end

switch model(1)
    case 1; string = [string '-C'];
end

switch model(16)
    case 1; % Nothing
    case 4;
    otherwise; string = [string 's'];
end

end