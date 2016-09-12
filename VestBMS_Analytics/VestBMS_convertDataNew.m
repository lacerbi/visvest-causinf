%VESTBMS_CONVERTDATA Convert Kalpana's data in a CueBMS-compatible format.

% Description of data matrix
% X(:, 1) Absolute trial number
% X(:, 2) Empty
% X(:, 3) Visual noise level (1 Low, 2 Med, 3 High)
% X(:, 4) Trial response type (1 Visual, 2 Vestibular, 3 Categorical)
% X(:, 5) Vestibular stimulus position
% X(:, 6) Visual stimulus position
% X(:, 7) Number of stimuli (1 Unimodal, 2 Bimodal)
% X(:, 8) Vestibular response
% X(:, 9) Visual response
% X(:, 10) Categorical response

function [data,bigdata,X] = VestBMS_convertDataNew(flag)

if nargin < 1 || isempty(flag); flag = 1; end

if flag     % Use new (corrected) datasets
    filesuffix = '_causalinf_leftright';
    monkeyinverted = 0;
else
    filesuffix = '_cause_leftright';    
    monkeyinverted = 1; % Monkey data were inverted
end

% All subjects:
% 9 and 14 have only 2 levels of visual reliability (instead of 3)
% 5 and 14 only have no bimodal localization trials
% 7 and 15 are missing

subjs = [1:4 6 8 10:13 16 24 44 45];
ishuman = [1:4 6 8 10:13 16];
ismonkey = [24 44 45];
% bincenters = -45:2.5:45;
bincenters = [-45,-40,-35,-30:2.5:-2.5,-1.25:0.625:1.25,2.5:2.5:30,35,40,45];

for iSubj = 1:length(subjs)
    nid = subjs(iSubj)
    
    % Load unimodal data
    filename = [num2str(nid) '_leftright.mat'];
    if ~exist(filename,'file'); filename = [num2str(nid) filesuffix '.mat']; end
    temp = load(filename);
    if isfield(temp,'dataset')  % We like inconsistent naming
        dataset = temp.dataset;
    else
        dataset = temp.comp_dataset;
    end
    D = dataset(dataset(:,1) == 1 | dataset(:,1) == 2, :);
    
    % Load bimodal localization data
    filename = [num2str(nid) filesuffix '.mat'];
    temp = load(filename);
    if isfield(temp,'dataset')
        dataset = temp.dataset;
    else
        dataset = temp.comp_dataset;
    end
    
    % Load bimodal unity judgments data (human only)
    if any(nid == ishuman)
        filename = [num2str(nid) '_causalinf_cause.mat'];
        temp = load(filename);
        temp.dataset(:,1) = 3;  % Unity judgement bisensory data
        dataset = [dataset; temp.dataset];
    end
    
    D = [D; dataset(dataset(:,1) == 3 | dataset(:,1) == 4, :)];
    
    % Remove missed trials (0 timeout, 4 trial start button, 5 unknown error)
    D(D(:, 5) == 0 | D(:, 5) == 4 | D(:, 5) == 5 | D(:, 9) == 0 | D(:, 9) == 4 | D(:, 9) == 5,:) = [];    
        
    % Remove out of range trials (only subject 10)
    idx = abs(D(:,3) - D(:,2)) == 50; 
    D(idx, :) = [];
    if sum(idx) > 0; fprintf('Removing %d out of range trials for subject %d.\n', sum(idx), iSubj); end
    
    % Correct inverted assignement of responses in monkey bisensory data
    if monkeyinverted && any(nid == ismonkey)
        D(D(:,1) == 4, 5) = 3 - D(D(:,1) == 4, 5);
    end    
    
    % Create empty datas matrix
    nTrials = size(D,1);
    X{iSubj} = zeros(nTrials, 10);    
    X{iSubj}(:, 1) = 1:nTrials;             % Trial number
    X{iSubj}(:, 3) = D(:, 7);               % Visual noise level
    
    temp = zeros(nTrials,1);
    temp(D(:, 1) == 1) = 2;                 % Vestibular-only trials
    temp(D(:, 1) == 2) = 1;                 % Visual-only trials
    temp(D(:, 1) == 3) = 3;                 % Categorization trials
    temp(D(:, 1) == 4) = 2;                 % Vestibular bimodal trials
    
    if flag
        X{iSubj}(:, 2) = D(:,10);           % Session number
    end
    
    X{iSubj}(:, 4) = temp;                  % Trial response type
    
    X{iSubj}(:, 5) = D(:, 2);               % Vestibular stimulus position
    X{iSubj}(:, 6) = D(:, 3);               % Visual stimulus position
    X{iSubj}(:, 7) = 1 + (D(:, 1) >= 3);    % Number of stimuli
    
    D(D(:, 5) == 1, 5) = -1;                % Leftward response
    D(D(:, 5) == 2, 5) = 1;                 % Rightward response
    
    X{iSubj}(:,[8 9]) = NaN;
    f = X{iSubj}(:, 4) == 2;                % Vestibular trials
    X{iSubj}(f,8) = D(f,5);                 % Vestibular localization response
    f = X{iSubj}(:, 4) == 1;                % Visual trials
    X{iSubj}(f,9) = D(f,5);                 % Visual localization response
    
    X{iSubj}(:, 10) = D(:, 9);              % Categorical response
    
    
end

% Analytics options
options.robustfitflag = 1;
options.quickplotflag = 0;
options.bincenters = bincenters;
options.psycholeftright = 1;
options.psycholeftrightdelta = 1;
options.bindata = 1;
options.flatten = 1;

data = VestBMS_analytics(X,options);
for i = 1:length(data); data{i}.id = i; end
display([num2str(length(data)), ' datasets successfully converted.']);

% All users data
allX = [];
for i = 1:length(data); allX = [allX; X{i}]; end

if nargout > 1
    bigdata = VestBMS_analytics({allX},options);
    display('Mean subject dataset converted.');
end


end