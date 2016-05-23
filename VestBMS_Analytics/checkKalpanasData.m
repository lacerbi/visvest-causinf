
%CHECKKALPANASDATA Check consistency in Kalpana's datasets.
function checkKalpanasData(nid)


files = {'_leftright','_cause_leftright','_cause'};

for iFile = 1:length(files)
   filename = [num2str(nid) files{iFile} '.mat'];
   tempdata = load(filename);
   dataset = tempdata.dataset;
   
   % Store all trials
   for iType = 1:4
       D{iFile}{iType} = dataset(dataset(:,1) == iType, :);
       D{iFile}{iType}(isnan(D{iFile}{iType})) = Inf;
   end
end

% Check that all trials are equal (order can be different)
for iType = 1:4
   for iFile = 1:length(files)-1
       thisD = D{iFile}{iType};
       if isempty(thisD); continue; end
       for compare = iFile+1:length(files)
           C = D{compare}{iType};
           if isempty(C); continue; end           
           
           % Start comparison
           for iTrial = 1:size(thisD,1)
               % Find equal trial
               f = find(all(bsxfun(@eq, thisD(iTrial,:), C),2),1)
               if isempty(f); error('Datasets do not match!'); end
               C(f, :) = NaN;
           end
           
       end
   end
end

end