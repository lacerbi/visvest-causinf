function data = readDeWinkel2017data
%READDEWINKEL2017DATA Read and preprocess data from de Winkel et al. 2017.

% Original data can be downloaded from 
% https://doi.org/10.1371/journal.pone.0169676.s008
% 
% Reference:
% De Winkel, K. N., Katliar, M., & Bülthoff, H. H. (2017). Causal inference 
% in multisensory heading estimation. PloS one, 12(1), e0169676.

expFolder = {'Expt I','Expt II'};
nSubjs = [8 9];

fits = load('deWinkel2017originalfits');

% Datafiles of individual participants (pp) of experiment I and II
% 
% The columns in the datafiles hold the following information:
% Column 1. 'stim_nr' -presentation order of the stimuli;
% Column 2. 'modality' -sensory modality: '1' = visual-only, '2' = inertial-only, '3' = visual-inertial;
% Column 3. 'stim_type' -motion profile 1,2,3 correspond to 2,4,6s duration motion profiles;
% Column 4. 'vis_theta' -heading angle of visual stimulus in radians. 0 when there was no inertial stimulus;
% Column 5. 'ine_theta' -heading angle of inertial stimulus in radians. 0 when there was no visual stimulus;
% Column 6. 'response' -heading angle of response, in radians;
% Column 7. 'corrupted' -boolean variable indicating problems with data recording during the experiment.

for iExp = 1:2
    for iSubj = 1:nSubjs(iExp)
        temp = load([expFolder{iExp} filesep 'pp' num2str(iSubj) '.mat']);
        idx = iSubj + nSubjs(1)*(iExp-1);
        data{idx}.Experiment = iExp;
        data{idx}.Subject = iSubj;        
        mat = table2array(temp.data);
        % Remove corrupted trials (we are not going to analyze them)        
        data{idx}.Mat = mat(mat(:,7) == 0,:);
        data{idx}.Fields = temp.data.Properties.VariableNames;
        
        % Read maximum-likelihood fits obtained by deWinkel et al. 2017
        data{idx}.Params = fits.fitsUnisensory(idx,:);        
    end
end

save('deWinkel2017data.mat','data');

end