function CBIG_OYSC_unpair2SampletTest_with_FDR_example()
% CBIG_OYSC_unpair2SampletTest_with_FDR_example()
% 
% This function is the example function of the unpaired two sample t-test with FDR correction. 
% 
% Written by Wulan Naren and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% set path
curPath = mfilename('fullpath');
path_prefix = extractBefore(curPath,'/nwulan_research_code');
savepath_prefix = extractBefore(curPath,'/code');
data_path = '/data/Data_HelenGroup_elderlySC_youngSC/log_norm_SC';
savepath = fullfile(savepath_prefix, '/results/');

% set parameter
is_spectralRadius = 1; % if tuning structural connectivity matrix with spectral radius
alpha = 1; % dynamics of structural connectivity matrix
significance = 0.05; % significance level for FDR

% load data
load(fullfile(path_prefix, data_path)); 
dataSC.group1_individual_level = old_SC_ind;
dataSC.group2_individual_level = young_SC_ind;
dataSC.group1_group_level = old_SC_avg;
dataSC.group2_group_level = young_SC_avg;

% t-test function
output = CBIG_OYSC_unpair2SampletTest_with_FDR(dataSC, significance, savepath, is_spectralRadius, alpha);
l_sigedges = length(output.significant_index)
l_wedge = length(output.weaker_connection_index)
l_sedges = length(output.stronger_connection_index)

end

