
clear;

% % plot group-level FC of elderly and young
% curPath = mfilename('fullpath');
% pathPrefix = extractBefore(curPath, '/nwulan_research_code');
% dataPath= '/data/Data_HelenGroup_elderlySC_youngSC/FC_average';
% load(fullfile(pathPrefix, dataPath));
% 
% ch_plot_TY400_fc(old_fc_avg, 'TY126');
% ch_plot_TY400_fc(young_fc_avg, 'TY126');

% % plot sparse group-level SC of elderly and young
% load('/mnt/isilon/CSC1/Yeolab/Users/nwulan/data/Data_HelenGroup_elderlySC_youngSC/young_sparse_SC_ind_group.mat');
% ch_plot_TY400_fc(output_y.group_level_SC, 'TY126');
% saveas(gcf, '../results/sparse_group_level_SC_young.png');
% load('/mnt/isilon/CSC1/Yeolab/Users/nwulan/data/Data_HelenGroup_elderlySC_youngSC/old_sparse_SC_ind_group.mat');
% ch_plot_TY400_fc(output_o.group_level_SC, 'TY126');
% saveas(gcf, '../results/sparse_group_level_SC_old.png');

%% unpaired two sample t-test for sparse SC matrix 
% set path
curPath = mfilename('fullpath');
path_prefix = extractBefore(curPath,'/nwulan_research_code');
savepath_prefix = extractBefore(curPath,'/code');
data_path_o = '/data/Data_HelenGroup_elderlySC_youngSC/old_sparse_SC_ind_group.mat';
data_path_y = '/data/Data_HelenGroup_elderlySC_youngSC/young_sparse_SC_ind_group.mat';
savepath = fullfile(savepath_prefix, '/results/sparse_');

% set parameter
is_spectralRadius = 1; % if tuning structural connectivity matrix with spectral radius
alpha = 1; % dynamics of structural connectivity matrix
significance = 0.05; % significance level for FDR
cbar = [-10,10];

% load data
load(fullfile(path_prefix, data_path_o)); 
load(fullfile(path_prefix, data_path_y)); 
dataSC.group1_individual_level = output_o.masked_individual_level_SC;
dataSC.group2_individual_level = output_y.masked_individual_level_SC;
dataSC.group1_group_level = output_o.group_level_SC;
dataSC.group2_group_level = output_y.group_level_SC;

% t-test function
output = CBIG_OYSC_unpair2SampletTest_with_FDR(dataSC, significance, savepath, is_spectralRadius, alpha, cbar);
l_sigedges = length(output.significant_index)
l_wedge = length(output.weaker_connection_index)
l_sedges = length(output.stronger_connection_index)
