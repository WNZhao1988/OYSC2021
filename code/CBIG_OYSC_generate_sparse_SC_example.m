function CBIG_OYSC_generate_sparse_SC_example()
% CBIG_OYSC_generate_sparse_SC_example()
% 
% This function is the example function to generate sparse SC. 
% 
% Written by Wulan Naren and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% set path
curPath = mfilename('fullpath');
path_prefix = extractBefore(curPath,'/nwulan_research_code');
data_path = '/data/Data_HelenGroup_elderlySC_youngSC/log_norm_SC';
savepath = fullfile(path_prefix, '/data/Data_HelenGroup_elderlySC_youngSC/');

% set parameter
proportion = 1; 
% controling the ratio for selecting edges, edges who appear in at least 
% a proportion of subjects are retained 

% load data
load(fullfile(path_prefix, data_path)); 

[output_o] = CBIG_OYSC_generate_sparse_SC(old_SC_ind, proportion);
[output_y] = CBIG_OYSC_generate_sparse_SC(young_SC_ind, proportion);

save([savepath,'old_sparse_SC_ind_group.mat'], 'output_o');
save([savepath,'young_sparse_SC_ind_group.mat'],'output_y');
end