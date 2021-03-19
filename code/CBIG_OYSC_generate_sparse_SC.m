function [output] = CBIG__OYSC_generate_sparse_SC(individual_level_SC, proportion)
% CBIG__OYSC_generate_sparse_SC(masked_individual_level_SC, proportion)
% 
% This function is to generate sparse SC matrix by selecting edges that 
% appear in a specific proportion of subjects. 
% For group level SC, the averaging only considers the subjects with non-zero 
% entries.
% 
% Input: 
%      - individual_level_SC:
%        Containing individual level structral connectivity (SC) of a group of 
%        subjects
%      
%      - proportion: 
%        proportion for selecting edges, defalut = 1.0
% 
% 
% Output: 
%      - group_level_SC:
%      
%      - masked_individual_level_SC: 
%        sparse individual level SC
%        
%      - SC_mask: 
%        mask containing selected edges for the group
% Written by Wulan Naren and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% setting
if isempty(proportion) 
   proportion = 1;
end

[roi_r,roi_c,num] = size(individual_level_SC);
num_edges = roi_r*roi_c;
masked_individual_level_SC_vector = zeros(num_edges,num);
group_level_SC_vector = zeros(num_edges,1);
individual_level_SC_vector = reshape(individual_level_SC,[num_edges,num]);

%%----------------------------------------------------------------------------
% make statistics for each edge, the number of the edge appearing in all subjects  
zero_count_map = zeros(num_edges, 1);
for i = 1:num_edges
    zero_count_map(i,1) = sum(individual_level_SC_vector(i,:) ~= 0);
end

SC_mask_vector = zero_count_map;
% only select edges whose ratio of appearance is larger than a proportion
% obtain sparse individual level SC
SC_mask_vector(SC_mask_vector(:,1) < round(proportion*num),1) = 0; 
SC_mask_vector(SC_mask_vector(:,1) ~= 0,1) = 1; 


for i = 1:num
   masked_individual_level_SC_vector(:, i) = SC_mask_vector.*individual_level_SC_vector(:, i);
end

% calculate group level SC
for i = 1:num_edges
    group_level_SC_vector(i) = mean(nonzeros(masked_individual_level_SC_vector(i, :)));
end
group_level_SC_vector(isnan(group_level_SC_vector)) = 0;
%%----------------------------------------------------------------------------

% output
output.group_level_SC = reshape(group_level_SC_vector, [roi_r,roi_c]);
output.masked_individual_level_SC = reshape(masked_individual_level_SC_vector, [roi_r,roi_c,num]);
output.SC_mask = reshape(SC_mask_vector, [roi_r,roi_c]);

end
