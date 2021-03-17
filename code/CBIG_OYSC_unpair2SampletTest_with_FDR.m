function [output] = CBIG_OYSC_unpair2SampletTest_with_FDR(dataSC, significance, savepath, is_spectralRadius, alpha)
% CBIG_OYSC_unpair2SampletTest_with_FDR(data, is_spectralRadius, alpha)
% 
% This function is to calculate the unpaired two sample t-test with FDR correction. 
% 
% Input: 
%      - dataSC:
%        Containing individual level structral connectivity (SC) of group 1, 
%        individual level structral connectivity of group 2
%        group level structral connectivity of group 1
%        group level structral connectivity of group 2
%      
%      - significance: 
%        Significance level for FDR, defalut = 0.05
% 
%      - savepath:
%        Specifying saving path for results.
% 
%      - is_spectralRadius:
%        If tuning structural connectivity matrix with spectral radius, defalut = 0.
%        
%      - alpha: 
%        Dynamics of structural connectivity matrix, defalut = 0.
% 
% Output: 
%      - significant_index:
%        index of significant edges
%      
%      - weaker_connection_index: 
%        index of weaker connections in significant edges
% 
%      - strongerer_connection_index: 
%        index of stronger connections in significant edges
%        
%      - alpha: 
%        Dynamics of structural connectivity matrix, defalut = 0.
% Written by Wulan Naren and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% setting
if isempty(significance) 
   significance = 0.05;
end

if isempty(savepath) 
   savepath = '../results/';
end

if isempty(is_spectralRadius) 
   is_spectralRadius = 0;
end

if isempty(alpha) 
   alpha = 1;
end


[roi_r,roi_c,num1] = size(dataSC.group1_individual_level);
[roi_r,roi_c,num2] = size(dataSC.group2_individual_level);
num_edge = roi_r*roi_c;
P_original_vector = zeros(num_edge,1); % original p-value vector
group1_SC_ind_level_vector = reshape(dataSC.group1_individual_level, [num_edge,num1]); % convert matrix to vector
group2_SC_ind_level_vector = reshape(dataSC.group2_individual_level, [num_edge,num2]);
weaker_connection_index = [];  % store index for edges whose connection weight in group 1 is smaller that in group
stronger_connection_index = [];

%% tuning SC matrix with spectral radius (largest eigenvalue) and adjust it by dynamic factor alpha
if is_spectralRadius == 1
    group1_SC_ind_level_eigenvalue_tuning = zeros(roi_r,roi_c,num1);
    group2_SC_ind_level_eigenvalue_tuning = zeros(roi_r,roi_c,num2);

    for i = 1 : num1
       e1 = eig(dataSC.group1_individual_level(:,:,i));
       ed1 = sort(e1,'descend');
       group1_SC_ind_level_eigenvalue_tuning(:,:,i) = alpha*dataSC.group1_individual_level(:,:,i)./ed1(1);
    end

    for i = 1:  num2
       e2 = eig(dataSC.group2_individual_level(:,:,i));
       ed2 = sort(e2,'descend');
       group2_SC_ind_level_eigenvalue_tuning(:,:,i) = alpha*dataSC.group2_individual_level(:,:,i)./ed2(1);
    end
    
    clear e2 e1 ed2 ed1

    group1_SC_ind_level_vector = reshape(group1_SC_ind_level_eigenvalue_tuning, [num_edge,num1]);
    group2_SC_ind_level_vector = reshape(group2_SC_ind_level_eigenvalue_tuning, [num_edge,num2]);

end
%% unpaired two sample t-test for each edge of structural connectivity in 
% group 1 and the corresponding edge in group 2 -> p-value matrix -> 
% FRD correction -> -log10 -> display


% ttest
for i = 1:num_edge
   [h,p] = ttest2(group1_SC_ind_level_vector(i,:), group2_SC_ind_level_vector(i,:)); % >
   P_original_vector(i) = p;
end

% convert vector of lower triangular matrix
temp_matrix = reshape(P_original_vector, [roi_r,roi_c]);
temp_matrix = tril(temp_matrix);
mask  = tril(true(size(temp_matrix))); % 0, 1 logical mask matrix
P_original__lowtri_vector = temp_matrix(mask).';  % A' and A.' same for real number matrices and arrays
clear mask

% FDR
[significant_index, final_threshold] = FDR(P_original__lowtri_vector, significance);   % FDR correction 
% for Bonferroni correction -> select 400+ significant edges, too strong     

% convert vector to symmetrix matrix and perform log function
P_original__lowtri_vector_FDR = P_original__lowtri_vector;
P_original__lowtri_vector_FDR(setdiff(1:end,significant_index)) = 1; % insignificant edges
P_log__lowtri_vector_FDR = -log10(P_original__lowtri_vector_FDR(:));

%% among significant different edges calculated above, find SC edges of old people so that weight of these edges is 
% smaller and great than that of young people 

%smaller
group1_SC_group_level_lowtri = tril(dataSC.group1_group_level);
mask1  = tril(true(size(group1_SC_group_level_lowtri))); % 0, 1 logical mask matrix
group1_SC_group_level_lowtri_vector = group1_SC_group_level_lowtri(mask1).';

group2_SC_group_level_lowtri = tril(dataSC.group2_group_level);
mask2  = tril(true(size(group2_SC_group_level_lowtri))); % 0, 1 logical mask matrix
group2_SC_group_level_lowtri_vector = group2_SC_group_level_lowtri(mask2).';
clear mask1 mask2

for i=1:length(significant_index)
    if group1_SC_group_level_lowtri_vector(significant_index(i)) < group2_SC_group_level_lowtri_vector(significant_index(i))
       weaker_connection_index = [weaker_connection_index; significant_index(i)]; 
    else group1_SC_group_level_lowtri_vector(significant_index(i)) > group2_SC_group_level_lowtri_vector(significant_index(i))
       stronger_connection_index = [stronger_connection_index; significant_index(i)];  
    end
end


P_log__lowtri_vector_FDR(weaker_connection_index) = -P_log__lowtri_vector_FDR(weaker_connection_index);
P_log_full = tril(ones(roi_r));  
P_log_full(P_log_full == 1) = P_log__lowtri_vector_FDR; % convert to full matrix (N by N)
P_log_full_sym = P_log_full + P_log_full'; % convert to symmtric matrix

%% visualization
[P_log_full_sym_reorder] = ch_plot_TY400_fc(P_log_full_sym, 'TY126');
caxis([-10 10]);
axs2 = struct(gca);
cb2 = axs2.Colorbar;
cb2.Label.String = ' ';
set(gcf,'color', 'w');
saveas(gcf,[savepath 'log10(p-value).png']);

% histogram 
figure(4);
weak_connection_edge = P_original__lowtri_vector_FDR(1:end,weaker_connection_index);
histogram(weak_connection_edge);
set(gcf,'color', 'w');
saveas(gcf,[savepath 'histogram_weakerweight.png']);

% histogram 
figure(5);
strong_connection_edge = P_original__lowtri_vector_FDR(1:end,stronger_connection_index);
histogram(strong_connection_edge);
set(gcf,'color', 'w');
saveas(gcf,[savepath 'histogram_strongweight.png']);
% heatmap(P_log_full,'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ", 'xlabel','x','ylabel','x','title',sprintf('x'));

%% output
output.significant_index = significant_index;
output.weaker_connection_index = weaker_connection_index;
output.stronger_connection_index = stronger_connection_index;
end

