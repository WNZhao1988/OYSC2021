function [corr_mat, major_grid2_out, major_grid1_out] = ch_plot_TY400_fc(corr_mat, parcellation, drop_rois, isSubplot, subplotVec, myPost,cl, figName, saveDir, clmap)
%  [corr_mat major_grid2_out major_grid1_out] = ch_plot_TY400_fc(corr_mat, parcellation, drop_rois, isSubplot, supplotVec, myPost, cl, figName, saveDir, clmap)
%  
%  Updated version of CH's and Jesisca's FC matrix plotting scripts.  
%  
%  Output 
%       corr_mat: Reordered fc matrix
%       major_grid2_out, major_grid1_out: Ploting parameters (lines separating major and sub networks)
%         
%  Inputs
%       corr_mat: Numeric array. The fc matrix. Can be the full matrix (a warning will be issued), 
%           or the matrix with invalid ROIs already removed. 
%        
%           A sample_fc.mat is accompanied for testing
% 
%       parcellation: String. Keyword that calls the default ROI-fc matrix structure is set within in the script
%           using my_grid().
%           Optional Strings are: 'TY126', 'TY144', 'TY400', and 'TY430'
% 
%       drop_rois: Numeric (row) vector. Indices of invalid ROIs. FC matrix and gridlines 
%           will be adjusted after dropping rows and columns according to
%           drop_rois. Unless -2 is included, the order of drop_rois should
%           be consistent with the raw FC matrix (i.e., usually both are
%           BEFORE any reordering).
%
%           Set to -1 to export new roi orders and nothing else (no reordering of FC matrix, no dropping roi allowed) 
%           Set to [-1 n1 n2 ... nZ] to export new roi orders with dropping roi allowed
%           Note 1: corr_mat can be empty
%           Note 2: A usually harmless warning is issued for -1. 
%           Note 3: A variable 'roi_order' will be generated in the Workspace
%           Set to -2 to disable reordering
%           Set to [-2 n1 n2 ... nZ] to diable reordering but allow ROIs to
%           be dropped. In this case the FC matrix and n1...nZ are usually
%           AFTER some prior reordering already.
%
%       isSubplot: Numeric. Simple shorthand for plotting multiple matrices into
%           subplot(). e.g., 133 stands for subplot(133). Note it is very
%           restrictive subplotting.
%           
%           0 means no subplot
%           -1 means not plotting at all. Useful for creating the reordered fc matrix
%
%       subplotVec (beta): Numeric (row) vector. If you want something more than subplot(123), 
%           provide the subplot vector here. The last digit in isSubplot will be ignored. 
%           Omit it or give [] to skip.
%
%       myPost: 4-element numeric vector. If you have specific [left bottom width height] to use, add it here.
%               Set it to 0 or [] to skip.
%
%       cl: Numeric vector [a b]; colorbar limits. If empty or not provided, the limits are 
%           derived automatically. Can be skipped.
%
%       figName, saveDir: Inputs for saving the output as .fig and .tiff (require matlab2015 or 
%           later); Pass with empty string '', or skip if not needed.
%     
%       clmap: custom color map. The one used by the lab and past publications is now embeded
%           in the script. Can be skipped.
%     
%  Note:
%  *  To use other parcellations/ ROIs, just add new definitions in my_grid(),
%     and provide the right parcellation and drop_rois when calling the function. 
%     color map is also embedded into the script
% 
%  ** Please use matlab2015 or above if you want valid TIFF output
% 
%  *** Written based on the ordering of
%  /mnt/HZLHD2/Data7/Members/Eric/rsFC_2017/data/Schaefer2018_400ROIParcellation_MNI152_2mm
%  (Older; shared by TY personally, uploaded by QX
%  Note that its order is the same but the img may (or may not) differ from
%  /mnt/HZLHD2/Data7/Members/Eric/rsFC_2017/data/Schaefer2018_400ROIParcellation_MNI152_2mm_reordered
%  (Newer; downloaded from GitHub)
%
%  History:
%  ch_plot_TY126(corr_mat, figName, saveDir)
%  based on TY126plot_new_order_js
%  use js order and colorset(gcf,'Colormap',rbmap2)
%  include reorder
%  Includes only cleaned reorder (after removing unwanted ROIs)
% 
%  Modified by Eric for 400 ROIs (no subcortical) on 25 Sep 2017
%  Modified for 430 ROIs (new ThomasYeo + subcortical) sometime in Nov 2017
%%
FontSize=10;
FontName='Arial';

% Grid Lines; OHBM abstract is .1, .5, and 1
cellGrid=0.0; % original 0.1
thinGrid=1; % original 1
majorGrid=3; % original 3

isEye = 1; % turn diagonals to 0 (1 = yes; 0 = no; e.g., ISC)?
% Overridden if all diagonals have the same value

% misc
ti = [0];% if you want to reduce white space; not very efficient
cpos = [];% size of colorbar; not very efficient
% myColorLabel = 'Bootstrap ratio'; % label of color bar
myColorLabel = 'Connectivity Strength'; % label of color bar
%%
if ~exist('drop_rois','var') || isempty(drop_rois)
   % drop_rois = [38 40 42 43 81  100 178  179  180,  235, ...
   %     237   238   239   240   241   242   243   280   281   284   299   302   303   380   381]; % SWMi rest-task reconfig: task
drop_rois = [0];
% elseif 
%    % drop_rois = [38 40 42 43 81  100 178  179  180,  235, ...
%    %     237   238   239   240   241   242   243   280   281   284   299   302   303   380   381]; % SWMi rest-task reconfig: task
% drop_rois = [0];
end

if ~exist('clmap','var')    
    clmap='rbmap2';
end
    

if ~exist('figName','var')
    figName='TY400 Sample Plot';
end

if exist('saveDir','var')
    
    if ~isempty(saveDir) && ~exist(saveDir,'dir')
        mkdir(saveDir);
    end

    fileName=[saveDir regexprep(figName,' ','_')];
end

% shortcut to get the roi_reorder
if numel(corr_mat) == 1 & corr_mat == -1
    my_grid(parcellation, -1);
    return
end

% load color scheme
% if exist('/mnt/HZLServer/Data2/Program_HZL/Scripts_CH/functions/plots/corr_mat_colorscale_bluemod.mat','file')
%     load('/mnt/HZLServer/Data2/Program_HZL/Scripts_CH/functions/plots/corr_mat_colorscale_bluemod.mat')
% elseif exist('X:/Data2/Program_HZL/Scripts_CH/functions/plots/corr_mat_colorscale_bluemod.mat');
%     load('X:/Data2/Program_HZL/Scripts_CH/functions/plots/corr_mat_colorscale_bluemod.mat')
% else
%     load('/Volumes/HZLServer/Data2/Program_HZL/Scripts_CH/functions/plots/corr_mat_colorscale_bluemod.mat')   
% end
%

% Provide the color map inside this function
rbmap2 = my_rbmap(); % empty = scale is always black at zero, 1 = auto scale

% load reorder
% if exist('/mnt/HZLServer/Data2/Program_HZL/Scripts_CH/functions/plots/fileOrder_to_js.mat','file')
%     load('/mnt/HZLServer/Data2/Program_HZL/Scripts_CH/functions/plots/fileOrder_to_js.mat')
% else
%     load('X:/Data2/Program_HZL/Scripts_CH/functions/plots/fileOrder_to_js.mat')
% end
% if exist('Z:\Data6\SYNC\Stats\SYNC_TCRP\FCScripts\fileOrder_to_142_cleaned.mat')
%     load('Z:\Data6\SYNC\Stats\SYNC_TCRP\FCScripts\fileOrder_to_142_cleaned.mat',...
%     'new_order_js','major_grid1','major_grid2','major_grid','major_label'); 
% else
%     load('/mnt/HZLHD2/Data6/DCFC_2015/Stats/DCFC2015_DFCwj/EEG_fMRI/cor_EEG_fMRI/scripts/fileOrder_to_js_cleaned.mat',...
%     'new_order_js','major_grid1','major_grid2','major_grid','major_label'); 
% end
%% Adjust order, labels and grids given drop rois, no effect on the FC matrix
try
    [major_grid1, major_grid2, major_grid, major_label, new_order_js] = my_grid(parcellation, drop_rois); % All outputs are ADJUSTED for omitted ROIs
catch
    warning('Variables ''fc_grids'', ''roi_order'' created; are you asking only to generate roi_order?');
    return
end
major_grid2_out = major_grid2; % output grid info for plotting outside this function
major_grid1_out = major_grid1;
% trim fc matrix and ROI information to fit the project-specific fc matrix (fewer than
% default nROIs)
% new_order_TY_sc=1:size(corr_mat,1);
% new_order_TY_sc(1:length(new_order_js))=new_order_js(1:length(new_order_js));

% reorder fc matrix according to plotting scheme
%[corr_mat]=ch_reorder(new_order_TY_sc,corr_mat);
%[corr_mat]=ch_reorder(new_order_js,corr_mat);

% replace by this?
%% Apply the new order, labels, and grids to the FC matrix
[corr_mat] = eric_reorder(new_order_js, corr_mat, drop_rois);

if isEye
    corr_mat(logical(eye(size(corr_mat,1)))) = 0;
end

if ~exist('isSubplot','var')
    isSubplot = 0;
end

% Plot corr_mat using imagesc
if isSubplot ~= -1 % skip plot if -1
    if ~isSubplot
        h1=figure; clf;
    else
        plotStr = num2str(isSubplot);
        if ~exist('subplotVec','var') || isempty(subplotVec)
            if strcmpi(plotStr(3),'1')
                h1 = figure; clf;
            else
                h1 = gcf;
            end
            subplot(isSubplot);
        else % if you want to stack with other plots (beta version)
            if ~ismember(1, subplotVec) % assume figure already created
                h1 = gcf;
            else
                h1 = figure; clf;
            end
            subplot(str2num(plotStr(1)),str2num(plotStr(2)), subplotVec);
        end
    end
    
    imagesc(corr_mat);
  
    
    if ~isempty(rbmap2) % do not want to use Matlab auto scale
        bb=['set(gcf,''Colormap'',' clmap ');'];
        eval(bb);
    end
    
    % set(gcf,'Colormap',rbmap2);
    % set(gcf,'Colormap',Cool);
    
    % Generate thin grid lines (white)
    [xline, yline, ymaj] = generateline(size(corr_mat,1));
    if cellGrid > 0
        patch(xline, yline,'w', 'edgecolor', 'w', 'Linewidth', cellGrid, 'EdgeAlpha', 0.2);
        patch(yline, xline,'w', 'edgecolor', 'w', 'Linewidth', cellGrid, 'EdgeAlpha', 0.2);
    end
    
    % new: minimize white space based on the current position?
    ax = gca;
    outerpos = get(ax,'OuterPosition');
    if isempty(ti)
        ti = get(ax,'TightInset');
    end
    if numel(ti) == 4
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        set(ax,'Position', [left bottom ax_width ax_height]);
    end
    
    % set plot position
    xlim(gca,[1 size(corr_mat, 1)]);
    ylim(gca,[1 size(corr_mat, 1)]);
    postn = get(gca, 'Position');
    postn(2) = 0.15; % shift it up
    if exist('myPost','var') && numel(myPost) == 4
        postn = myPost; % manually define position and size of matrix
    end
    set(gca, 'Position', postn);
   
    
    % Set colorbar
    isPlot = 0;
    if ~isSubplot
        colorbarWhere = 'side';
        isPlot = 1;
    else
        colorbarWhere = 'bottom';
        if ~exist('subplotVec','var') || isempty(subplotVec)
            if strcmpi(plotStr(2),plotStr(3)) % plot color bar only when it is the last panel
                isPlot = 1;
            end
        else
            colorbarWhere = 'side';
            isPlot = 1; % when a vector is provided, just plot
        end
    end
    if isPlot && strcmpi(colorbarWhere, 'bottom')
        hcol=colorbar('peer',gca,'SouthOutside');
        cpos=get(hcol,'Position');
        cpos(4)=cpos(4)/4; % Halve the thickness
        cpos(3)=cpos(3)*0.75; % Reduce length
        cpos(1)=cpos(1) + 0.05; % + .1 Move it to the center
        cpos(2)=cpos(2) -.05;%+ 0.12; % Move it down outside the plot
        set(hcol,'Position',cpos);
        set(get(hcol, 'xlabel'),'String',myColorLabel);
    elseif isPlot && strcmpi(colorbarWhere, 'side')
        hcol=colorbar('peer',gca,'EastOutside');
        if isempty(cpos)
            cpos=get(hcol,'Position');
        end
        if numel(cpos) == 4
            % flip the dimensions
            cpos(3)=cpos(3)/4; % scale thickness
            old4 = cpos(4);
            cpos(4)=cpos(4)*.8; % Reduce length
            cpos(1)= max(postn(1) + postn(3)*.85,0);%,postn(3)+postn(1)+.01); % Move it to the side more
            cpos(2)=cpos(2)+abs(old4-cpos(4))/2; % Move it up; keep it center
            set(hcol,'Position',cpos);
        end
        set(get(hcol, 'ylabel'),'String',myColorLabel);
    end
    
    % Set color limit
    collim = max(max(abs(corr_mat)));
    scalelim = [-1*collim, 1*collim];
    
    if scalelim(1)==scalelim(2)
        scalelim(2)=scalelim(2)+1;
    end
    
    
    % corr_mat((corr_mat==0))=nan;
    % max(corr_mat(:))
    % min(corr_mat(:))
    % scalelim = [min(corr_mat(:)), max(corr_mat(:))];
    
    % scalelim=[-10,10]; %%%
    % scalelim=[3,15]; % color scale for without GSR, RSFC plot
    % scalelim=[-1,1];
    
    % override auto color limit
    if exist('cl','var') && ~isempty(cl)
        scalelim=cl;
    end
    
    set(gca, 'CLim', scalelim);
    
    % make it square
    axis equal;
    grid off;
    axis([-5 size(corr_mat, 1)+5.5 -5 size(corr_mat, 1)+5.5]);
    set(gca, 'Visible', 'off');
    set(gcf, 'color', 'white');
    
    
    % this is thin line, for sub-networks
    %major_grid1 = [6  15  28 39 67 86 100];
    
    % this is thick line, for major networks
    %major_grid2 = [ 24  50  54  78  92  102  108  112  114];
    
    % add lines of networks
    patch(xline(:,major_grid1), yline(:,major_grid1),'w', 'edgecolor', 'w', 'Linewidth', thinGrid,'EdgeAlpha', 0.4);
    patch(yline(:,major_grid1), xline(:,major_grid1),'w', 'edgecolor', 'w', 'Linewidth', thinGrid,'EdgeAlpha', 0.4);
    
    patch(xline(:,major_grid2), yline(:,major_grid2),'w', 'edgecolor', 'w', 'Linewidth', majorGrid,'EdgeAlpha', 0.8);
    patch(yline(:,major_grid2), xline(:,major_grid2),'w', 'edgecolor', 'w', 'Linewidth', majorGrid,'EdgeAlpha', 0.8);
    
    
    % location of label, the first label is for block 0-12, 2nd for 12-24, last
    % is for 124-126
    %major_grid = [6 15 24 28 39 50 54 67 78 86 92 100 102 108 112 114];
    %major_label= {'Default C','Default B','Default A','Cont C', 'Cont B','Cont A','Limbic', ' SalVenAttn B', ' SalVenAttn A','DorsAttn B','DorsAttn A','SomMot B','SomMot A','VisPeri','VisCent', 'TempPar','Subcortical'};
    
    x_label=major_label;
    y_label=major_label;
    
    major_grid1b=[0, major_grid];
    major_grid2b=[major_grid, size(corr_mat,1)];
    for i2=1:length(major_label);
        label_location(i2)=ceil((major_grid2b(i2)-major_grid1b(i2))/2)+major_grid1b(i2);
    end
    
    xtick=label_location; % location of x label
    ytick=label_location;
    
    tx=text(xtick,(size(corr_mat,1)+1.5)*ones(1,length(xtick)),x_label);
    set(tx,'HorizontalAlignment','right','VerticalAlignment','top', 'Rotation', 90, 'FontSize',FontSize,'FontWeight','bold');
    set(gca,'XTickLabel','','XTick',xtick);
    
    yx=text((-0.5)*ones(1,length(ytick)),ytick,y_label);
    set(yx,'HorizontalAlignment','right','VerticalAlignment','middle', 'Rotation', 0, 'FontSize',FontSize,'FontWeight','bold');
    set(gca,'yTickLabel','','YTick',ytick);
    
    set(gca, 'fontname', FontName, 'fontweight','bold');
    
    
    % fine adjustment for saving figures
    set(gcf, 'position' ,[10 10 900 600]);
    figpos = get(gcf, 'position');
    %set(gcf, 'papersize', figpos(3:4).*100)
    set(gcf, 'papersize', figpos(3:4), 'paperpositionmode','auto','inverthardcopy','off');
    
    if exist('fileName','var')
        if ~isempty(fileName)
            saveas(h1,fileName,'fig');
            print(h1,fileName,'-dtiff','-r200');
            close(gcf);
        end
    end
    
end
end


function [x,y, ymaj] = generateline(n)
x = 1.5:1:n;
x = [ x; x; repmat(nan,1,(n-1)) ];
y = [ 0.5 n+0.5 nan ].';
y = repmat(y,1,(n-1));

ymaj = [ -5 n+5.5 nan]';
% For short grid (same as the figure boundary)
% ymaj = [ 0.5 n+0.5 nan]';
ymaj = repmat(ymaj,1,(n-1));
end

function [ rB ] = ch_reorder( A,B )
% Reorder the matrix

if length(A) ~= size(B, 1)
    error('the matrix size is not right :)')
end


for i=1:length(A)
    for i2=1:length(A)
        rB(i,i2)=B(A(i),A(i2));
    end
end

end

function [corr_mat] = eric_reorder(new_order_js,corr_mat, drop_rois)
if max(new_order_js) > size(corr_mat,1)
    % expand corr_mat
    tmp_mat = []; n = 1;
    for i = 1:max(new_order_js) % rows first
        if ~ismember(i, drop_rois)
            tmp_mat(i,:) = corr_mat(n,:); % assign rows to the expanded matrix
            n = n + 1;
        else
            tmp_mat(i,1:size(corr_mat,2)) = nan;
        end
    end
    for j = 1:max(new_order_js) % insert columns
        if ismember(j, drop_rois)
            % split the matrix, insert a column
            t1 = tmp_mat(:,1:j-1);
            t3 = tmp_mat(:,j:end);
            t2 = nan(size(tmp_mat,1),1);
            tmp_mat = [t1 t2 t3];
        end
    end
    corr_mat = tmp_mat(new_order_js, new_order_js); %
elseif max(new_order_js) == size(corr_mat,1) % no ROI dropped, i.e., equal number of indices and matrix rows/columns
    corr_mat = corr_mat(new_order_js, new_order_js);
else
    error('Input functional matrix is smaller than expected');
end
if numel(unique(diag(corr_mat))) == 1 % e.g. ISC would have meaningful diagonals
    for i = 1:size(corr_mat,1) % remove diagnonal
        corr_mat(i,i) = nan;
    end
end
end

function [rbmap2] = my_rbmap(varargin)
if nargin == 0
    % directly provided color map array (always black at the middle of the scale) instead of reading a separate file
    rbmap2 = [ 0.6000    1.0000    1.0000
        0.5250    1.0000    1.0000
        0.4500    1.0000    1.0000
        0.3750    1.0000    1.0000
        0.3000    1.0000    1.0000
        0.2250    1.0000    1.0000
        0.1500    1.0000    1.0000
        0.0750    1.0000    1.0000
        0    1.0000    1.0000
        0    0.9375    1.0000
        0    0.8750    1.0000
        0    0.8125    1.0000
        0    0.7500    1.0000
        0    0.6875    1.0000
        0    0.6250    1.0000
        0    0.5625    1.0000
        0    0.5000    1.0000
        0    0.4375    1.0000
        0    0.3750    1.0000
        0    0.3125    1.0000
        0    0.2500    1.0000
        0    0.1875    1.0000
        0    0.1250    1.0000
        0    0.0625    1.0000
        0         0    1.0000
        0         0    0.9333
        0         0    0.8667
        0         0    0.8000
        0         0    0.7333
        0         0    0.6667
        0         0    0.6000
        0         0    0.5333
        0         0    0.4667
        0         0    0.4000
        0         0    0.3333
        0         0    0.2667
        0         0    0.2000
        0         0    0.1333
        0         0    0.0667
        0         0         0
        0         0         0
        0.0385         0         0
        0.0769         0         0
        0.1154         0         0
        0.1538         0         0
        0.1923         0         0
        0.2308         0         0
        0.2692         0         0
        0.3077         0         0
        0.3462         0         0
        0.3846         0         0
        0.4231         0         0
        0.4615         0         0
        0.5000         0         0
        0.5385         0         0
        0.5769         0         0
        0.6154         0         0
        0.6538         0         0
        0.6923         0         0
        0.7308         0         0
        0.7692         0         0
        0.8077         0         0
        0.8462         0         0
        0.8846         0         0
        0.9231         0         0
        0.9615         0         0
        1.0000         0         0
        1.0000    0.1000         0
        1.0000    0.2000         0
        1.0000    0.3000         0
        1.0000    0.4000         0
        1.0000    0.5000         0
        1.0000    0.6000         0
        1.0000    0.7000         0
        1.0000    0.8000         0
        1.0000    0.9000         0
        1.0000    1.0000         0
        1.0000    1.0000    0.2000
        1.0000    1.0000    0.4000
        1.0000    1.0000    0.6000];
else
    rbmap2 = [];
end
end

function varargout = my_grid(nROIs, drop_rois)
% step 1 defines full ROI matrix
switch nROIs
    case 'TY126'
        major_grid1 = [6 15  28 39 67 86  100  ];
        major_grid2 = [24 50 54 78 92 102 108 112 114]; % omit the last = nROI
        major_grid = sort([major_grid1 major_grid2]); %16/17 networks
        network_labels = {'Default C', 'Default B', 'Default A', 'Cont C', 'Cont B', 'Cont A', 'Limbic', ...
            ' SalVentAttn B', ' SalVentAttn A', 'DorsAttn B', 'DorsAttn A', 'SomMot B', 'SomMot A', 'VisPeri', ...
            'VisCent', 'TempPar', 'Subcortical'};
        % reordered ROI index
        roi_order = [[54:56 111:113],[49:53 107:110],[45:48 102:106],[43:44 100:101],[37:42 95:99],[31:36 90:94], ...
            [29:30 88:89],...
            [23:28 81:87],[18:22 75:80],[14:17 71:74],[11:13 68:70],[7:10 64:67],[6 63], ...
            [3:5 60:62],[1:2 58:59],[57 114],[115:2:125 116:2:126]];
    case 'TY144'
        % no roi_reorder because ROI has been ordered in the mask files;
        % may as well make a dummy structure here in the future.
        % check.
        % Revised for full 144 ROIs on 18Jan2018
        major_grid1 = [6 15 28 39 67 86 100 108]; % major_grid1; subnetworks thin lines
        major_grid2 = [24 50 54 78 92 102 112 114]; % major_grid2; major network thick lines
        major_grid = sort([major_grid1 major_grid2]); % major_grid 7 networks x 2 hemispheres?
        network_labels  = { 'Default C', 'Default B', 'Default A', 'Cont C', 'Cont B', 'Cont A', 'Limbic', ...
            ' SalVentAttn B', ' SalVentAttn A', 'DorsAttn B', 'DorsAttn A', 'SomMot B', 'SomMot A', 'VisPeri', ...
            'VisCent', 'TempPar', 'Subcortical'};
        roi_order = [ 54    55    56   111   112   113    49    50    51    52    53   107   108   109   110    45    46    47    48 ...
            102   103   104   105   106    43    44   100   101    37    38    39    40    41    42    95    96    97    98 ...
            99    31    32    33    34    35    36    90    91    92    93    94    29    30    88    89    23    24    25 ...
            26    27    28    81    82    83    84    85    86    87    18    19    20    21    22    75    76    77    78 ...
            79    80    14    15    16    17    71    72    73    74    11    12    13    68    69    70     7     8     9 ... 
            10    64    65    66    67     6    63     3     4     5    60    61    62     1     2    58    59    57   114 ...
            115   116   117   118   119   120   121   122   123   124   125   126   127   128   129   130   131   132   133 ...
            134   135   136   137   138   139   140   141   142   143   144]; 
    case 'TY400'
        % grid line positions; not ROI index, but vector index
        %major_grid1 = [7 34 64 84 101 129 152 172 196 228 253 283 317 349 372 390]; % networks A B C thin lines
        %major_grid2 = [13 45 79 91 116 140 164 181 215 240 267 298 337 360 384 400]; % main network thick lines
        major_grid1 = [13 45 91 116 181 240 298];
        major_grid2 = [79 140 164 215 267 337 360 384]; % omit the last = nROI
        major_grid = sort([major_grid1 major_grid2]); %16/17 networks
        network_labels = {'Default C', 'Default B', 'Default A', 'Cont C', 'Cont B', 'Cont A', 'Limbic', ...
            ' SalVentAttn B', ' SalVentAttn A', 'DorsAttn B', 'DorsAttn A', 'SomMot B', 'SomMot A', 'VisPeri', ...
            'VisCent', 'TempPar'};
        % reordered ROI index
        roi_order = [ 188   189   190   191   192   193   194   385   386   387   388   389   390   167   168   169   170   171   172 ...
            173   174   175   176   177   178   179   180   181   182   183   184   185   186   187   374   375   376   377 ...
            378   379   380   381   382   383   384   149   150   151   152   153   154   155   156   157   158   159   160 ...
            161   162   163   164   165   166   358   359   360   361   362   363   364   365   366   367   368   369   370 ...
            371   372   373   144   145   146   147   148   351   352   353   354   355   356   357   134   135   136   137 ...
            138   139   140   141   142   143   336   337   338   339   340   341   342   343   344   345   346   347   348 ...
            349   350   121   122   123   124   125   126   127   128   129   130   131   132   133   325   326   327   328 ...
            329   330   331   332   333   334   335   109   110   111   112   113   114   115   116   117   118   119   120 ...
            313   314   315   316   317   318   319   320   321   322   323   324   101   102   103   104   105   106   107 ...
            108   304   305   306   307   308   309   310   311   312    86    87    88    89    90    91    92    93    94 ...
            95    96    97    98    99   100   285   286   287   288   289   290   291   292   293   294   295   296   297 ...
            298   299   300   301   302   303    73    74    75    76    77    78    79    80    81    82    83    84    85 ...
            273   274   275   276   277   278   279   280   281   282   283   284    60    61    62    63    64    65    66 ...
            67    68    69    70    71    72   259   260   261   262   263   264   265   266   267   268   269   270   271 ...
            272    44    45    46    47    48    49    50    51    52    53    54    55    56    57    58    59   244   245 ...
            246   247   248   249   250   251   252   253   254   255   256   257   258    25    26    27    28    29    30 ...
            31    32    33    34    35    36    37    38    39    40    41    42    43   224   225   226   227   228   229 ...
            230   231   232   233   234   235   236   237   238   239   240   241   242   243    13    14    15    16    17 ...
            18    19    20    21    22    23    24   213   214   215   216   217   218   219   220   221   222   223     1 ...
            2     3     4     5     6     7     8     9    10    11    12   201   202   203   204   205   206   207   208 ...
            209   210   211   212   195   196   197   198   199   200   391   392   393   394   395   396   397   398   399 ...
            400];
    case 'TY430' % subcortical ROIs appended
        % grid line positions; not ROI index, but vector index
        %major_grid1 = [7 34 64 84 101 129 152 172 196 228 253 283 317 349 372 390]; % networks A B C thin lines
        %major_grid2 = [13 45 79 91 116 140 164 181 215 240 267 298 337 360 384 400]; % main network thick lines
        major_grid1 = [13 45 91 116 181 240 298];
        major_grid2 = [79 140 164 215 267 337 360 384 400]; % omit the last = nROI
        major_grid = sort([major_grid1 major_grid2]); %16/17 networks
        network_labels = {'Default C', 'Default B', 'Default A', 'Cont C', 'Cont B', 'Cont A', 'Limbic', ...
            ' SalVentAttn B', ' SalVentAttn A', 'DorsAttn B', 'DorsAttn A', 'SomMot B', 'SomMot A', 'VisPeri', ...
            'VisCent', 'TempPar', 'Subcortical'};
        % reordered ROI index
        roi_order = [ 188   189   190   191   192   193   194   385   386   387   388   389   390   167   168   169   170   171   172 ...
            173   174   175   176   177   178   179   180   181   182   183   184   185   186   187   374   375   376   377 ...
            378   379   380   381   382   383   384   149   150   151   152   153   154   155   156   157   158   159   160 ...
            161   162   163   164   165   166   358   359   360   361   362   363   364   365   366   367   368   369   370 ...
            371   372   373   144   145   146   147   148   351   352   353   354   355   356   357   134   135   136   137 ...
            138   139   140   141   142   143   336   337   338   339   340   341   342   343   344   345   346   347   348 ...
            349   350   121   122   123   124   125   126   127   128   129   130   131   132   133   325   326   327   328 ...
            329   330   331   332   333   334   335   109   110   111   112   113   114   115   116   117   118   119   120 ...
            313   314   315   316   317   318   319   320   321   322   323   324   101   102   103   104   105   106   107 ...
            108   304   305   306   307   308   309   310   311   312    86    87    88    89    90    91    92    93    94 ...
            95    96    97    98    99   100   285   286   287   288   289   290   291   292   293   294   295   296   297 ...
            298   299   300   301   302   303    73    74    75    76    77    78    79    80    81    82    83    84    85 ...
            273   274   275   276   277   278   279   280   281   282   283   284    60    61    62    63    64    65    66 ...
            67    68    69    70    71    72   259   260   261   262   263   264   265   266   267   268   269   270   271 ...
            272    44    45    46    47    48    49    50    51    52    53    54    55    56    57    58    59   244   245 ...
            246   247   248   249   250   251   252   253   254   255   256   257   258    25    26    27    28    29    30 ...
            31    32    33    34    35    36    37    38    39    40    41    42    43   224   225   226   227   228   229 ...
            230   231   232   233   234   235   236   237   238   239   240   241   242   243    13    14    15    16    17 ...
            18    19    20    21    22    23    24   213   214   215   216   217   218   219   220   221   222   223     1 ...
            2     3     4     5     6     7     8     9    10    11    12   201   202   203   204   205   206   207   208 ...
            209   210   211   212   195   196   197   198   199   200   391   392   393   394   395   396   397   398   399 ...
            400 401:430];
end

if numel(drop_rois) == 1 && drop_rois == -1 % rigid; no dropping of rois
    assignin('base','roi_order',roi_order);
    return
elseif numel(drop_rois) > 1 && drop_rois(1) == -1% flexible, allow dropping rois, alternative to eric_reorder
    idx = ~ismember(roi_order, drop_rois(2:end)); % get rois that are not dropped, still with their original order ID
    roi_order2 = roi_order(idx); % drop rois
    
    % Twisted way to get the adjusted labels (could have merged this
    % section with the other if loop below
    % convert ROI idx in drop_rois into vector index
    drop_rois_vec_idx = [];
    for idrop = drop_rois(:)'
        drop_rois_vec_idx = [drop_rois_vec_idx find(roi_order == idrop)]; % record the new positions
    end
    drop_rois_vec_idx = sort(drop_rois_vec_idx);
    % shift grid lines accordingly; update grid1 and grid2; shift by the number
    % of indices smaller than the grid line index
    adj1 = []; adj2 = [];
    for igrid = 1:numel(major_grid1) % these are vector indices
        adj1 = [adj1 sum(major_grid1(igrid) >= drop_rois_vec_idx)];
    end
    for igrid = 1:numel(major_grid2) % these are vector indices
        adj2 = [adj2 sum(major_grid2(igrid) >= drop_rois_vec_idx)];
    end
    fc_grids.major_grid1= [major_grid1 - adj1];
    fc_grids.major_grid2 = [major_grid2 - adj2];
    fc_grids.major_grid = sort([fc_grids.major_grid1, fc_grids.major_grid2]);
    fc_grids.network_labels = network_labels; % n
    assignin('base','roi_order',roi_order2);
    assignin('base','fc_grids',fc_grids);
    return
end

if numel(drop_rois) == 1 && drop_rois == -2 % no reordering, no ROI dropping
    varargout(1) = {major_grid1};
    varargout(2) = {major_grid2};
    varargout(3) = {sort([varargout{1}, varargout{2}])};
    varargout(4) = {network_labels};
    varargout(5) = {1:numel(roi_order)}; % no reordering, so it is just 1,2,...n
else
    if numel(drop_rois) > 1 && drop_rois(1) == -2 % no reordering, allow ROI dropping
        roi_order = 1:numel(roi_order); % no reordering
        drop_rois = drop_rois(2:end);
    end
    % step 2: trimmed ROIs, adjust grid indices
    idx = ~ismember(roi_order, drop_rois); % find out the new position of the to-be-dropped ROIs
    % convert ROI idx in drop_rois into vector index
    drop_rois_vec_idx = [];
    for idrop = drop_rois(:)'
        drop_rois_vec_idx = [drop_rois_vec_idx find(roi_order == idrop)]; % record the new positions
    end
    drop_rois_vec_idx = sort(drop_rois_vec_idx);
    % shift grid lines accordingly; update grid1 and grid2; shift by the number
    % of indices smaller than the grid line index
    adj1 = []; adj2 = [];
    for igrid = 1:numel(major_grid1) % these are vector indices
        adj1 = [adj1 sum(major_grid1(igrid) >= drop_rois_vec_idx)];
    end
    for igrid = 1:numel(major_grid2) % these are vector indices
        adj2 = [adj2 sum(major_grid2(igrid) >= drop_rois_vec_idx)];
    end
    varargout(1) = {[major_grid1 - adj1]};
    varargout(2) = {[major_grid2 - adj2]};
    varargout(3) = {sort([varargout{1}, varargout{2}])};
    varargout(4) = {network_labels}; % network labels usually don't change; unless a whole network is dropped (does not work for this right now)
    varargout(5) = {roi_order(idx)}; % remove omitted ROIs from roi_order
end
end

