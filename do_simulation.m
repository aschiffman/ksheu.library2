% TNF SIMULATION
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 8*60;

global START_TIME END_TIME NUM_CELLS

% Transcription factor vector
% single cell TFs trajectories
% data_name = 'WT_10ngTNF';
% data_name = 'ikbamut_10ngTNF';
% data_name = 'lps_100ng';
% data_name = 'cpg_100nM';
% data_name = 'pic_50ug';

names = {'TNF10ng_762', 'P3CSK4100ng_547','CpG100nM_610', 'LPS100ng_756','PIC50ug_566','aKO_TNF3.3ng', 'aKO_LPS33ng', 'ikbamut_10ngTNF'};
% names = {'aKO_TNF3.3ng', 'aKO_LPS33ng'};
% names = {'WT_10ngTNF', 'ikbamut_10ngTNF'};
%%
for j = 1:length(names)
data_name = char(names(j));
data = load(strcat('F://enhancer_dynamics/nfkb_trajectories_08142019/nfkb_dynamics_',data_name,'.mat'));
% data = load('F://enhancer_dynamics/nfkb_trajectories_20190302/2KO_TNF.mat');
data = cell2struct(struct2cell(data), {'trajectories'});
if istable(data.trajectories)
    data = table2array(data.trajectories);
else
    data = data.trajectories;
end
data_smooth = smoothrows(data);
data_smooth(any(isnan(data_smooth), 2), :) = []; %remove NaN rows

datazero = vec2mat(data_smooth(:, 1), 1);
subtract = repmat(datazero, 1, size(data,2));
data_smooth = data_smooth - subtract; %subtract the first column to normalize
data_smooth(data_smooth<0) = 0; %take neg. values to be 0
maxA = max(data_smooth, [], 2);
[~, index] = sort(maxA);
data_smooth = data_smooth(index, :);

START_TIME =0;
END_TIME = 480;
NUM_CELLS = size(data_smooth,1); 
sim.time = START_TIME:END_TIME;
% tf = data_smooth(5,1:96); %cut to 8hrs
% plot(sim.time*5, interp1( 1:length(tf), tf, START_TIME:END_TIME,'nearest'));
% xlabel('time (min)')
% ylabel('[TF]')	
% ylim([0 inf]);
% f = fit( transpose(time(:,1:96)*5), transpose(tf), 'linearinterp');
% plot( f, transpose(time(:,1:96)*5), transpose(tf) )

% Starting Conditions
    initvalues = zeros(15,1);
    initvalues(1,1) = 1;    %E0
    initvalues(2,1) = 0;    %E1

    initvalues(3,1) = 0;   %E2
    initvalues(4,1) = 0;   %E3

    initvalues(5,1) = 0;   %E4
    initvalues(6,1) = 0;   %E5
    initvalues(7,1) = 0;   %E6
    initvalues(8,1) = 0;   %E7
    initvalues(9,1) = 0;   %E8
    initvalues(10,1) = 0;   %E9
    initvalues(11,1) = 0;   %E10
    initvalues(12,1) = 0;   %E11
    initvalues(13,1) = 0;   %E12
    initvalues(14,1) = 0;   %E13
    initvalues(15,1) = 0;   %E14
    

output_enhancer = zeros(NUM_CELLS,481); %container to store the outputs
r = (1:NUM_CELLS); 
for i = r
    disp(i);
    tf = data_smooth(i,:); %all already cut to 8hrs 
    time = linspace(START_TIME, END_TIME, length(tf));
    [tsim1, results1] = ode15s(@(t,x) chromatinOde(t, x, time,{}, tf),[START_TIME END_TIME],initvalues);
    
%   [t_sim, x_sim] = ode15s(@(t,y)'chromatinOde', v.SIM_TIME,starting_vals,ode_opt,v);

    output = transpose(interp1(tsim1,results1,START_TIME:END_TIME, 'linear'));
    output_enhancer(i,:) = output(15,:); %store all the outputs

end

enhancer{1,1} = output_enhancer;
% 
% mat = enhancer{1,1};
% mat = sortrows(mat, 100);
% surf(mat);
% shading interp 
% colorbar

% output_enhancer = sortrows(output_enhancer,[10 90]);
subplot(1,length(names),j);
imagesc(output_enhancer);
title(char(names(j)));
colorbar
save(strcat('output_enhancer_',data_name,'.mat'), 'output_enhancer')
end
%%
%replot input 
names = {'TNF10ng_762', 'P3CSK4100ng_547','CpG100nM_610', 'LPS100ng_756','PIC50ug_566',};
% names = {'TNF10ng_762','aKO_TNF3.3ng', 'LPS100ng_756','aKO_LPS33ng'};
% names = {'TNF10ng_762', 'ikbamut_10ngTNF'};
for j = 1:length(names)
    data_name = char(names(j));
    data = load(strcat('F://enhancer_dynamics/nfkb_trajectories_08142019/nfkb_dynamics_',data_name,'.mat'));
    data = cell2struct(struct2cell(data), {'trajectories'});
    if istable(data.trajectories)
        data = table2array(data.trajectories);
    else
        data = data.trajectories;
    end
    data_smooth = smoothrows(data);
    data_smooth(any(isnan(data_smooth), 2), :) = []; %remove NaN rows

    datazero = vec2mat(data_smooth(:, 1), 1);
    subtract = repmat(datazero, 1, size(data,2));
    data_smooth = data_smooth - subtract; %subtract the first column to normalize
    data_smooth(data_smooth<0) = 0; %take neg. values to be 0
    maxA = max(data_smooth, [], 2); %absolute row max
    [~, index] = sort(maxA);
    data_smooth = data_smooth(index, :);
    save(strcat('model_input_nfkb_dynamics_',data_name,'.mat'), 'data_smooth')

%     data_smooth = sortrows(data_smooth, {'maxA'}); %max within first 100 mins

    % Plot heatmap, new colormap
%     mod_colormap1 = [123 59 48]/255;
%     mod_colormap2 = [17 103 177]/255;  44,46,67; 54,94,130
%     mod_colormap = [mod_colormap2; [220 221 217]/255; mod_colormap1];
%     mod_colormap = [[transpose(linspace(54,220)) transpose(linspace(94,221)) transpose(linspace(130,217))];
%     [transpose(linspace(220, 123)) transpose(linspace(221, 59)) transpose(linspace(217,48))]]/255;
 
    subplot(1,length(names),j);
    data = data_smooth;
    clims = [0 4];
    imagesc(data, clims);
    loadcolormaps;
    colormap(colormaps.byr);
    
%     colormap(mod_colormap);
    title(char(names(j)));
    shading interp 
    colorbar

end
%%
%replot output
names = {'TNF10ng_762', 'P3CSK4100ng_547','CpG100nM_610', 'LPS100ng_756','PIC50ug_566',};
% names = {'TNF10ng_762', 'ikbamut_10ngTNF'};
for j = 1:length(names)
    data_name = char(names(j));
    data = load(strcat('output_enhancer_',data_name,'.mat'));
    subplot(1,length(names),j);
    data = data.output_enhancer;
    imagesc(data);
    mod_colormap = [[transpose(linspace(0,220)) transpose(linspace(0,221)) transpose(linspace(128,217))];
    [transpose(linspace(220, 178)) transpose(linspace(221, 34)) transpose(linspace(217,34))]]/255;
    colormap(mod_colormap);
    title(char(names(j)));
    colorbar

end
%%
%plot summary curves, column means
names = {'TNF10ng_762', 'P3CSK4100ng_547','CpG100nM_610', 'LPS100ng_756','PIC50ug_566',};
% names = {'TNF10ng_762', 'ikbamut_10ngTNF'};
for j = 1:length(names)
    data_name = char(names(j));
    data = load(strcat('output_enhancer_',data_name,'.mat'));
    subplot(1,length(names),j);
    data = data.output_enhancer;
    % sdev = std(data);
    % plot(smoothrows(mean(data)));
    stdshade(data, 0.2);
    ylim([0 1]);
    xlim([0 480]);
    title(char(names(j)));

end

%%
%plot boxplot of maximums
names = {'TNF10ng_762', 'P3CSK4100ng_547','CpG100nM_610', 'LPS100ng_756','PIC50ug_566',};
for j = 1:length(names)
    data_name = char(names(j));
    data = load(strcat('output_enhancer_',data_name,'.mat'));
    subplot(1,length(names),j);
    data = data.output_enhancer;

    test =max(data,[],2);
    disp(median(test));
    violin(test, 'b', 0.2);
    ylim ([0 1]);
    title(char(names(j)));
    hold on
    x=repmat(1:1,length(test),1);
    scatter(x(:),test(:),10,'filled','jitter','on','jitterAmount',0.15);
end

% load('F://enhancer_dynamics/model_v1/output_enhancer_lps_100ng.mat');
% plot(smoothrows(mean(output_enhancer)));
% ylim([0 100]);