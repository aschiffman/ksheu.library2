
% TNF SIMULATION
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 8*60;

global START_TIME END_TIME 


START_TIME =0;
END_TIME = 480;
% % END_TIME = 900;
% Transcription factor vector


% names = {'nfkb_curves_TNF10ng', 'nfkb_curves_PAM3CSK100ng', 'nfkb_curves_CpG330nM',  'nfkb_curves_LPS100ng','nfkb_curves_pic50ug'};
% names = {'nfkb_curves_TNF10ng', 'nfkb_curves_PAM3CSK100ng', 'nfkb_curves_CpG330nM'};
% names = {'nfkb_curves_TNF10ng','nfkb_curves_TNF_ikbamm'};
% names = {'nfkb_oscillatory','nfkb_nonoscillatory', 'nfkb_oscillatory_hiamp', 'nfkb_nonoscillatory_hiamp'};
names = {'nfkb_oscillatory','nfkb_nonoscillatory'};
% names = {'nfkb_oscillatory_2xtotalactivity','nfkb_persistent_2xtotalactivity'}; %use END_TIME=900 if this TF sim


output_container = zeros(12, 2);
for j = 1:length(names)
    data_name = char(names(j));
    data = load(strcat(data_name,'.mat'));
    data = cell2struct(struct2cell(data), {'nfkb_curves'});
    % data = (data.nfkb_curves)*30;
    
    data_use = (data.nfkb_curves)/4; %times 8 to get on the same scale as real data

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

    %% scan different parameters
    n=1;
    for forward_factor = [1.0 1.1 1.2 1.3 1.4 1.5] 
        disp(forward_factor);
    for reverse_factor = [0.5 0.6 0.7 0.8 0.9 1.0]
    
    k1 =10;
    ratio = 7.5;
    
    p_mod = [
    1 1 k1 % k1
    3 1 k1*forward_factor % k1
    5 1 k1*forward_factor^2 % k1
    7 1 k1*forward_factor^3 % k1
    9 1 k1*forward_factor^4% k1
    11 1 k1*forward_factor^5 % k1
    13 1 k1*forward_factor^6% k1
    15 1 k1*forward_factor^7 % k1
    17 1 k1*forward_factor^8 % k1
    19 1 k1*forward_factor^9 % k1
    21 1 k1*forward_factor^10 % k1
    23 1 k1*forward_factor^11 % k1
    25 1 k1*forward_factor^12 % k1
    27 1 k1*forward_factor^13 % k1
    
    2 1 k1*ratio % k1*ratio
    4 1 k1*ratio*reverse_factor % k1*ratio
    6 1 k1*ratio*reverse_factor^2 % k1*ratio
    8 1 k1*ratio*reverse_factor^3 % k1*ratio
    10 1 k1*ratio*reverse_factor^4% k1*ratio
    12 1 k1*ratio*reverse_factor^5 % k1*ratio
    14 1 k1*ratio*reverse_factor^6% k1*ratio
    16 1 k1*ratio*reverse_factor^7 % k1*ratio
    18 1 k1*ratio*reverse_factor^8 % k1*ratio
    20 1 k1*ratio*reverse_factor^9 % k1*ratio
    22 1 k1*ratio*reverse_factor^10 % k1*ratio
    24 1 k1*ratio*reverse_factor^11 % k1*ratio
    26 1 k1*ratio*reverse_factor^12 % k1*ratio
    28 1 k1*ratio*reverse_factor^13 % k1*ratio
    ];
%%
%     for Hill = [0.01 0.05 0.1 0.5 1 1.5 2 2.5 3 4 5 6]
%         disp(Hill);
%     for Kd = [0.01 0.1 0.5 0.8 1 1.2 2 5 10]/32 
%         disp(Kd);
%     p_mod = [
%     1 3 Kd % Kd1
%     3 3 Kd % Kd1
%     5 3 Kd % Kd1
%     7 3 Kd % Kd1
%     9 3 Kd % Kd1
%     11 3 Kd % Kd1
%     13 3 Kd % Kd1
%     15 3 Kd % Kd1
%     17 3 Kd % Kd1
%     19 3 Kd % Kd1
%     21 3 Kd % Kd1
%     23 3 Kd % Kd1
%     25 3 Kd % Kd1
%     27 3 Kd % Kd1
%     
%     1 2 Hill % Hill
%     3 2 Hill % Hill
%     5 2 Hill % Hill
%     7 2 Hill % Hill
%     9 2 Hill % Hill
%     11 2 Hill % Hill
%     13 2 Hill % Hill
%     15 2 Hill % Hill
%     17 2 Hill % Hill
%     19 2 Hill % Hill
%     21 2 Hill % Hill
%     23 2 Hill % Hill
%     25 2 Hill % Hill
%     27 2 Hill % Hill
%     ];


    tf = transpose(data_use); %cut to 8hrs
    time = linspace(START_TIME, END_TIME, length(tf));
    [tsim1, results1] = ode15s(@(t,y) chromatinOde_pmod(t, y, time,{}, tf, p_mod),[START_TIME END_TIME],initvalues);

    output = transpose(interp1(tsim1,results1,START_TIME:END_TIME, 'linear'));

    %plot single simulation
    if strcmp(data_name, 'nfkb_oscillatory')==1
        plot(output(15,:), 'color',[0 0 1]); %blue
    else
        plot(output(15,:), 'color',[1 0.5 0]); %orange
    end
    xlim ([0 480]);
    ylim([0 1]);
    hold on;

    max_enhancer = max(output(15,:));
    output_container(n,j) = max_enhancer;
    n=n+1;
    end
    end
end
%%
%plot output_container for  reverse and forwrad factor sweep
output_container(:,3) = output_container(:,2)./output_container(:,1);
i=1;
for forward_factor = [1.0 1.1 1.2 1.3 1.4 1.5]
    disp(forward_factor);
    for reverse_factor = [0.5 0.6 0.7 0.8 0.9 1.0]
        disp(reverse_factor);
        output_container(i,4) = forward_factor;
        output_container(i,5) = reverse_factor;
        i= i+1;
    end
end

Data = xyz2grid(output_container(:,4), output_container(:,5), output_container(:,3));
figure;
imagesc(Data);
colorbar;

%%
%plot output_container for Hill and Kd sweep
output_container(:,3) = output_container(:,2)./output_container(:,1);
i=1;
for Hill = [0.01 0.05 0.1 0.5 1 1.5 2 2.5 3 4 5 6]
        disp(Hill);
    for Kd = [0.01 0.1 0.5 0.8 1 1.2 2 5 10]/32 
        disp(Kd);
        output_container(i,4) = Hill;
        output_container(i,5) = Kd;
        i= i+1;
    end
end

Data = xyz2grid(output_container(:,4), output_container(:,5), output_container(:,1));
figure;
imagesc(Data);
colorbar;

%%
%plot fold change in max chromatin opening for osc. vs non-osc
output_container(:,3) = [0.1 1 2 5 10 20 30 40 50 60];
figure;
plot(output_container(:,3), (output_container(:,2)./output_container(:,1)));

xlabel('model parameter');
ylabel('fold-change (non-osc/osc)');
%%
% paired boxplot of max chromatin opening
N = length(output_container);
x = output_container(:,1);
y = output_container(:,2);
all = [x;y];
g=[1*ones(N,1);2*ones(N,1)];

figure;
boxplot(all, g);
hold on
scatter(g(:),all(:),30,'filled','jitter','on','jitterAmount',0.15);
% Plot lines between corresponding pairs
for k = 1 : length(x)
  plot([x(k), y(k)], ...
    'rs-', 'LineWidth', .5, 'MarkerSize', 5);
  hold on
end
ylabel('Max chromatin opening');

% line(repmat([(1:2).';NaN], [N,1]), ...
%   reshape(measures(1:N,[1:2, 1]).', [], 1), ...
%   'Color', 0.7*[1 1 1], 'Marker', '.', 'MarkerSize', 10);

%%
%plot simTFs

names = {'nfkb_oscillatory','nfkb_nonoscillatory'};
figure;
for j = 1:length(names)
    data_name = char(names(j));
    data = load(strcat(data_name,'.mat'));
    data = cell2struct(struct2cell(data), {'nfkb_curves'});
    data = (data.nfkb_curves)*1;
    
    
    plot(data);
    xlim ([0 480]);
    ylim ([0 1.2]);
    hold on;
    
end