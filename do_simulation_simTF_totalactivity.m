%% ------ simulated TF vectors - osc vs non-osc



%%  SIMULATION
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 8*60;

global START_TIME END_TIME 


START_TIME =0;
END_TIME = 480;
% % END_TIME = 900;



% names = {'nfkb_curves_TNF10ng', 'nfkb_curves_PAM3CSK100ng', 'nfkb_curves_CpG330nM',  'nfkb_curves_LPS100ng','nfkb_curves_pic50ug'};
% names = {'nfkb_curves_TNF10ng', 'nfkb_curves_PAM3CSK100ng', 'nfkb_curves_CpG330nM'};
% names = {'nfkb_curves_TNF10ng','nfkb_curves_TNF_ikbamm'};
% names = {'nfkb_oscillatory','nfkb_nonoscillatory', 'nfkb_oscillatory_hiamp', 'nfkb_nonoscillatory_hiamp'};
% names = {'nfkb_oscillatory','nfkb_nonoscillatory'};
% names = {'nfkb_oscillatory_2xtotalactivity','nfkb_persistent_2xtotalactivity'}; %use END_TIME=900 if this TF sim

output_container = zeros(17, 2);
figure;
n=1;
for j = linspace(.75, 7, 17)
    
    % Transcription factor vector
    pulse_on = j; % (hrs)
    pulse_off = .75; % (hrs)
    tf = zeros(1, END_TIME-START_TIME);
    start_idx = 60*pulse_off;
    tf(round(start_idx):round(start_idx+60*pulse_on)) = 1;
    plot(tf) 
    ylim([0 1.2]);
    hold on;
    
    m=1;
    for k = [0.1, 0.15 0.2, 0.3, 0.4, 0.5, 1.0] %linspace(0.5, 1, 2) %0.5 or 1
        disp(k);
        data_use = (tf)*8*k; %times 8 to get on the same scale as real data (max 0.25nM)

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

        tf = transpose(data_use); %cut to 8hrs
        time = linspace(START_TIME, END_TIME, length(tf));
        [tsim1, results1] = ode15s(@(t,y) chromatinOde(t, y, time,{}, tf),[START_TIME END_TIME],initvalues);

        output = transpose(interp1(tsim1,results1,START_TIME:END_TIME, 'linear'));

        %plot single simulation
%         plot(output(15,:));
%         xlim ([0 480]);
%         ylim([0 1]);
%         hold on;

    max_enhancer = max(output(15,:));
    output_container(n,m) = max_enhancer;
    m = m+1;
    end

    n=n+1;
end

%%
%plot output_container, max chromatin opening
output_container(:,m+1) = start_idx+60*linspace(.75, 7, 17);
figure;
plot(output_container(:,m+1), output_container(:,1));
hold on;
plot(output_container(:,m+1), output_container(:,2));
hold on;
plot(output_container(:,m+1), output_container(:,3));
hold on;
plot(output_container(:,m+1), output_container(:,4));
hold on;
plot(output_container(:,m+1), output_container(:,5));
hold on;
plot(output_container(:,m+1), output_container(:,6));
hold on;
plot(output_container(:,m+1), output_container(:,7));
hold off

legend('non-oscillatory amp0.1','non-oscillatory amp0.15','non-oscillatory amp0.2','non-oscillatory amp0.3','non-oscillatory amp0.4','non-oscillatory amp0.5','non-oscillatory amp1.0')
% legend('non-oscillatory amp0.5','non-oscillatory amp0.6','non-oscillatory amp0.8','non-oscillatory amp1.0','non-oscillatory amp1.2')
% legend('non-oscillatory amp0.5','non-oscillatory amp1')

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