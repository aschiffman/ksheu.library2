
% TNF SIMULATION
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 650;

global START_TIME END_TIME 


START_TIME = 0;
END_TIME = 650;
% % END_TIME = 900;
% Transcription factor vector


% names = {'nfkb_curves_TNF10ng', 'nfkb_curves_PAM3CSK100ng', 'nfkb_curves_CpG330nM',  'nfkb_curves_LPS100ng','nfkb_curves_pic50ug'};
% names = {'nfkb_curves_TNF10ng', 'nfkb_curves_PAM3CSK100ng', 'nfkb_curves_CpG330nM'};
% names = {'nfkb_curves_TNF10ng','nfkb_curves_TNF_ikbamm'};
% names = {'nfkb_oscillatory','nfkb_nonoscillatory', 'nfkb_oscillatory_hiamp', 'nfkb_nonoscillatory_hiamp'};
% names = {'nfkb_oscillatory','nfkb_nonoscillatory'};
names = {'nfkb_oscillatory_ext','nfkb_nonoscillatory_ext'};
% names = {'nfkb_oscillatory_2xtotalactivity','nfkb_persistent_2xtotalactivity'}; %use END_TIME=900 if this TF sim

output_container = zeros(21, 2);
num_amplitudes = 21*length(names);
output_unmod = zeros(num_amplitudes, 651);
output_mod = zeros(num_amplitudes, 651);
for j = 1:length(names)
    data_name = char(names(j));
    data = load(strcat(data_name,'.mat'));
    data = cell2struct(struct2cell(data), {'nfkb_curves'});
    

    
    n=1;
    for k = linspace(0, 1, 21) % cycle through amplitudes
%         k = 1;
        disp(k);
        data_use = (data.nfkb_curves)*8*k; %times 8 to get on the same scale as real data (max 0.25nM)

        
        % Starting Conditions
        initvalues = zeros(30,1);
        initvalues(1,1) = 1;    %E0 - closed
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
        initvalues(15,1) = 0;   %E14 - open
        initvalues(16,1) = 0;   %Eh0 - closed, modified
        initvalues(17,1) = 0;   %Eh1
        initvalues(18,1) = 0;   %Eh2
        initvalues(19,1) = 0;   %Eh3
        initvalues(20,1) = 0;   %Eh4
        initvalues(21,1) = 0;   %Eh5
        initvalues(22,1) = 0;   %Eh6
        initvalues(23,1) = 0;   %Eh7
        initvalues(24,1) = 0;   %Eh8
        initvalues(25,1) = 0;   %Eh9
        initvalues(26,1) = 0;   %Eh10
        initvalues(27,1) = 0;   %Eh11
        initvalues(28,1) = 0;   %Eh12
        initvalues(29,1) = 0;   %Eh13
        initvalues(30,1) = 0;   %Eh14 - open, modified

        tf = transpose(data_use); %cut to 8hrs
        time = linspace(START_TIME, END_TIME, length(tf));
        [tsim1, results1] = ode15s(@(t,y) chromatinOde(t, y, time,{}, tf),[START_TIME END_TIME],initvalues);

        output = transpose(interp1(tsim1,results1,START_TIME:END_TIME, 'linear'));
        if j==1
            output2 = output;
        end       

        %store single simulation
        output_unmod(j*n,:) = output(15,:);
        output_mod(j*n,:) = output(30,:);

        %store all values, amplitude = 1
        %TODO: generalize this for more than two names
        if k==1
            if j==1
                output_osc = output; % each row is a time simulation of a state 1-30, see SPECIES in chromatinInitialize.m
            elseif j==2
                output_nonosc = output;
            end
        end
        
    max_enhancer = max(output(15,:));
    max_enhancer_mod = max(output(30,:));
    output_container(n,j) = max_enhancer;
    output_container_mod(n,j) = max_enhancer_mod;
    n=n+1;
    end

end

%plot single simulation, mmodified open fraction over time
figure;
for k = 1:size(output_mod,1)
    plot(output_mod(k,:)); %15=E_0 open, 30=E_h0 open
    xlim ([0 650]);
    ylim([0 1]);
    hold on;
end
hold off;
xlabel('time (min)');
ylabel('fraction (E_{h0})');

%plot single simulation, unmmodified open fraction over time
figure;
for k = 1:size(output_unmod,1)
    plot(output_unmod(k,:)); %15=E_0 open, 30=E_h0 open
    xlim ([0 650]);
    ylim([0 1]);
    hold on;
end
hold off;
xlabel('time (min)');
ylabel('fraction (E_{0})');


%%

% % Plot single simulation, mult states, amplitude = 1
% figure;
% hold on;
% plot(output_osc(1,:)); %open, unmod, non osc
% plot(output_osc(16,:)); %closed, mod, non osc
% plot(output_nonosc(1,:)); %closed, unmod, non osc
% plot(output_nonosc(16,:)); %closed, mod, non osc
% hold off;
% xlim ([0 480]);
% ylim([0 1]);
% hold on;
% 
% hold off;
% xlabel('time (min)');
% ylabel('fraction at amp=1');
% legend ("E_{14} non-osc","E_{h14} non-osc", "E_{14} osc", "E_{h14} osc");
% 
% % Plot single simulation, mult states, amplitude = 1
% figure;
% hold on;
% 
% plot(output_osc(15,:)); %open, unmod, non osc
% plot(output_osc(30,:)); %open, mod, non osc
% plot(output_nonosc(15,:)); %open, unmod, non osc
% plot(output_nonosc(30,:)); %open, mod, non osc
% hold off;
% xlim ([0 480]);
% ylim([0 1]);
% hold on;
% 
% hold off;
% xlabel('time (min)');
% ylabel('fraction at amp=1');
% legend ("E_0 non-osc","E_{h0} non-osc", "E_0 osc","E_{h0} osc");

% % Area single simulation, mult states, amplitude = 1
% figure;
% area(transpose(output_osc));
% colormap winter;
% xlim ([0 480]);
% ylim([0 1]);
% xlabel('time (min)');
% ylabel('(osc) fraction at amp=1');
% 
% % Area single simulation, mult states, amplitude = 1
% figure;
% area(transpose(output_nonosc));
% xlim ([0 480]);
% ylim([0 1]);
% xlabel('time (min)');
% ylabel('(non-osc) fraction at amp=1');

% 
% %plot output_container, max chromatin opening
% output_container(:,3) = linspace(0,1,21);
% figure;
% plot(output_container(:,3), output_container(:,1));
% hold on;
% plot(output_container(:,3), output_container(:,2));
% hold off
% 
% xlabel('amplitude (fold)');
% ylabel('fraction (E_0)');
% legend('oscillatory','non-oscillatory')
% 
% % plot max methylated chromatin opening
% output_container_mod(:,3) = linspace(0,1,21);
% figure;
% plot(output_container_mod(:,3), output_container_mod(:,1));
% hold on;
% plot(output_container_mod(:,3), output_container_mod(:,2));
% hold off
% 
% xlabel('amplitude (fold)');
% ylabel('fraction (E_{h0})');
% legend('oscillatory','non-oscillatory')
% 
% % plot max total chromatin opening
% output_container(:,3) = linspace(0,1,21);
% output_container_total = output_container;
% output_container_add = [output_container_mod(:,1:2),(zeros(size(output_container_mod,1),1))];
% output_container_total = output_container_total + output_container_add;
% 
% figure;
% plot(output_container_total(:,3), output_container_total(:,1));
% hold on;
% plot(output_container_total(:,3), output_container_total(:,2));
% hold off
% 
% xlabel('amplitude (fold)');
% ylabel('fraction (E_{h0} + E_0)');
% legend('oscillatory','non-oscillatory')

%plot max methylated chromatin opening vs reaction rate constant
% TODO: maybe continue to read in reaction rates as currently structured,
% but iterate over rates near reaction rate to get plot. Have to think
% about this a little more.

% plot Eh_0 fraction at various rxn rates vs. time, with amplitude=1
% TODO

% %%
% %plot fold change in max chromatin opening for osc. vs non-osc
% output_container(:,3) = linspace(0,1,21);
% figure;
% plot(output_container(:,3), output_container(:,2)/output_container(:,1));
% 
% xlabel('amplitude (fold)');
% ylabel('fold-change (non-osc/osc)');

%plot simTFs

figure;
for j = 1:length(names)
    data_name = char(names(j));
    data = load(strcat(data_name,'.mat'));
    data = cell2struct(struct2cell(data), {'nfkb_curves'});
    data = (data.nfkb_curves)*1;
    
    
    plot(data);
    xlim ([0 650]);
    ylim ([0 1.2]);
    hold on;
    
end
