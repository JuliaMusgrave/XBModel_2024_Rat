
% This script is set up to facilitate fitting of various permutations of
% the linearised model to experimental data for both Pi and ATP as done in
% the paper "Analysis of metabolite and strain effects on cardiac 
% cross-bridge dynamics using model linearisation techniques

% Author: Julia Musgrave
% Date: Oct 2023

clear

load('final_fit.mat','xs')
% cell array containing an initial guess (currently using previous best
% parameters, in 16x4 array)
init_xs=xs;


%% Loading data for fitting

% load experimental data for fitting

% experimental frequencies stored in 'freqs' array
% remaining data stored in a 7xn cell array of mixed data types, called 'data'. 
% col 1 contains row descriptions and cols 2-n contain the data for n-1 
% different metabolite conditions.
% Rows:
% 1: 'Solution' - char arrays with description of the metab condition
% 2: 'concentrations' - 1x2 array of ATP and Pi concentrations (mM)
% 3: 'mean CM' - mean complex modulus data, a complex double array matching
% the size of 'freqs' (MPa)
% 4: 'EM SE' - SEM of the real component of the complex moduli, array
% matching the size of 'freqs' (MPa)
% 5: 'VM SE' - SEM of the imaginary component of the complex moduli, array
% matching the size of 'freqs' (MPa)
% 6: 'mean F0' - mean steady-state active stress (kPa)
% 7: 'F0 SE' - SEM of the steady-state active stress (kPa)
load('rat_data.mat','data','freqs')
no_conds=width(data)-1;

%% Model fitting

% general set up
no_ps = 14; % max no of parameters
lower_bounds  =   [0.5  0.5  0.5  0.1   0.5   0.5  0.005 1  2e3    0    1    1   0.1  0.05];
upper_bounds =    [200  200  200  200   200   10    0.5  5  1e5    45  1000 1000  5    5  ];
idx_keys={'1' '{-1}' '2' '{-2}' '3'};
curent_key='None';
col_keys = {'(A) Pi:M1,ATP:M1', '(B) Pi:M2,ATP:M1','(C) Pi:M1,ATP:M2', '(D) Pi:M2,ATP:M2'};

% fitting to different model perms (strain on each pair of rates x 4 metabolite
% dependencies)
met_assign=1:4; % metabolite dependencies to consider (up to 4 options)
% 1=(A)  2=(B)  3=(C)  4=(D) - see col_keys above

strain_loc=1:16; % strain dependencies to consider (up to 16 options)

% begin iterating through the different permutations
for s=strain_loc  % going through each of the perms of strains

    sd=find_sd_array(s);

    for md=met_assign  % going through the metabolite permutations
    
    % fit active CM for all concentrations to linear model
    guess = init_xs{s,md};
    ub=upper_bounds;
    lb=lower_bounds;

    % removing range for parameters not being fitted
    if s<=6;  ub(12)=0; lb(12)=0; guess(12)=0; end  % only one strain dep
    if s==1;  ub(11)=0; lb(11)=0; guess(11)=0; end   % no strain dep at all
    if mod(md,2);  ub(14)=0; lb(14)=0; guess(14)=0; end % M1 for Pi
    if md<=2; ub(13)=0; lb(13)=0; guess(13)=0; end % M1 for Pi

    true_no_params(s,md) = sum(ub>0);
    
    % fitting function (see bottom section)
    fun=@(x)XBmodel_APfit(x, sd, md, freqs, data);
    
    options=optimoptions('particleswarm','HybridFcn','patternsearch',...
        'InitialSwarmMatrix',guess,'Display','none','SwarmSize',1000); 
    % options without a initial parameter guess
    %options=optimoptions('particleswarm','HybridFcn','patternsearch','Display','none','SwarmSize',1000);
    
    [x,OBJ]=particleswarm(fun,no_ps,lb,ub,options);
    
    % storing optimal objective function and parameters in arrays
    OBJs(s,md)=OBJ; 
    xs{s,md}=x;
    end
    
    % just for correct table labelling purposes
    if s>1
    if s<=6
        curent_key=['k_' idx_keys{s-1}];
    else
        i=find(sd);
        curent_key=['k_' idx_keys{i(1)} ' and k_' idx_keys{i(2)}];
    end
    end
    row_keys(s)={curent_key};  

end
% Normalised RMSE (as % of range) table (Table 3 in manuscript)
Table3=array2table(round(OBJs/0.556*100,2),'RowNames',row_keys,'VariableNames',col_keys);

% Data for Figure 5, % improvement from simplest model
RMSE_improve=array2table((OBJs(1,1)-OBJs)/OBJs(1,1)*100,'RowNames',row_keys,'VariableNames',col_keys);


%% plotting best fit found
blue=[0 0.447 0.741];
red=[0.85 0.325 0.098];
purple=[0.494 0.184 0.556];
green=[0.466 0.674 0.188];

[minO,minI]=min(OBJs,[],'all');
min_row=mod(minI,16)+16*(mod(minI,16)==0);
md=ceil(minI/16);
xmin=xs{min_row,md};
sd=find_sd_array(min_row);

% plotting details
EM_pos=[0.075 0.15 0.4 0.8];
VM_pos=[0.575 0.15 0.4 0.8];
colour={[0 0 0],red,blue,green,purple};
figure('Position',[687,365,813,335])

for i=1:no_conds % plotting each metabolite
    [~,Yam]=XBmodel_2024_linear_perms(xmin,sd,md,data{2,i+1});

    subplot('Position',EM_pos)
    errorbar(freqs,real(data{3,i+1}),data{4,i+1},'.','LineStyle','none','MarkerSize',12,'Color',colour{i})
    hold on
    ax=gca;
    ax.XScale='log';
    semilogx(logspace(-1,2,100),real(Yam),'LineWidth',1,'Color',colour{i})

    subplot('Position',VM_pos)
    errorbar(freqs,imag(data{3,i+1}),(data{5,i+1}),'.','LineStyle','none','MarkerSize',12,'Color',colour{i})
    hold on
    ax=gca;
    ax.XScale='log';
    semilogx(logspace(-1,2,100),imag(Yam),'LineWidth',1,'Color',colour{i})

end

% tidying figure
subplot('Position',EM_pos)
legend('','Baseline','','0.1 mM ATP','','1 mM ATP','','0 mM Pi','','5 mM Pi','Location','northwest')
xlim([0.1 99])
xticks([0.1 1 10 99])
xticklabels({'0.1' '1' '10' '100'})
xlabel('Frequency (Hz)')
ylabel('Elastic Modulus (MPa)')
subplot('Position',VM_pos)
xlim([0.1 99])
xticks([0.1 1 10 99])
xticklabels({'0.1' '1' '10' '100'})
xlabel('Frequency (Hz)')
ylabel('Viscous Modulus (MPa)')

%% Functions

% Multiple condition fitting function
function OBJ=XBmodel_APfit(x,sd,md,freqs,data)
% Inputs:   x: vector of parameters to optimise
%           sd: strain dependence array (see find_sd_array)
%           md: code for metabolite dependence   
%           freqs: frequencies the CM was collected at
%           data: cell array with data stored in format specified in above
%           script
% Output:   OBJ: objective value (averaged across all conditions)

    OBJ=0; 
    for i=1:width(data)-1
        Y=data{3,i+1};
        F0=data{6,i+1};
        Fe=data{7,i+1};
        mets=data{2,i+1};
        OBJi=XBmodel_2024_linear_perms(x,sd,md,mets,Y,freqs,F0,Fe);
        OBJ=OBJ+OBJi;
    end
    OBJ=OBJ/i;
end

% Converting strain dependence number (table row) to logical array
function sd=find_sd_array(s)
% returns a 1x5 array where 1s represent locations with strain dependence
% [k1 k-1 k2 k-2 k3]
    
    sd=zeros(1,5);
    if s<=1
        return
    elseif s<=6
        sd(s-1)=1;
    elseif s<=10
        sd([1 s-5])=1;
    elseif s<=13
        sd([2 s-8])=1;
    elseif s<=15
        sd([3 s-10])=1;
    else
        sd([4 5])=1;
    end

end