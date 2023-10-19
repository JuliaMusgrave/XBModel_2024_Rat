%

% This script can be run to plot each the figures in "Analysis of metabolite 
% and strain effects on cardiac cross-bridge dynamics using model 
% linearisation techniques"
%
% Contains examples of calling XBmodel_2024_linear_perms and
% XBmodel_2024_Rat (the ode version of the final model)
%
% uses rat_data.mat and final_fit.mat as well as
% XBmodel_2024_linear_perms.m, XBmodel_2024_Rat.m, XBmodel_2022_linear.m
% and Rice_style_Fredev.m

% Author: Julia Musgrave
% Date: October 2023

% reset (can run this section again once you're done)
set_all_interpreters('none') % will set interpreters back to normal
clear
close all

%% Set up

% run this set-up section first to ensure all plots work
set_all_interpreters('latex')
% https://au.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
blue=[0 0.447 0.741];
red=[0.85 0.325 0.098];
purple=[0.494 0.184 0.556];
green=[0.466 0.674 0.188];

%% Figure 1 (General linearisation breakdown)

% Using Musgrave 2022 model and parameters
model=@XBmodel_2022_linear;
x= [1.82, 16.78, 6.57, 99.99, 0.5, 99582, 4.18, 0.13, 1.13, 4.59];
[~,Yam,HxB_comp,HxC_comp,HC_comp,HC_Lcomp,HC_Scomp]=model(x);

% plotting details
EM_pos=[0.075 0.15 0.4 0.8];
VM_pos=[0.575 0.15 0.4 0.8];

figure('Name','General Linearisation (Fig 1)','Position',[687,365,813,335])

two_panel_CMplot(logspace(-1,2,100),Yam,'-',2, [0 0 0])
two_panel_CMplot(logspace(-1,2,100),HxB_comp,'-',2, blue)
two_panel_CMplot(logspace(-1,2,100),HxC_comp,'-',2, red)
two_panel_CMplot(logspace(-1,2,100),HC_comp,'-',2, green)
two_panel_CMplot(logspace(-1,2,100),HC_Lcomp,'-.',1.5, green)
two_panel_CMplot(logspace(-1,2,100),HC_Scomp,'--',1.5, green)

subplot('Position',EM_pos)
text(0.04,max(ylim),'A','FontSize',16,'FontWeight','bold')

l=legend('Combined response','$H_{\mathrm{xB}}$','$H_{\mathrm{xC}}$','$H_\mathrm{C}$','length','strain',...
    'Location','northwest','AutoUpdate','Off','Box','Off');
l.ItemTokenSize=[20,18];

subplot('Position',VM_pos)
text(0.04,max(ylim),'B','FontSize',16,'FontWeight','bold')

two_panel_CMplot(logspace(-1,2,100),Yam,'-',2, [0 0 0])

%% Figure 3 (CM rat data with errors)

% load data
load('rat_data.mat')

% defining the indices where means and se's are stored
i_m=3;
i_re=4;
i_ie=5;
i_x_top=[2 4 3]; % ATP variation
i_x_bottom=[2 5 6];   % Pi variation

% plotting details
colour={[0 0 0],red,blue};
labels_top={'5 mM ATP','1 mM ATP','0.1 mM ATP'};
labels_bottom={['1 mM ' P_i],['0 mM ' P_i], ['5 mM ' P_i]};
[f,A,B,C,D]=four_panel_plot('ATP/Pi data (Fig 3)');

for j=1:3
    x=i_x_top(j);

    set(f,'CurrentAxes',A)
    errorbar(freqs,real(data{i_m,x}),data{i_re,x},'.','LineStyle','none','MarkerSize',14,'Color',colour{j})

    set(f,'CurrentAxes',B)
    errorbar(freqs,imag(data{i_m,x}),(data{i_ie,x}),'.','LineStyle','none','MarkerSize',14,'Color',colour{j})

    x=i_x_bottom(j);

    set(f,'CurrentAxes',C)
    errorbar(freqs,real(data{i_m,x}),data{i_re,x},'.','LineStyle','none','MarkerSize',14,'Color',colour{j})

    set(f,'CurrentAxes',D)
    errorbar(freqs,imag(data{i_m,x}),(data{i_ie,x}),'.','LineStyle','none','MarkerSize',14,'Color',colour{j})

end

% plot admin
axs={'A','B','C','D'};
for i=1:4
    ax=eval(axs{i});
    set(f,'CurrentAxes',ax)
    xlim([0.1 99])
    xticks([0.1 1 10 99])
    xticklabels({'0.1' '1' '10' '100'})
    text(0.035,max(ylim),axs{i},'FontSize',16,'FontWeight','bold')
end
set(f,'CurrentAxes',A)
legend(labels_top,'Location','northwest','EdgeColor','none')
set(f,'CurrentAxes',C)
legend(labels_bottom,'Location','northwest','EdgeColor','none')

% significance markers
annotation("line",[0.2673 0.2673],[0.9467 0.8738])
annotation("line",[0.2573 0.2673],[0.8738 0.8738])
annotation("line",[0.2573 0.2673],[0.9467 0.9467])
annotation("line",[0.2573 0.2673],[0.9103 0.9103])
annotation("textbox",[0.264 0.8892 0.0309 0.0414],'String','*','EdgeColor','none','FontSize',15)

annotation("line",[0.2336 0.2336],[0.432 0.3964])
annotation("line",[0.2236 0.2336],[0.432 0.432])
annotation("line",[0.2236 0.2336],[0.3964 0.3964])
annotation("textbox",[0.2295 0.3941 0.0309 0.0414],'String','*','EdgeColor','none','FontSize',15)

%% Figure 4 (generalised effect of ATP and Pi)

model=@XBmodel_2022_linear;
x= [1.82, 16.78, 6.57, 99.99, 0.5, 99582, 4.18, 0.13, 1.13, 4.59];

[~,Yb]=model(x);
maxY=max(real(Yb));
fs=logspace(-1,2,100);

n=10;
ATPs=logspace(-0.05,1.3,n);
Pis=logspace(-0.7,0.85,n);

[f,A,B,C,D]=four_panel_plot('ATP/Pi effect (Fig 4)');

%plotting mass of gradient colours
for i=1:n

    xA=x; xA(3)=x(3)*ATPs(i)/5; % accounting for 5 mM [ATP] in baseline
    [~,YA]=model(xA);

    xP=x; xP(4)=x(4)*Pis(i);
    [~,YP]=model(xP);

    Ci=[0.71-0.71*i/n 0.83-0.66*i/n 0.91-0.64*i/n];

    set(f,'CurrentAxes',A)
    semilogx(fs,real(YA)/maxY,'-','Color',Ci)

    set(f,'CurrentAxes',B)
    semilogx(fs,imag(YA)/maxY,'-','Color',Ci)

    set(f,'CurrentAxes',C)
    semilogx(fs,real(YP)/maxY,'-','Color',Ci)

    set(f,'CurrentAxes',D)
    semilogx(fs,imag(YP)/maxY,'-','Color',Ci)

end

%plotting baseline + annotating
set(f,'CurrentAxes',A)
semilogx(fs,real(Yb)/maxY,'k-','LineWidth',2)
ylabel('Elastic Modulus')
ylim([-0.4 3.1])
annotation('textarrow', [0.33,0.2], [0.62,0.82],'String','decreasing [ATP]')
set(A,'FontSize',11)
text(0.035,max(ylim),'A','FontSize',16,'FontWeight','bold')

set(f,'CurrentAxes',B)
semilogx(fs,imag(Yb)/maxY,'k-','LineWidth',2)
ylabel('Viscous Modulus')
annotation('textarrow', [0.84,0.68], [0.61,0.93],'String','decreasing [ATP]')
set(B,'FontSize',11)
text(0.035,max(ylim),'B','FontSize',16,'FontWeight','bold')

set(f,'CurrentAxes',C)
semilogx(fs,real(Yb),'k-','LineWidth',2)
ylabel('Elastic Modulus')
ylim([-0.3 3.5])
annotation('textarrow', [0.41,0.34], [0.14,0.35],'String',['decreasing [' P_i ']'])
set(C,'FontSize',11)
text(0.035,max(ylim),'C','FontSize',16,'FontWeight','bold')

set(f,'CurrentAxes',D)
semilogx(fs,imag(Yb),'k-','LineWidth',2)
ylabel('Viscous Modulus')
ylim([-0.2 1.6])
annotation('textarrow', [0.92,0.88], [0.14,0.43],'String',['decreasing [' P_i ']'])
set(D,'FontSize',11)
text(0.035,max(ylim),'D','FontSize',16,'FontWeight','bold') % hard-coded


%% Figure 5: bar graph of best improvements

% load improvement in RMSE values table
load('final_fit.mat','RMSE_improve')

% sorting RMSE improvement %s from best to worst
sorted=RMSE_improve;
sorted{:,5}=max(table2array(RMSE_improve),[],2);
sorted=sortrows(sorted,"Var5");

figure('Units', 'normalized' ,'OuterPosition', [0.1,0.3,0.3,0.53],'Name', 'RMSE improvements (Fig 5)')
purp_bar=purple*1.2;
colororder([blue; red; green; purp_bar]);

barh(table2array(sorted(:,1:4)),0.6,'grouped','EdgeColor','flat')

yticklabels(add_dollars(sorted.Properties.RowNames))
yticks(1:height(sorted));
ylabel('Rate with strain dependence')
xlabel('Improvement in RMSE (\%)')
set(gca,'YDir','reverse')
xlim([-2 50])
annotation('arrow',[0.309 0.2],[0.76 0.76])
box off
l=legend({[P_i ': M1, ATP: M1'],[P_i ': M2, ATP: M1'],[P_i ': M1, ATP: M2'],...
    [P_i ': M2, ATP: M2']});
title(l,'Metabolite binding method')


%% Figure 6 (best vs simplest fit)

% load final parameters
load('final_fit.mat','xs')

model=@XBmodel_2024_linear_perms;

sd = [ 0 0 1 1 0];
md=4;
x_b=xs{12,md};
x_simple=xs{1,1};

% load relevant data
load('rat_data.mat')
fs=logspace(-1,2,100);

[f,A,B,C,D]=four_panel_plot('Best vs simplest fit (Fig 6)');

% plotting details
colour={0,[0 0 0],blue,red,blue,red};

for i=[3 4 2]

    set(f,'CurrentAxes',A)
    [~,Yam_b]=model(x_b,sd,md,data{2,i});
    [~,Yam_simple]=model(x_simple,zeros(1,5),1,data{2,i});
    
    errorbar(freqs,real(data{3,i}),data{4,i},'.','LineStyle','none','MarkerSize',12,'Color',colour{i})
    semilogx(fs,real(Yam_b),'-' ,'Color',colour{i})
    semilogx(fs,real(Yam_simple),'--' ,'Color',colour{i})

    set(f,'CurrentAxes',B)
    errorbar(freqs,imag(data{3,i}),(data{5,i}),'.','LineStyle','none','MarkerSize',12,'Color',colour{i})
    semilogx(fs,imag(Yam_b),'-' ,'Color',colour{i})
    semilogx(fs,imag(Yam_simple),'--' ,'Color',colour{i})

end

for i=[5 2 6]

    set(f,'CurrentAxes',C)
    [~,Yam]=model(x_b,sd,md,data{2,i});
    [~,Yam_simple]=model(xs{1,1},zeros(1,5),1,data{2,i});
    
    errorbar(freqs,real(data{3,i}),data{4,i},'.','LineStyle','none','MarkerSize',12,'Color',colour{i})
    semilogx(logspace(-1,2,100),real(Yam),'-' ,'Color',colour{i})
    semilogx(logspace(-1,2,100),real(Yam_simple),'--' ,'Color',colour{i})

    set(f,'CurrentAxes',D)
    errorbar(freqs,imag(data{3,i}),(data{5,i}),'.','LineStyle','none','MarkerSize',12,'Color',colour{i})
    semilogx(logspace(-1,2,100),imag(Yam),'-' ,'Color',colour{i})
    semilogx(logspace(-1,2,100),imag(Yam_simple),'--' ,'Color',colour{i})

end

set(f,'CurrentAxes',A)
text(0.035,max(ylim),'A','FontSize',16,'FontWeight','bold')
legend('','0.1 mM ATP','','','1 mM ATP','','','5 mM ATP','','Location','northwest')

set(f,'CurrentAxes',B)
text(0.035,max(ylim),'B','FontSize',16,'FontWeight','bold')

set(f,'CurrentAxes',C)
text(0.035,max(ylim),'C','FontSize',16,'FontWeight','bold')
legend('',['0 mM ' P_i],'','',['1 mM ' P_i],'','',['5 mM ' P_i],'','Location','northwest')

set(f,'CurrentAxes',D)
text(0.035,max(ylim),'D','FontSize',16,'FontWeight','bold')

%% Figure 7 (final fit with errors)

% load final parameters
load('final_fit.mat','final_params')

sd = [ 0 0 0 1 1];
md=4;
x=final_params;

% load relevant data
load('rat_data.mat')

% plotting details
EM_pos=[0.075 0.15 0.4 0.8];
VM_pos=[0.575 0.15 0.4 0.8];
colour={[0 0 0],red,blue,green,purple};

include_error=1;

figure('Name','Final fit fig (Fig 7)')
for i=1:5
    [~,Yam]=XBmodel_2024_linear_perms(x,sd,md,data{2,i+1});
    if ~include_error
    two_panel_CMplot(freqs,data{3,i+1},'.',2, colour{i})
    else
    subplot('Position',EM_pos)
    errorbar(freqs,real(data{3,i+1}),data{4,i+1},'.','LineStyle','none','MarkerSize',12,'Color',colour{i})
    hold on
    ax=gca;
    ax.XScale='log';

    subplot('Position',VM_pos)
    errorbar(freqs,imag(data{3,i+1}),(data{5,i+1}),'.','LineStyle','none','MarkerSize',12,'Color',colour{i})
    hold on
    ax=gca;
    ax.XScale='log';

    end
    two_panel_CMplot(logspace(-1,2,100),Yam,'-' ,1,colour{i})
end
set(gcf,'Position',[687,365,813,335])

subplot('Position',EM_pos)
text(0.04,max(ylim),'A','FontSize',16,'FontWeight','bold')
l=legend('','Baseline','','0.1 mM ATP','','1 mM ATP','',['0 mM ' P_i],'',['5 mM ' P_i],'Location','northwest');
l.ItemTokenSize=[20,18];

subplot('Position',VM_pos)
text(0.04,max(ylim),'B','FontSize',16,'FontWeight','bold')


%% Figure 8: Force redev - general effect of ATP and Pi

% load('final_fit.mat')
% x=final_params;


figure('Name', 'Force redev (Fig 8)')

w = 0.375;
h=0.4;
x1=0.1;
x2=0.6;
y1=0.575;
y2=0.1;

A= [0.1 1 5];
for i=1:3
    [F,t]=Rice_style_Fredev(@XBmodel_2024_Rat, x, [A(i) 1]);
    subplot('Position',[x1 y1 w h])
    hold on
    plot(t,F,'LineWidth',1)
    subplot('Position',[x1 y2 w h])
    hold on
    plot(t, (F-min(F))/(max(F)-min(F)),'LineWidth',1)
end
subplot('Position',[x1 y1 w h])
l=legend('0.1 mM ATP','1 mM ATP', '5 mM ATP','Location','southeast');
l.ItemTokenSize=[20,18];
ylabel('Stress (kPa)')
set(gca,'FontSize',11)
text(-0.2,max(ylim),'A','FontSize',16,'FontWeight','bold')


subplot('Position',[x1 y2 w h])
l=legend('0.1 mM ATP','1 mM ATP', '5 mM ATP','Location','southeast');
l.ItemTokenSize=[20,18];
ylabel('Stress (normalised)')
xlabel('Time (s)')
set(gca,'FontSize',11)
text(-0.2,max(ylim),'C','FontSize',16,'FontWeight','bold')

P= [1e-6 1 5];

for i=1:3
    [F,t]=Rice_style_Fredev(@XBmodel_2024_Rat, x, [5 P(i)]);
    subplot('Position',[x2 y1 w h])
    hold on
    plot(t,F,'LineWidth',1)
    subplot('Position', [x2 y2 w h])
    hold on
    plot(t, (F-min(F))/(max(F)-min(F)),'LineWidth',1)
end
subplot('Position', [x2 y1 w h])
Pi_legend={['0 mM ' P_i],['1 mM ' P_i], ['5 mM ' P_i]};
l=legend(Pi_legend,'Location','southeast');
l.ItemTokenSize=[20,18];
ylabel('Stress (kPa)')
set(gca,'FontSize',11)
text(-0.2,max(ylim),'B','FontSize',16,'FontWeight','bold')


subplot('Position', [x2 y2 w h])
l=legend(Pi_legend,'Location','southeast');
l.ItemTokenSize=[20,18];
ylabel('Stress (normalised)')
xlabel('Time (s)')
set(gca,'FontSize',11)
text(-0.2,max(ylim),'D','FontSize',16,'FontWeight','bold')

%% Figure A (breakdowns of final fit)

% load final parameters
load('final_fit.mat')

model=@XBmodel_2024_linear_perms;
mets={[5 1], [0.1 1], [1 1], [5 1e-6], [5 5]};

sd = [ 0 0 0 1 1];
md=4;
x=final_params;

% plotting details
colour={[0 0 0],blue,red,green,purple};
w=0.4;
h=0.175;
labels={'Combined','$H_\mathrm{xB}$','$H_\mathrm{xC}$','$H_\mathrm{CL}$','$H_\mathrm{CS}$'};
fs=logspace(-1,2,100);

figure('Units', 'normalized' ,'OuterPosition',  [0.3, 0.05, 0.4, 0.95],'Name','Breakdown of different transfer functions (Fig ');

for i=1:5
    [~,Yam,HxB_comp,HxC_comp,HC_comp,HC_Lcomp,HC_Scomp]=model(x,sd,md,mets{i});
    comps = {Yam, HxB_comp,HxC_comp,HC_Lcomp,HC_Scomp};
    for j=1:5
        yp=1-0.19*j;
        subplot('Position',[0.085 yp w h])
        semilogx(fs,real(comps{j}),'Color',colour{i},'LineWidth',2)
        hold on
        subplot('Position',[0.58 yp w h])
        semilogx(fs,imag(comps{j}),'Color',colour{i},'LineWidth',2)
        hold on

    end
end
[~,Yam,HxB_comp,HxC_comp,HC_comp,HC_Lcomp,HC_Scomp]=model(x,sd,md,mets{1});
    comps = {Yam, HxB_comp,HxC_comp,HC_Lcomp,HC_Scomp};

for j=1:5
    yp=1-0.19*j;
    subplot('Position',[0.085 yp w h])
    semilogx(fs,real(comps{j}),'Color',colour{1},'LineWidth',2)
    yma=max(ylim);
    set(gca,'FontSize',11,'XLim',[0.1 100])
    xticks([0.1 1 10 100])
    ylabel('Elastic Mod (MPa)','FontSize',12)
    text(0.03,yma,char(63+2*j),'FontSize',16,'FontWeight','bold')
    text(0.13,yma-0.1*(yma-min(ylim)),labels{j},'FontSize',12);
    if j==5 
    xlabel('Frequency (Hz)','FontSize',12); 
    xticklabels({'0.1' '1' '10' '100'})
    else
    xticklabels({})
    end

    subplot('Position',[0.58 yp w h])
    ylabel('Viscous Mod (MPa)','FontSize',12)
    if j==1
    legend('Baseline','0.1 mM ATP','1 mM ATP',['0 mM ' P_i],['5 mM ' P_i],'Location','northwest','AutoUpdate','Off')
    end
    semilogx(fs,imag(comps{j}),'Color',colour{1},'LineWidth',2)
    set(gca,'FontSize',11,'XLim',[0.1 100])
    xticks([0.1 1 10 100])
    
    text(0.03,max(ylim),char(64+2*j),'FontSize',16,'FontWeight','bold')
    if j==5 
    xlabel('Frequency (Hz)','FontSize',12); 
    xticklabels({'0.1' '1' '10' '100'})
    else
    xticklabels({})
    end

end


%% functions

% plotting and formatting of a 2 panel complex modulus plot
function []=two_panel_CMplot(fs,Y,style,thickness,colour)
    w = 0.4;
    h=0.8;
    
    subplot('Position',[0.075 0.15 w h])
    if nargin<=4
        semilogx(fs,real(Y),style,'LineWidth',thickness,'MarkerSize',10)
    else
        semilogx(fs,real(Y),style,'LineWidth',thickness,'MarkerSize',10,'Color',colour)
    end
    ylabel('Elastic Modulus (MPa)')
    xlim([0.1 100])
    xticklabels({'0.1' '1' '10' '100'})
    xlabel('Frequency (Hz)')
    hold on
    set(gca,'Fontsize',11)
    box off 
    
    set(gcf, 'Position',  [500, 300, 1000, 400])
    
    subplot('Position',[0.575 0.15 w h])
    if nargin<=4
        semilogx(fs,imag(Y),style,'LineWidth',thickness,'MarkerSize',10)
    else
        semilogx(fs,imag(Y),style,'LineWidth',thickness,'MarkerSize',10,'Color',colour)
    end
    ylabel('Viscous Modulus (MPa)')
    xlim([0.1 100])
    xticklabels({'0.1' '1' '10' '100'})
    xlabel('Frequency (Hz)')
    hold on
    set(gca,'Fontsize',11)
    box off

end

% set up for a four panel plot
function[f,A,B,C,D]=four_panel_plot(name)

f=figure('Name',name, 'Units', 'normalized' ,'OuterPosition', [0.3, 0.3, 0.4, 0.5]);
w = 0.4;
h=0.4;

subplot('Position',[0.075 0.575 w h])
hold on
ylabel('Elastic Modulus (MPa)','FontSize',12)
A=gca;
set(A,'Xscale','log','XLim',[0.1 101],'FontSize',11)
xticks([0.1 1 10 100])
xticklabels({'0.1' '1' '10' '100'})
box off

subplot('Position',[0.575 0.575 w h])
hold on
ylabel('Viscous Modulus (MPa)','FontSize',12)
B=gca;
set(B,'Xscale','log','XLim',[0.1 101],'FontSize',11)
xticklabels({'0.1' '1' '10' '100'})
box off

subplot('Position',[0.075 0.1 w h])
hold on
ylabel('Elastic Modulus (MPa)','FontSize',12)
C=gca;
set(C,'Xscale','log','XLim',[0.1 101],'FontSize',11)
xticklabels({'0.1' '1' '10' '100'})
xlabel('Frequency (Hz)','FontSize',12)
box off

subplot('Position',[0.575 0.1 w h])
hold on
ylabel('Viscous Modulus (MPa)','FontSize',12)
D=gca;
set(D,'Xscale','log','XLim',[0.1 101],'FontSize',11)
xticklabels({'0.1' '1' '10' '100'})
xlabel('Frequency (Hz)','FontSize',12)
box off
end

% Adds appropriate $s to row names for latex interpreting
function [latexed]=add_dollars(row_names)
for i=1:length(row_names)
    row=row_names{i};
    if ~strcmp(row,'None')
    new_row='$';
    ks=strfind(row,'k');
    if length(ks)>1
    spaces=strfind(row,' ');
    j=spaces(1);
    k=ks(2);
    new_row(2:length(row)+3)=[row(1:j-1) '$' row(j:k-1) '$' row(k:end)];
    else
    new_row(2:length(row)+1)=row(1:end);
    end
    new_row(length(new_row)+1)='$';
    row=new_row;
    end
    latexed{i}=row;
end
end

% Latex formating for Pi
function [Pi]=P_i()
    Pi='P$_{\mathrm{i}}$';
end


% Turns all plots to latex
function []=set_all_interpreters(interpreter)
%Inputs: 'none' for normal
%        'latex' for latexy
    list_factory = fieldnames(get(groot,'factory'));
    index_interpreter = find(contains(list_factory,'Interpreter'));
    for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,interpreter);
    end

end
