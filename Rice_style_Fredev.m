%

% Simulates force redevelopment in a model as following rapid length step
% in the same way as described by Rice 2008
% Redevelopment is found by speeding up detachment
% and lowering attachment - 500x each. 
%
%
% Inputs:
%       - odes: function handle for the ODE model
%       - params: fitting parameters used in the model
%       - L_final: % of optimal muscle length following the step (eg. 0.95)
% Outputs:
%       - F: output stress of the model as force is redeveloped (kPa)
%
% Author: Julia Musgrave
% Date last updated: Aug 2023

function [F,t]=Rice_style_Fredev(odes,params,mets)

% ICs
yi=odes();
L0=2.2; % experimental muscle length

%solving ode to steady state (at optimal length)
tspan=[0 1];
prev=0;
curr=100;
options=odeset('RelTol',1e-6,'Abstol',1e-6,'MaxStep',0.001);
while abs((prev-curr)/curr)>1e-5
[t,y]=ode15s(@(t,y)odes(t,y,L0,params,mets),tspan,yi,options);
y0=y(end,:);
prev=curr;
[~,curr]=odes(t(end),y0,L0,params,mets);
end

% simulating the force drop
drop_params=params;
drop_params(1)=params(1)/500;
drop_params([2 5])=params([2 5])*500;
drop_params(3:4)=params(3:4)*10;
options=odeset('RelTol',1e-6,'Abstol',1e-6,'MaxStep',0.001);
[~,ypre]=ode15s(@(t,y)odes(t,y,L0,drop_params,[1 1]),[0 0.002],y0,options);


% simulating force recovery
[t,y]=ode15s(@(t,y)odes(t,y,L0,params,mets),[0 1],ypre(end,:),options);
[~,F]=odes(t,y,L0,params,mets);

end
