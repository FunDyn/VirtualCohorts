function FC = FC_MFM(SC)

% ----------------------------------------------------------------------- %
% <<<<<<<<<<<<< code for non-linear SC-to-FC completion >>>>>>>>>>>>>>>>> %
% ----------------------------------------------------------------------- %

% -> This function takes as input the empirical structure connectivity SC
%    and simulate the MFM functional connectivity FC 

% -> input: SC --> empirical structural connectivity
% -> output: FC --> non-linear virtual FC (FC_MFM)


rng('shuffle')

% local parameters of Wong-Wang neural mass model
gamma = .641/1000;
JN = .2609;

% global parameters of Wong-Wang neural mass model
W = .9;
I0 = .32;
G = 2;
tao = 20;

n = size(SC,2);          % determine the number of regions
noise = .001;            % noise value

% the simulation time parameters
dt = .1;
tmax = 101000;
tspan = dt:dt:tmax;
time = length(tspan);

NA = zeros(time,n); % preallocate the neural activity NA

ts(:,1)=.3*rand(n,1)+.2;  % initial condition

for k=1:time
    x = I0 + JN*(W*ts + G*(SC*ts));
    r = local_MFM(x);
    ts = ts + dt*(-ts/tao+(1-ts)*gamma.*r) + (sqrt(2*dt)*noise*randn(n,1));
    ts(ts>1) = 1;
    ts=ts.*(ts>=0);
    NA(k,:)=ts';
end
clear y x r k


NA=NA(1:10:end,:);   % downsampling of simulated neural activity

for region=n:-1:1
    B = BOLD_MFM(NA(:,region));
    Bold(:,region) = B;
end

Bold=Bold(1:100:end,:);     % downsampling of simulated bold

FC = corr(Bold);        % return FC_MFM as output


end

