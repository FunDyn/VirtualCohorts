function SC = SC_MFM(FC)

% ----------------------------------------------------------------------- %
% <<<<<<<<<<<<< code for non-linear FC-to-SC completion >>>>>>>>>>>>>>>>> %
% ----------------------------------------------------------------------- %

% -> This function takes as input the empirical functional connectivity FC
%    and simulate the MFM structure connectivity SC

% -> input: TS --> empirical time-series
% -> output: SC --> non-linear virtual SC (SC_MFM)


rng('shuffle')

n = size(FC,2);          % determine the number of regions

FCemp = FC/max(max(FC));   % normalize empirical input FC

% build SC_0: intial random matrix
sc0 = randn(n);
sc0 = sc0-sc0.*eye(n);
sc0 = abs(sc0);
sc0 = (sc0+sc0')/2;
sc0 = sc0/max(max(sc0));

% train the sc0
for iter=1:1000

        FCsim = FC_MFM(sc0);
        FCsim = FCsim/max(max(FCsim));
        D = FCemp-FCsim;
        cc = corrcoef(FCemp,FCsim)
        
        sc0 = sc0+0.001.*D;
        sc0(sc0<0) = 0;
        sc0 = (sc0+sc0')/2;
        sc0 = sc0/max(max(sc0));
        
        if cc(2) > .5
            SC = sc0; % return SC_MFM as output
            break
        end
end
end

