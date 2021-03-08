function SC = SC_SLM(TS)

% ----------------------------------------------------------------------- %
% <<<<<<<<<<<<<< code for linear FC-to-SC completion >>>>>>>>>>>>>>>>>>>> %
% ----------------------------------------------------------------------- %

% -> This function takes as input the empirical time-series TS and infere 
%    the SLM structure connectivity SC

% -> input: TS --> empirical TS, columns are regions, rows are time points
% -> output: SC --> linear virtual SC (SC_SLM)


n = size(TS,2);         % determine the number of regions

COV = cov(TS);          % compute the covariance matrix from time-series

X = -inv(COV);          % evaluate the SLM structure connectivity

X = X - X.*eye(n);      % remove its diagonal

X(X<0) = 0;             % set negative values to zero

X = X/norm(X);          % normalize SC matrix

SC = X;                 % return the output


end

