function FC = FC_SLM(SC)

% ----------------------------------------------------------------------- %
% <<<<<<<<<<<<<< code for linear SC-to-FC completion >>>>>>>>>>>>>>>>>>>> %
% ----------------------------------------------------------------------- %

% -> This function takes as input the empirical structure connectivity SC 
%    and infere the SLM functional connectivity FC

% -> input: SC --> empirical structural connectivity
% -> output: FC --> linear virtual FC (FC_SLM)


n = size(SC,2);         % determine the number of regions

SC = SC - SC.*eye(n);   % remove the diagonal and normalize the input 
SC = SC/norm(SC);       

G_ref = 0.83;           % the reference value of global scale 

X = G_ref*SC;           % coupling matrix for the linear system

Xcov = -inv(X);         % virtual covariance matrix

FC = zeros(n);          % preallocate the output


for i=1:n
    for j=1:n
        FC(i,j)=Xcov(i,j)/sqrt(Xcov(i,i)*Xcov(j,j));  % compute the SLM FC
    end
end


end