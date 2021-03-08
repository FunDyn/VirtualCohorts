function H=local_MFM(I)

% ----------------------------------------------------------------------- %
% <<<<<<<<<<<<<<< complementary code for MFM algorithm >>>>>>>>>>>>>>>>>> %
% ----------------------------------------------------------------------- %

% local parameters of Wong-Wang neural mass model in Hz
a=270.0;
b=108.0;
c=0.154; 

x=a*I-b;

if(any(abs(c*x)<1e-6))
    indices = abs(c*x)<1e-6;
    H = zeros(length(I),1);
    H(indices) = .5*x(indices)+1.0/c;
    H(~indices) = x(~indices)./(1-exp(-c*x(~indices)));
else
    H=x./(1-exp(-c*x));
end
