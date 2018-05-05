function [LRM_j,LRM_uni] = LRM_growthrate(LRM,p,restriction)

% tbeucler - 2017
% Computes column-integrated growth rate for point-wise perturbation
% and uniform perturbation in s^-1
%%% INPUTS
% LRM: Linear response matrix [day^-1] in the form LRM_{ji} !!!
% where i is the response level and j is the perturbation level
% p: Pressure [Pa or hPa]
% restriction: Restriction domain for LRM [index array]
%%% OUTPUTS
% LRM_j: Col-int growth rate for point-wise perturbation at level j [day^-1]
% LRM_uni: Col-int growth rate for uniform perturbation [day^-1]

res=restriction; 
% Case where the restriction domain for LRM starts at 1
if restriction(1)==1, restriction=restriction(2:end); end
% General case
L_restriction = numel(restriction);
LRM_j = zeros(L_restriction,1);
% Sum (Mij*delPi from i=1 to i=N)/dePj
for ires=1:L_restriction
    ind = restriction(ires);
    LRM_j(ires)=sum((LRM(ind,restriction)).*...
        (p(restriction+1)-p(restriction-1)))/(p(ind+1)-p(ind-1));
end
% Sum (Mij*delPi from i,j=1,1 to i,j=N,N)/Sum(delPi from i,j=1,1 to i,j=N,N)
LRM_uni = sum(sum((LRM(ind,restriction)).*...
    (p(restriction+1)-p(restriction-1))))/...
    sum(p(restriction+1)-p(restriction-1));
% Case where the restriction domain for LRM starts at 1
if res(1)==1, LRM_j = cat(1,LRM(1,1)+0.5*sum((LRM(1,restriction)).*...
        (p(restriction+1)-p(restriction-1)))/(p(2)-p(1)),LRM_j);
end