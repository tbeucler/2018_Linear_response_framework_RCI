%% RCI_lead_eig_real_part
% tbeucler - 4/11/2018
% Computes the leading eigenvalue real part and the complex modulus
% of the corresponding eigenvector for a given matrix M

function [lambda,eigenvector] = RCI_lead_eig_real_part(M)

[V,D]=eig(M); % V has the eigenvectors and D is the diagonal matrix 
% containing the eigenvalues
[l,i]=sort(real(diag(D))); lambda = l(end);
eigenvector = ((real(V(:,i(end)))).^2+(imag(V(:,i(end)))).^2).^0.5;

end