function [ A_B ] = differentiate( A,B )
% differentiate A_B = delta(A)/delta(B)
% along their first dimension
% Assumes two dimensions for now

S_A = size(A);S_B = size(B);

if S_A ~= S_B
    disp('A and B must be the same size');
    return
end

A_B = zeros(S_A(1),S_A(2));

for i2 = 1:S_A(2)
    
    A_B(1,i2) = (A(2,i2)-A(1,i2))/(B(2,i2)-B(1,i2));
    for i1 = 2:(S_A(1)-1)
        A_B(i1,i2) = (A(i1+1,i2)-A(i1-1,i2))/(B(i1+1,i2)-B(i1-1,i2));
    end
    A_B(S_A(1),i2) = (A(S_A(1),i2)-A(S_A(1)-1,i2))/(B(S_A(1),i2)-B(S_A(1)-1,i2));
    
end

end

