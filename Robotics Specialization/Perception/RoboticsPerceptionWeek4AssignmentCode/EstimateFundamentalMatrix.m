function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2

x1 =  [x1 ones(size(x1,1),1)];
x2 =  [x2 ones(size(x2,1),1)];
T1 = normalization_transform(x1);
T2 = normalization_transform(x2);

% Normalize the correspondences
x1 = (T1*x1')';
x2 = (T2*x2')';

n = size(x1, 1);
% A = zeros(n,9);
% for i=1:n
%     A(i,:)= [x1(i,1)*x2(i,1), x1(i,1)*x2(i,2), x1(i,1), x1(i,2)*x2(i,1), x1(i,2)*x2(i,2), x1(i,2), x2(i,1), x2(i,2), 1]; 
% end

A = kron(x2, [1,1,1]).*[x1, x1, x1];

[U1, S1, V1] = svd(A);
F_unconditioned = reshape(V1(:,end),[3,3])';

[U,S,V] = svd(F_unconditioned);
S(end, end)=0;
F_conditioned = U * S * V';

% Denormalizing the Fundamental Matrix
F = T2'*F_conditioned*T1;
F = F/F(end,end);
end

function [T] = normalization_transform(x)
 % Find mean of the interest points
 mu = mean(x);
 % Find the distance of the interest points from their mean
 difference = x-repmat(mu,size(x,1),1);
 distance = sqrt(difference(:,1).^2+difference(:,2).^2);
 % Normalize all the interest points using the mean and distance
 distance_mean = mean(distance);
 scale = sqrt(2)/distance_mean;
 x_translate = -scale*mu(1);
 y_translate = -scale*mu(2);
 % Find the transformation for normalizing points
 T = [scale 0 x_translate; 0 scale y_translate ; 0 0 1];
end
