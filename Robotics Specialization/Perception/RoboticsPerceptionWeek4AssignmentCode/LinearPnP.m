function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly

% focalLength    = [568.9961, 568.9884];
% imageSize      = [960, 1280];
% principalPoint = [643.2106, 477.9828];
% intrinsics = cameraIntrinsics(focalLength,principalPoint,imageSize)
% [R,T] = estimateWorldCameraPose(x, X, intrinsics)

n = size(x, 1);
A=zeros(3*n,12);
% for i=1:n
%     A(2 * i - 1, :) = [X(i, 1) X(i, 2) X(i, 3) 1 0 0 0 0 -x(i, 1)*X(i, 1) -x(i, 1)*X(i, 1) -x(i, 1)*X(i, 3) -x(i, 1)];
%     A(2 * i, :)     = [0 0 0 0 X(i, 1) X(i, 2) X(i, 3) 1 -x(i, 2)*X(i, 1) -x(i, 2)*X(i, 1) -x(i, 2)*X(i, 3) -x(i, 2)];
% end

x = [x,ones(size(x,1),1)];
X = [X,ones(size(X,1),1)];
% Convert image points to normalized image coordinates
xn = (inv(K)*x')';
%xn = x/K;
xn = xn./xn(:, end);
u = xn(:, 1);
v = xn(:, 2);
n = size(X, 1);
A = [zeros(n,4), -X, v.*X
      X, zeros(n,4), -u.*X
      -v.*X, u.*X, zeros(n,4)];
 
% for i = 1 : n
%      Xt = X(i,:);
%      j = (i-1)*3 + 1;
%      A(j:j+2,:) = [zeros(1,4), -Xt, v(i)*Xt;
%                    Xt, zeros(1,4), -u(i)*Xt;
%                    -v(i)*Xt, u(i)*Xt, zeros(1,4)];
% end

[U,S,V] = svd(A);
P = reshape(V(:,end), 4, 3)';


% % Estimate C.
% P1 = P * sign(det(P(1:3, 1:3)));
% [U,S,V] = svd(P1);
% C1 = V(1:3,end) / V(end, end);
% 
% % Estimating K and R by RQ decomposition.
% [K1, R1] = rq(P1(1:3,1:3));
% % Enforce positive diagonal of K by changing signs in both K and R.
% D = diag(sign(diag(K1)));
% K1 = K1*D;
% R1 = D*R1;
% K1 = K1/K1(end,end);
% %Determine t from the estimated C.
% T1 = -R1*C1;

R = P(:,1:3);
[U,S,V] = svd(R);
R = U*V';
T = P(:,4)./S(1,1);
if det(R) < 0
    R =- R;
    T =- T;
end
C = -R' * T;

end

function [R,Q] = rq(M)
[Q,R] = qr(rot90(M,3));
R = rot90(R,2)';
Q = rot90(Q);
end