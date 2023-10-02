function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
%% Nonlinear_Triangulation
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     x
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations 

P1 = K * R1 * [eye(3) -C1];
P2 = K * R2 * [eye(3) -C2];
P3 = K * R3 * [eye(3) -C3];
syms x y z 
x_3d = [x;y;z;1];
x_c1 = P1*x_3d;
x_c1n = x_c1./x_c1(end);
x_c2 = P2*x_3d;
x_c2n = x_c2./x_c2(end);
x_c3 = P3*x_3d;
x_c3n = x_c3./x_c3(end);
    
x_c = [x_c1n(1:2);x_c2n(1:2);x_c3n(1:2)];
J = jacobian(x_c,[x,y,z]);

%X = [X0,ones(size(X0,1),1)];
X = X0;
b = [x1,x2,x3];

for i = 1:size(X,1)
    err = 1;
    iteration = 0;
    x_c_r = eval(subs(x_c,[x,y,z],X(i,:)));
    while err > 1e-05 && iteration < 100
        J_r = eval(subs(J,[x,y,z],X(i,:)));  
        delta_x = (J_r'*J_r)\J_r'*(b(i,:)'-x_c_r);
        X(i,:) = X(i,:) + delta_x';
        x_c_r_old = x_c_r;
        x_c_r     = eval(subs(x_c,[x,y,z],X(i,:)));
        err = norm(x_c_r - x_c_r_old)/norm(x_c_r_old);
        iteration = iteration + 1;
    end
end

end
