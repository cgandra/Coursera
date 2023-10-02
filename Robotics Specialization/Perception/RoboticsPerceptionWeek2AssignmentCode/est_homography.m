function [ H ] = est_homography(video_pts, logo_pts)
% est_homography estimates the homography to transform each of the
% video_pts into the logo_pts
% Inputs:
%     video_pts: a 4x2 matrix of corner points in the video
%     logo_pts: a 4x2 matrix of logo points that correspond to video_pts
% Outputs:
%     H: a 3x3 homography matrix such that logo_pts ~ H*video_pts
% Written for the University of Pennsylvania's Robotics:Perception course

% YOUR CODE HERE
H = [];
vp = video_pts;
lp = logo_pts;

A = zeros(8, 9);
for i = 1:length(video_pts)    
    A(2 * i - 1, :) = [-vp(i,1), -vp(i,2), -1, 0, 0, 0, vp(i,1)*lp(i,1), vp(i,2)*lp(i,1), lp(i,1)];
    A(2 * i, :) = [0, 0, 0, -vp(i,1), -vp(i,2), -1, vp(i,1)*lp(i,2), vp(i,2)*lp(i,2), lp(i,2)];
end

[U, S, V] = svd(A);
h = V(:,end);
H = reshape(h,[3,3])';
H = H./H(3, 3);
end

