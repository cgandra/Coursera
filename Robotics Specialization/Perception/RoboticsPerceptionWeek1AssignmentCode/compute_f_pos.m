function [ f, pos ] = compute_f_pos( d1_ref, d2_ref, H1, H2, ratio, f_ref )
%% Compute camera focal length and camera position to achieve centain ratio between objects
%
% In this function, we focus on two objects: object A with height H1 and
% d1_ref as distance to camera in 3D world, object B with height H2 and
% d2_ref as distance to camera in 3D world.
% We will keep the size of object A in image the same as before while
% adjusting the size of object B in image.
%
% Input:
% - d1_ref: distance of the first object
% - d2_ref: distance of the second object
% - H1: height of the first object in physical world
% - H2: height of the second object in physical world
% - ratio: ratio between two objects in image coordinate (h1/h2)
% - f_ref: 1 by 1 double, previous camera focal length
% Output:
% - f: 1 by 1, camera focal length
% - pos: 1 by 1, camera position on z axis

% YOUR CODE HERE
% ratio = h1'/h2'
% But you don't have control over h1' and h2' directly 
% you have control over your position (which I call p in the diagram, pos in the HW) and your focal length f. 
% By similar triangles you have
% h1'/f = H1/(d1_ref-p)
% h2'/f = H2/(d2_ref-p)
% h2'*ratio = h1';
% H2*ratio/(d2_ref-p) = H1/(d1_ref-p)
% (H2*ratio)*(d1_ref - pos) = H1*(d2_ref - pos)
% pos = (H2*ratio*d1_ref-H1*d2_ref)/(H2*ratio-H1)

% (f * H1) / (d1_ref - pos) = (f_ref * H1) /d1_ref
% f = f_ref*(d1_ref - pos)/d1_ref

pos = ((H2*ratio*d1_ref)-(H1*d2_ref))/((H2*ratio)-H1);
f = f_ref*(d1_ref-pos)/d1_ref;

end

