function route = GradientBasedPlanner (f, start_coords, end_coords, max_its)
% GradientBasedPlanner : This function plans a path through a 2D
% environment from a start to a destination based on the gradient of the
% function f which is passed in as a 2D array. The two arguments
% start_coords and end_coords denote the coordinates of the start and end
% positions respectively in the array while max_its indicates an upper
% bound on the number of iterations that the system can use before giving
% up.
% The output, route, is an array with 2 columns and n rows where the rows
% correspond to the coordinates of the robot as it moves along the route.
% The first column corresponds to the x coordinate and the second to the y coordinate

[gx, gy] = gradient (-f);

%%% All of your code should be between the two lines of stars.
% *******************************************************************
start_coords = double(start_coords);
end_coords = double(end_coords);

route = start_coords;
its = 1;
sc = start_coords;
while(its<=max_its && norm(end_coords-sc)>=2)
    scr = round(sc);
    temp = double([gx(scr(2), scr(1)), gy(scr(2), scr(1))]);
    sc = sc+(temp/norm(temp));
    route = [route;sc];
    its = its+1;
end
% *******************************************************************
end
