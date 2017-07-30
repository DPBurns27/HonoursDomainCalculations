function [blLat, blLong, dx, dy, nx, ny, X, Y] = ctmDomain(resolution, width, height, cLat, cLong)

% ctmDomain: Takes in some information about a domain and returns the
% format of the domain needed by CTM

% Inputs:

% resolution    =   grid size 
% width         =   width of the domain
% height        =   height of the domain
% cLat          =   centre lattitude of the domain
% cLong         =   centre longitude of the domain

% Outputs:

% blLat     =   The lattitude of the bottom left grid square's centre 
% blLong    =   The longitude of the bottom left grid square's centre 
% dx        =   The length of each grid square (units?)
% dy        =   The height of each grid square (units?)
% nx        =   The number of horizontal grid squares
% ny        =   The number of vertical grid squares

dx = round(resolution*3.2,2,'significant');
dy = round(resolution*3.2,2,'significant');

% I'm scaling down the domain here so that i try and miss the singularities
% at the corners, but I don't think it is enough to just reduce the grid
% points by a constant value, i think it had to be some function of the
% grid width, maybe i can just make it scale with nx

nx = ceil(width*0.9/dx);
ny = ceil(height*0.9/dx);

if nx < 50
   nx = 50;
   ny = 50;
   dx = round(width*0.9/nx,2,'significant');
   dy = round(height*0.9/ny,2,'significant');
end

% The full length of the domain is nx*dx, but the length between the end
% nodes is (nx-1)*dx

blLong = cLong - dx*(nx - 1)/2;
blLat = cLat - dx*(nx - 1)/2;

% ok now actually make the grid from these values

x = blLong+dx*(0:nx);
y = blLat+dx*(0:ny);
[X,Y] = meshgrid(x,y);



