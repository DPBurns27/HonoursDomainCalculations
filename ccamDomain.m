function [finalLengthd, finalLengthk, finalResolution, gridClean, schmidtClean] = ccamDomain(centreLat,width,resolution)
% ccamDomain: This script takes in a number of variables relating to the
% creation of a domain that will be used in ccam. It outputs the necessary
% variables to enter into ccam, along with the changes made to the starting
% values so that the values given to ccam are nice (rounded)

% Inputs: 
% centreLat     =   lattitude of domains centre point
% centreLong    =   longitude of domains centre point
% width         =   width of the domain in km
% resolution    =   resolution of the domain in km

% Outputs:

% finalLengthd      =   final length of domain in degrees
% finalLengthk      =   final length of domain in km
% finalResolution   =   final resolution of domain in km
% gridClean         =   grid size selected (multiple of 24)
% schmidtClean      =   rounded schmidt number
% widthdeg          =   the length in degrees before converting to ccam
%                       style

% find phi, the conversion factor between km and degs

radiusEarth = 6371; %km

% This equation comes from a bunch of maths I did to expose the dependance
% of the degree to km conversion on longitude. I considered the circle
% cross sections at the equator and at the centre point, and thus the
% conversion between their radii
phi = 2*pi*radiusEarth*cos(centreLat*pi/180)/360;

% Ok now we need to take the width and then convert back through this
% conversion factor to get the degrees
widthdeg = width/phi;

% Next lets figure out the variables that ccam wants

% ccam needs only a few variables to create it's domain
% long, lat: the centre point of the domain
% schmidt: the schmidt number of the domain, related to the size of the
% front face of the cube
% gridsize: number of grid points horizontally or vertically

% So then you pick a domain of whatever size and position you want and
% figure out the size, and schmidt number for it, you have to work back
% from the resolution and grid size to find the schmidt number

% Marcus has given me 2 equations for resolution and length in degrees
% to calculate the 2 variables ccam needs schmidt and gridsize
% Here I have rearranged them to expose the variables needed by ccam

schmidt = widthdeg/90;
gridsize = phi*widthdeg/resolution;

% But Marcus says that arbitrary values for gridsize and schmidt number may
% result in bad performance by ccam so lets find the right values for them

% preffered gridsizes
prefGrid = [0, 48, 72, 96, 144, 192, 288, 384, 576, 768, 1152, 1536];

% better to round the schmidt number
schmidtClean = round(schmidt, 1, 'significant');


% take gridsize away from a vector containing the ok grid sizes, then find
% the smallest one
 [M,I] = min(abs(prefGrid - gridsize));

%maybe I should just take the first positive value so I always try for a
%finer resolution, i also think that im going the wrong way here, i should
%be rounding down
 
% Is the chosen gridsize quite close (within 10%?) to the calculated gridsize?
if ((prefGrid(I) - gridsize) < 0) && (M > (prefGrid(I+1) - prefGrid(I))*0.2)
    % Its not close enough so shift it up to the next one 
    gridClean = prefGrid(I+1);
else
    % All good lets use it
    gridClean = prefGrid(I);
end
    
% Ok now that we have decided on some new values for schmidt and gridsize,
% what are my new values for length and resolution?

finalLength = 90*schmidtClean;
finalResolution = phi*finalLength/gridClean;

finalLengthk = 90*phi*schmidtClean;
finalLengthd = 90*schmidtClean;

% I dunno maybe ive approached this wrong, perhaps I should be generating
% vectors of these variables and trying to find the combination that gives
% the closest value to the original resolution



