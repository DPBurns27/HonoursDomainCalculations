function [Lat, Long] = squareDomain(centreLat, centreLong, widthdeg, subdivisions)

% Creates a pair of vectors containing a square calculated from it's centre
% point and width

% Inputs:

% centreLat     = lattitude of the squares centre 
% centreLong    = longitude of the squares centre
% widthdeg      = width of the square in degrees
% subdivision   = the number of subdivisions along each axis

% Calculates the longitude points    
topLong = linspace((centreLong-widthdeg/2),(centreLong+widthdeg/2),subdivisions+2);
rightLong = linspace(centreLong+widthdeg/2,centreLong+widthdeg/2, subdivisions);
botLong = fliplr(topLong);
leftLong = linspace(centreLong-widthdeg/2,centreLong-widthdeg/2, subdivisions);

% Calculates the lattitude points
topLat = linspace(centreLat+widthdeg/2, centreLat+widthdeg/2, subdivisions+2);
rightLat = linspace(centreLat+widthdeg/2, centreLat-widthdeg/2, subdivisions+2);
rightLat = rightLat(2:end-1);
botLat = linspace(centreLat-widthdeg/2, centreLat-widthdeg/2, subdivisions+2);
leftLat = fliplr(rightLat);

% Stitches the Lat and Long points together to form the square
Long = [topLong, rightLong, botLong, leftLong];
Long = [Long,Long(1)];
Lat = [topLat, rightLat, botLat, leftLat];
Lat = [Lat,Lat(1)];
