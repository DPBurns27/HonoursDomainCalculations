% Script to generate the numbers needed for the first CTM run over
% brisbane

close all
clc
clear

%Location off the coast near Bowen
%-19.624822, 148.355649

%% Initialise values

% Initial domains
names = {'gbr1','gbr2','gbr3'};
resolution = [80;9;3];
initLengthk = [6000,2000,1000];
centreLat = -19.62;
centreLong = 148.36;
numDomains = length(resolution);

% Plot axis buffer
buffer = 10; %deg

% Size of plot ticks
ticksize = 10; %deg

% Produce google maps files?
gmaps = false;

% Dump out lines for scripts?
scriptlines = false;

%% Calculate length in degrees

radiusEarth = 6371; %km

% This equation comes from a bunch of maths I did to expose the dependance
% of the degree to km conversion on longitude. I considered the circle
% cross sections at the equator and at the centre point, and thus the
% conversion between their radii
phi = 2*pi*radiusEarth*cos(centreLat*pi/180)/360;

% Ok now we need to take the length and then convert back through this
% conversion factor to get the degrees

initLengthd = initLengthk./phi;

%% Load in the high res coastline data

% Opens up the gshhs file to import high definition coastal data: _i is
% medium res, _f is full resolution
latlim = [centreLat-max(initLengthd)/2-buffer centreLat+max(initLengthd)/2+buffer];
lonlim = [centreLong-max(initLengthd)/2-buffer centreLong+max(initLengthd)/2+buffer];
S = gshhs('gshhs_i.b', latlim, lonlim);
% Pulls out the different levels of the map ie, coastline, rivers, etc
levels = [S.Level];
% Picks out just the coastal data
L1 = S(levels == 1);

%% Plot the coastline, the GBR and the domain centre

figure

% Sets the projection type, and axis limits of the map
axesm('mercator', 'MapLatLimit', latlim, 'MapLonLimit',lonlim)
%axesm('robinson', 'MapLatLimit', latlim, 'MapLonLimit',lonlim)
gridm; mlabel; plabel

% Changes the distance between subdivisions of the grid
setm(gca, 'MLineLocation', ticksize, 'MLabelLocation', ticksize, 'PLineLocation', ticksize, 'PLabelLocation', ticksize)
        
% Plots the coastal data on the map
geoshow([L1.Lat], [L1.Lon], 'Color', 'black')
        
% Draw the GBR in for reference
[GBRlat,GBRlong] = GBRCoords();
% Cut off the coastline section of the reef
GBRlong = GBRlong(1:8);
GBRlat = GBRlat(1:8);

% Plots the GBR on the map
geoshow(GBRlat, GBRlong, 'Color', 'red')

% Plots the centre point
cpointPlot = geoshow(centreLat, centreLong, 'DisplayType', 'Point', 'Marker', '+', 'Color', 'black', 'MarkerEdgeColor', 'black');


%% Plot the initial domains
for i = 1:numDomains
    
    % Calculate the square shape
    [initSquareLat, initSquareLong] = squareDomain(centreLat, centreLong, initLengthd(i), 0);
    % Plot the square
    initPlot = geoshow(initSquareLat, initSquareLong, 'DisplayType', 'Line', 'Color', 'red');
    
    % Produce shape for google maps
    if gmaps == true
        % Calculate the square shape but with many subdivisions
        [Llat, Llong] = squareDomain(centreLat, centreLong, initLengthd(i), 20);
        % output a kml file for use in google earth
        % need to change the directory its going to documents/matlab
        shape = geoshape(Llat,Llong);
        %kmlwrite(filename,S,'Name',S.Name,'FaceColor',colors,'EdgeColor','k')
        kmlwrite('InitialDomain',shape,'Color','r')
    end
end

%% Calculate the ccam domains from the initial domains

% Predeclare the output variables for ccamDomain
ccamLengthd = zeros(numDomains,1);
ccamLengthk = zeros(numDomains,1);
ccamResolution = zeros(numDomains,1);
gridsizeClean = zeros(numDomains,1);
schmidtClean = zeros(numDomains,1);

for i = 1:numDomains
    
    % Calculates the appropriate ccam values from the initial domain values
    [ccamLengthd(i), ccamLengthk(i), ccamResolution(i), gridsizeClean(i), schmidtClean(i)] = ...
        ccamDomain(centreLat,initLengthk(i),resolution(i));

    % Calculate the ccam resolution in degrees
    ccamResolutiond = ccamLengthd./gridsizeClean;

    % Calculate the square shape
    [ccamSquareLat, ccamSquareLong] = squareDomain(centreLat, centreLong, ccamLengthd(i), 0);
    % Plot the square
    ccamPlot = geoshow(ccamSquareLat, ccamSquareLong, 'DisplayType', 'Line', 'Color', 'blue');
    
    % Produce shape for google maps
    if gmaps == true
        % Calculate the square shape but with many subdivisions
        [Llat, Llong] = squareDomain(centreLat, centreLong, ccamLengthd(i), 20);
        % Output a kml file for use in google earth
        % need to change the directory its going to documents/matlab
        shape = geoshape(Llat,Llong);
        kmlwrite('InitialDomain',shape,'Color','b')
    end
    
    %% Dump out ccam script lines
    
    if scriptlines == true
        % Dumps out the values that need to be changed in the ccam script
        disp('--------------------------');
        disp(names{i});
        
        minlat = centreLat - ccamLengthd(i)/2;
        maxlat = centreLat + ccamLengthd(i)/2;
        minlong = centreLong - ccamLengthd(i)/2;
        maxlong = centreLong + ccamLengthd(i)/2;
        
        Sminlat=sprintf('minlat=%s                                ;# output min latitude (degrees)',num2str(minlat));
        disp(Sminlat)
        Smaxlat=sprintf('maxlat=%s                                ;# output max latitude (degrees)',num2str(maxlat));
        disp(Smaxlat)
        Sminlong=sprintf('minlon=%s                               ;# output min longitude (degrees)',num2str(minlong));
        disp(Sminlong)
        Smaxlong=sprintf('maxlon=%s                               ;# output max longitude (degrees)',num2str(maxlong));
        disp(Smaxlong)
        
        disp(strcat('name=ccam_',num2str(resolution(i)),'km'));
        disp(strcat('nproc=',num2str(gridsizeClean(i))));
        ntaskspnode = ceil(gridsizeClean(i)/10)*10;
        disp(strcat('#SBATCH --ntasks-per-node=',num2str(ntaskspnode)));
        
        disp(strcat('domain=',num2str(gridsizeClean(i)),'_',num2str(centreLong),'_',num2str(centreLat),'_',num2str(schmidtClean(i))));
        
        if i == 1
            disp('dmode=0');
            disp('bcdom=ccam_eraint_');
            disp('bcdir=/datastore/wwwusers/244528/eraint')
        else
            disp('dmode=2');
            disp(strcat('bcdom=ccam_',num2str(resolution(i-1)),'km'));
            disp(strcat('bcdir=$insdir/scripts/brisruns',num2str(gridsizeClean(i-1)),'/OUTPUT'));
        end
    end
    
end


%% Calculate the ctm domains from the ccam domains

% Predeclare the output variables for ctmDomain
blLat = zeros(numDomains,1);
blLong = zeros(numDomains,1);
dx = zeros(numDomains,1);
dy = zeros(numDomains,1);
nx = zeros(numDomains,1);
ny = zeros(numDomains,1);

for i = 1:numDomains
    
    % Calculate the ctm domain from the ccam domain
    [blLat(i), blLong(i), dx(i), dy(i), nx(i), ny(i), X2, Y2] = ctmDomain(ccamResolutiond(i), ccamLengthd(i), ccamLengthd(i), centreLat, centreLong);
    
    % Convert the output of ctmdomain to the values needed for squareDomain
    ctmLengthd = dx(i)*nx(i);
    ctmcentreLong = blLong(i) - dx(i)/2 + ctmLengthd/2;
    ctmcentreLat = blLat(i) - dy(i)/2 + ctmLengthd/2;
    
    % Calculate the square shape
    [ctmSquareLat, ctmSquareLong] = squareDomain(ctmcentreLat, ctmcentreLong, ctmLengthd, 0);
    % Plot the square
    ctmPlot = geoshow(ctmSquareLat, ctmSquareLong, 'DisplayType', 'Line', 'Color', 'green');
    
    % Produce shape for google maps
    if gmaps == true
        % Calculate the square shape but with many subdivisions
        [Llat, Llong] = squareDomain(ctmcentreLat, ctmcentreLong, ctmLengthd, 20);
        % Output a kml file for use in google earth
        % need to change the directory its going to documents/matlab
        shape = geoshape(Llat,Llong);
        kmlwrite('InitialDomain',shape,'Color','g')
    end
    
    %% Dump out lines for ctm script files
    
    if scriptlines == true
        
        % Output for the j_prep file
        disp('--------------------------');
        
        disp('Position of SW corner of each grid: xlong,ylat');
        disp([num2str(blLong(i)), ' ', num2str(blLat(i))]);
        
        disp('EACH GRID HORIZONTAL DIMENSION: NX,NY');
        disp([num2str(nx(i)), ' ', num2str(ny(i))]);
        
        % Maybe i need the full length of the resolution?
        disp('EACH GRID HORIZONTAL RESOLUTION: HRESX,HRESY (deg) , for *NC naming');
        disp([num2str(dx(i)), ' ', num2str(dy(i))]);
        
        % Output for the emission scripts
        disp('--------------------------');
        % Remember that this value is actually in degrees not metres so
        % times by 100 as well as 1000, could possibly use phi instead
        disp(strcat('set airshed=', names{i}));
        disp(strcat('set cellsize=',num2str(dx(i)*100*1000),'m'));
        disp(strcat('set xw=',num2str(blLong(i))));
        disp(strcat('set ys=',num2str(blLat(i))));
        disp(strcat('set dxy=',num2str(dx(i))));
        disp(strcat('set nx=',num2str(nx(i))));
        disp(strcat('set ny=',num2str(ny(i))));
            
    end
    
end

% add in a legend
legend([initPlot, ccamPlot, ctmPlot, cpointPlot],{'Initial Domain','CCAM Domain','CTM Domain', '(-19.62, 148.36)'}, 'Location', 'southeast')

% Tightens up the axis to whatever is plotted
tightmap




