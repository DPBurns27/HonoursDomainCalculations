function [gbrlong, gbrlat] = GBRCoords()

% Get the full earth coast coords (loads in as vectors lat and long)
load coast;

% Eastern points of the gbr in degrees
% These values are taken from http://www.gbrmpa.gov.au/__data/assets/pdf_file/0019/4906/mp_009_full.pdf
eastgbrlat = [10,49,48; 12,59,55; 14,59,55; 17,29,55; 20,59,54; 24,29,54; 24,29,54];
eastgbrlong = [145,00,04; 145,00,04; 146,00,04; 147,00,04; 152,55,04; 154,00,04; 152,6,00];

% Eastern points of the gbr in decimal
eastbgrlatdec = (eastgbrlat(:,1) + eastgbrlat(:,2)./60 + eastgbrlat(:,3)./60./100).*(-1);
eastbgrlongdec = eastgbrlong(:,1) + eastgbrlong(:,2)./60 + eastgbrlong(:,3)./60./100;

% Get the western points of gbr from the world map
refcoast = lat <= eastbgrlatdec(1) & lat >= eastbgrlatdec(end)...
    & long >= 142.4 & long <= eastbgrlongdec(end);
westgbrlat = lat(refcoast);
westgbrlong = long(refcoast);

% Combine the east and west points of gbr
gbrlat = [flipud(eastbgrlatdec); westgbrlat; eastbgrlatdec(end)];
gbrlong = [flipud(eastbgrlongdec); westgbrlong; eastbgrlongdec(end)];