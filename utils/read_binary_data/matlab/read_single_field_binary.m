% -
%
% SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
% SPDX-License-Identifier: MIT
%
% -
%
% setting up some parameters
%
precision = 'double';     % precision of the real-valued data
r0 = [0.,0.,0.];    % domain origin
%
% read geometry file
%
geofile  = "geometry.out";
data = dlmread(geofile);
ng   = data(1,:);
%
% read grid
%
f   = fopen('grid_x.bin');
grid_x = fread(f,[ng(1),4],precision);
fclose(f);
xp = r0(1) + grid_x(:,3)'; % centered  x grid
xu = r0(1) + grid_x(:,4)'; % staggered x grid
f   = fopen('grid_y.bin');
grid_y = fread(f,[ng(2),4],precision);
fclose(f);
yp = r0(2) + grid_y(:,3)'; % centered  y grid
yv = r0(2) + grid_y(:,4)'; % staggered y grid
f   = fopen('grid_z.bin');
grid_z = fread(f,[ng(3),4],precision);
fclose(f);
zp = r0(3) + grid_z(:,3)'; % centered  z grid
zw = r0(3) + grid_z(:,4)'; % staggered z grid
%
% read checkpoint binary file
%
filenamei = input("Name of the binary file written by CaNS (e.g. vex_fld_0000000.bin)]: ")
if isempty(filenamei)
end
iskipx      = input("Data saved every (ix, iy, iz) points. Value of ix? [1]: ")
if isempty(iskipx)
    iskipx = 1
end
iskipy      = input("Data saved every (ix, iy, iz) points. Value of iy? [1]: ")
if isempty(iskipy)
    iskipy = 1
end
iskipz      = input("Data saved every (ix, iy, iz) points. Value of iz? [1]: ")
if isempty(iskipz)
    iskipz = 1
end
iskip       = [iskipx,iskipy,iskipz]
n           = floor((ng-1)./iskip)+1
f = fopen(filenamei);
fld = fread(f,prod(n),precision);
fclose(f);
if numel(fld) ~= prod(n)
    error("expected %d values for the requested skip, found %d",prod(n),numel(fld))
end
data = reshape(fld,[n(1),n(2),n(3)]);
