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
filenamei = input("Name of the binary restart file written by CaNS [fld.bin]: ");
if isempty(filenamei)
    filenamei = "fld.bin";
end
data       = zeros([ng(1),ng(2),ng(3),4]); % u(:,:,:),v(:,:,:),w(:,:,:),p(:,:,:)
fldinfo    = zeros([2,1]);
f = fopen(filenamei);
for p = 1:4
    data(:,:,:,p) = reshape(fread(f,ng(1)*ng(2)*ng(3),precision),[ng(1),ng(2),ng(3)]);
end
fldinfo(:) = fread(f,2,precision);
fclose(f);
%
% store data in arrays
%
u = data(:,:,:,1);
v = data(:,:,:,2);
w = data(:,:,:,3);
p = data(:,:,:,4);
time  =       fldinfo(1);
istep = round(fldinfo(2));
