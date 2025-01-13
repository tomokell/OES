%***********************************************************************
%Method to convert vessel locations from mm (to 4.d.p) into pixel
%locations
%***********************************************************************

function [pixves] = convert_vessels_mm2pix(vesloc,pad_size)

% Work out the FOV/2 and make the dimensions of the pixel matrix equal to the FOV
half_FOV = pad_size/2;
vessel_matrix_dim = 2*half_FOV;
%create a vector to represent the x and y axes that define the vessel
%coordinates in distance units (up to 4d.p.)
vector = -half_FOV:0.0001:half_FOV;
%find the number of distance values that can map to one pixel
values = size(vector,2)/vessel_matrix_dim;
%code to convert mm coordinates into pixel coordinates (search for values
%within 5e-5 as for large numbers there is not necessarily an exact binary
%equivalent of the number)
for m = 1:size(vesloc,2)
    xcoord = find(abs(vector-vesloc(1,m))<5e-5);
    ycoord = find(abs(vector-vesloc(2,m))<5e-5);
    for n = 1:vessel_matrix_dim
        if xcoord <= n*values && xcoord > (n-1)*values
            x(m) = n;
        end
        if ycoord <= n*values && ycoord > (n-1)*values
            y(m) = n;
        end
    end
end
pixves = [x;y]; %these are the pixel coordinates of the vessels

return