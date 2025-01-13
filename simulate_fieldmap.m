%*******Code to simulate a fieldmap given linear/quadratic phase offset terms and an offset direction (x = 1 or y = 2)*******

function map = simulate_fieldmap(vesloc_FOV, linear_offset, linear_offset_dirn, quadratic_offset)

map = zeros(vesloc_FOV,vesloc_FOV);
x_vector = zeros(1,vesloc_FOV);
y_vector = zeros(vesloc_FOV,1);

%calculating an example phase offset to populate the fieldmap with
phase_vec = -pi:(2*pi/vesloc_FOV):pi;

x_vector(1,:) = phase_vec(1:vesloc_FOV);
y_vector(:,1) = phase_vec(1:vesloc_FOV);

for n = 1:vesloc_FOV
    if linear_offset_dirn == 1 && quadratic_offset == 0
        map(n,:) = x_vector;
    elseif linear_offset_dirn == 2 && quadratic_offset == 0
        map(:,n) = y_vector;
    elseif linear_offset_dirn == 0 && quadratic_offset == 1
        for m = 1:vesloc_FOV
            map(n,m) = sqrt(x_vector(m)^2 + y_vector(n)^2);
        end
    elseif linear_offset_dirn == 1 && quadratic_offset == 1
        for m = 1:vesloc_FOV
            map(n,m) = sqrt(x_vector(m)^2 + y_vector(n)^2) + x_vector(1,m);
        end
    elseif linear_offset_dirn == 2 && quadratic_offset == 1
        for m = 1:vesloc_FOV
            map(n,m) = sqrt(x_vector(m)^2 + y_vector(n)^2) + y_vector(n,1);
        end
    else
        error('Need to designate either 0, x = 1 or y = 2 for linear offset and/or quadratic offset = true')
    end
end

map = map+linear_offset;

return