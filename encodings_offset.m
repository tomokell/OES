% Calculate the encoding matrix using a realistic modulation function
% for a given set of vessel locations, encodings and off-resonance values

function [vessel_offset_enc] = encodings_offset(vesloc, frequency, phase, off_res)

load modmat_lookup

pp = zeros(size(vesloc,2));
%find the encoding value at the vessel locations
vessel_offset_enc = zeros(1,size(vesloc,2));

for d = 1:size(vesloc,2)
    pp(d) = frequency*vesloc(1:2,d) + phase + off_res(d);
    pp(d) = mod(pp(d)+pi,2*pi)-pi; % Wrap to -pi -> +pi
    vessel_offset_enc(d) = interp1(modmat_phase, modmat_uni_lookup,pp(d));
end

return
