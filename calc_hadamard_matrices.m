% Calculate and return a hadamard matrix for encoding via the OES

function [hadamard_complete, hadamard_for_OES, TagMode] = calc_hadamard_matrices(vesloc,had_enc,standard_enc,PCASL_enc,off_res_flag,off_res)

if had_enc == true
    [hadamard_for_OES, TagMode] = choose_hadamard(vesloc);
elseif standard_enc == true
    if size(vesloc,2) ~= 4
        error('standard_enc only works for 4 vessels')
    end
    hadamard_for_OES = [-1 -1 -1 -1; ...
                         1  1  1  1; ...
                        -1  1 -1  1; ...
                         1 -1  1 -1; ...
                        -1 -1  1  1; ...
                         1  1 -1 -1; ...
                        -1  1  1 -1; ...
                         1 -1 -1  1];
    TagMode = [0 1 2 2 2 2 2 2];
elseif PCASL_enc == true
    control_all = ones(1,size(vesloc,2));
    tag_all = -control_all;
    hadamard_for_OES = [tag_all; control_all];
    TagMode = [0 1];
else
    error('Desired type of ideal encodings must be chosen as true/false in the inputs to the function')
end

%creating the full Hadamard/encoding matrix (with static tissue etc) for
%display/comparison
static_tissue = ones(size(hadamard_for_OES,1),1);
hadamard_complete = [hadamard_for_OES,static_tissue];

%If there are off resonance values then adjust the ideal encoding scheme.
%NB: offres values need to be in radians...
if off_res_flag
    %turn the ideal scheme in to 0s and pis (pi phase at control location; the real phase we want for an ideal control)
    h_phi = (1+hadamard_for_OES)*(pi/2); 
    h_phi_offres = zeros(size(h_phi));
    %remove the off resonance phase that will accrue between pulses due to the field inhomogeneity
    for n = 1:size(h_phi,1)
        h_phi_offres(n,:) = h_phi(n,:) - off_res; 
    end
    %Turn the encodings into complex numbers (e^itheta)
    h_offres = exp(1i*h_phi_offres);
    hadamard_for_OES = h_offres;
end