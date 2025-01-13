%***********************************************************************
% Method to calculate the encoding SNR efficiency (cf. Wong 2007)
%***********************************************************************
function [SNR,meanSNR] = WongSNR(encodings)

Aplus = pinv(encodings);
N = size(Aplus,2);
SNRratio = zeros(1,size(Aplus,1));
    for n = 1:size(Aplus,1)
        AplusSumSquares = sum(Aplus(n,:).^2);
        SNRratio(n) = 1/sqrt(N*AplusSumSquares);
    end
SNR = SNRratio;

% T.O. Edit: only want the SNR efficiency of the blood components, not
% static tissue
%meanSNR = mean(SNRratio);
meanSNR = mean(SNRratio(1:end-1));

return