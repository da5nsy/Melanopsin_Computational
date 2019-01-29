function [vrhel_sq, S_RF_f] = vrhel_square(refs)

load sur_vrhel.mat sur_vrhel S_vrhel
if nargin == 0
    refs = [38, 15, 134, 137, 138, 65, 19, 24, 140, 26];
end
T_SRF = sur_vrhel(:,refs);

vrhel_sq = corr(T_SRF');

S_RF_f = SToWls(S_vrhel);

end