function [Px] = get_projection_lowerbound(x,lb)
%GET_PROJECTION_LOWERBOUND Summary of this function goes here
%   Detailed explanation goes here

Px = x;

if ~isempty(lb)
    Px = max(lb,Px);
end

end

