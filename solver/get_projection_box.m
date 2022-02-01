function [Px] = get_projection_box(x,lb,ub)
%GET_PROJECTION_BOX Summary of this function goes here
%   Detailed explanation goes here

Px = x;

if ~isempty(lb)
    Px = max(lb,Px);
end

if ~isempty(ub)
    Px = min(ub,Px);
end

end

