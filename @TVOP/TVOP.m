function  res = TVOP(Aleft)

%res = TVOP()
%
% Implements a spatial finite-differencing operator.
%
% (c) Michael Lustig 2007

if nargin >= 1
    res.Aleft = Aleft;
end

res.adjoint = 0;
res = class(res,'TVOP');



