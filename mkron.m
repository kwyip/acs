function O = mkron(varargin)
% modified from Ville Bergholm 2009
O = 1;
for k = 1:nargin
  O = kron(O, varargin{k});
end
