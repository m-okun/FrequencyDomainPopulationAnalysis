% converts input vector (of whatever orientation) to be a row vector
function x = torow(x)
if isempty(x)
  return
end
assert(isvector(x), 'input not a vector')

if size(x, 1) > 1
  x = x';
end