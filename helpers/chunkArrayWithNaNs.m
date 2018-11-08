% Input: a - vector where some entries might be NaN
%            or a matrix where some entire(!) columns are NaN
% Ouput  b - cell array, where elements are longest contiguous (i.e. without NaN) portions of a 
function b = chunkArrayWithNaNs(a)
b = {};
if isempty(a)
  return
elseif isnan(a(1,1))
  pos = find(~isnan(a(1,:)), 1);
  b = chunkArrayWithNaNs(a(:, pos:end));
  return
end

pos = find(isnan(a(1,:)), 1);
if isempty(pos)
  b{1} = a;
  return
end

b{1} = a(:, 1:pos-1);
pos2 = find(~isnan(a(1, pos+1:end)), 1);
if ~isempty(pos2)
  b = [b; chunkArrayWithNaNs(a(:, pos+pos2:end))];
end
