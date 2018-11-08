% Input: A - matrix of phases.
% 1st dim - different variables (e.g. units)
% 2nd dim - frequencies
% Will try to unwrap form both ends, and return the one that is better
% Furthermore will process each interval separately (interval is a part of the sequence not broken by NaNs)
function B = bestUnwrap(A, opt)

if nargin < 2 || ~isfield(opt, 'smooth')
  opt.smooth = false;
end

if isempty(A) || all(isnan(A(:)))
  B = [];
  return
elseif any(isnan(A(:)))
  B = A;
  for k = 1:size(A,1)
    I = find(isnan(A(k,:)));
    if isempty(I) % no NaNs in this row
      B(k,:) = bestUnwrap(A(k,:), opt);
      continue
    end
    pos = 1;
    for i = 1:numel(I)
      B(k, pos:I(i)-1) = bestUnwrap(A(k, pos:I(i)-1), opt);
      posPrev = find(~isnan(B(k,pos-1:-1:1)), 1); % position of the closest previous phases which isn't NaN
      if isempty(posPrev)
        pos = I(i)+1;
        continue
      end
      
      posPrev = pos - posPrev;
      % if shifting up/down by 2*pi brings us closer to the end of the previous chunk, do it
      if abs(diff([B(k, pos) + 2*pi, B(k, posPrev)])) < abs(diff([B(k, pos), B(k, posPrev)]))
        B(k, pos:I(i)-1) = B(k, pos:I(i)-1) + 2*pi;
      elseif abs(diff([B(k, pos) - 2*pi, B(k, posPrev)])) < abs(diff([B(k, pos), B(k, posPrev)]))
        B(k, pos:I(i)-1) = B(k, pos:I(i)-1) - 2*pi;
      end
      pos = I(i)+1;
    end % loop on segments without NaNs
    
    % Finally we need to handle the piece between last NaN and the end of the row
    B(k, I(end)+1:end) = bestUnwrap(A(k, I(end)+1:end), opt);
    if I(end)+1 > size(B, 2)
      continue
    end
    pos = I(end)+1;
    posPrev = find(~isnan(B(k,pos-1:-1:1)), 1);
    posPrev = pos - posPrev;
    if abs(diff([B(k, pos) + 2*pi, B(k, posPrev)])) < abs(diff([B(k, pos), B(k, posPrev)]))
      B(k, pos:end) = B(k, pos:end) + 2*pi;
    elseif abs(diff([B(k, pos) - 2*pi, B(k, posPrev)])) < abs(diff([B(k, pos), B(k, posPrev)]))
      B(k, pos:end) = B(k, pos:end) - 2*pi;
    end    
  end % loop on different rows 
  B = B - round(median(stripNans(B(:)))/(2*pi))*2*pi;
  return
end

% If we're here, there's no NaN business...
b1 = unwrap(A, [], 2);
b2 = zeros(size(b1));
b2(:,end:-1:1) = unwrap(A(:,end:-1:1), [], 2);
e1 = 0; e2 = 0;
for i = 1:size(A,1)
  for j = 1:i-1
    e1 = e1 + sum((b1(i,:) - b1(j,:)).^2);
    e2 = e2 + sum((b2(i,:) - b2(j,:)).^2);
  end
end
if e1 < e2
  B = b1;
else
  B = b2;
end
B = B - round(median(B(:))/(2*pi))*2*pi;
if opt.smooth
  B = smooth(B, 7);
end