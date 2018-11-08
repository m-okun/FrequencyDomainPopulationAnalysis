function varargout = stripNans(inMat)
%strips all rows containing nans from a matrix
%jwp
% usage:
%     [v1 v2 ... vN] = stripNans( [v1, v2,...vN] )
% 
% written:
% 	jwp 11/2002
% 
% see also:
% 	stripInfs, stripZeros

if isvector(inMat) % added by Mush, 3/2017
  inMat = torow(inMat)';
end

nansPerRow = sum(isnan(inMat),2);
outMat = inMat(find(nansPerRow==0),:);

if nargout>1
	for n = 1:size(outMat,2)
        varargout{n} = outMat(:,n);
	end
else
	varargout = {outMat};
end