function varargout = expandSizes(ignoredim, varargin)
nin=nargin-1;
siz = zeros(nin,1); 
% siz = cellfun(@size, varargin, 'uniformoutput', 0); 
% msiz = max(cellfun(@length, siz)); 

for ii = 1:nin
    si = size(varargin{ii}); 
    siz(ii, 1:length(si)) = si;
%     siz{ii}(end+1:msiz) = 1;
end
siz(siz==0)=1; 
% siz = vertcat(siz{:}); 

sizMax = max(siz, [],1); 
% varargout = cell(1,nin); 
for ii = 1:nin
    rep = sizMax.*(sizMax > siz(ii,:));
    rep(ignoredim) = 0; 
%     assert(all(siz(ii, rep>0) == 1), 'Sizes disagree on non-singleton dimensions!');
    if any(rep)
        rep(rep == 0) = 1; 
        varargout{ii} = repmat(varargin{ii}, rep); 
    else
        varargout{ii} = varargin{ii}; 
    end
end