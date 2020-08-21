function [ varargout ] = pararrayfun_progress( func, varargin )
%pararrayfun_progress Run arrayfun in parallel with progress meter
%
% Copyright (C) 2020, Jaqaman Lab - UT Southwestern 
%
% This file is part of timeCourseAnalysisPlus.
% 
% timeCourseAnalysisPlus is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% timeCourseAnalysisPlus is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with timeCourseAnalysisPlus.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Find where the parameters begin
arraySize = size(varargin{1});
for paramIdx=1:length(varargin)
    inSize = size(varargin{paramIdx});
    if(~all(inSize == arraySize))
        break;
    end
end
if(mod(length(varargin)-paramIdx,2)==0 ...
    && paramIdx > 1 ...
    && ischar(varargin{paramIdx-1}))
    % Remaining arguments for parameters is an odd number.
    % Previous argument could be a parameter name that
    %    happens to be the same size as other array input
    paramIdx = paramIdx - 1;
elseif(paramIdx == length(varargin))
    paramIdx = paramIdx + 1;
else
end

% Just convert arrays to cells for now and use parcellfun_progress
arrays = cellfun(@num2cell,varargin(1:paramIdx-1),'UniformOutput',false);

[varargout{1:nargout}] = parcellfun_progress(func, arrays{:}, varargin{paramIdx:end});

end
