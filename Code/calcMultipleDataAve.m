function [dataAve, inOutFlag] = calcMultipleDataAve(data, times, inOutFlag, aveInterval, shiftTime, ignoreIsolatedPts)
%Scatter plots given data sets in one plot and averages in aveInterval intervals
%
%SYNOPSIS [dataAve] = calcMultipleDataAve(data, times, inOutFlag, aveInterval)
%
%INPUT
%   data            : cell array of data sets (each set must be a column)
%   times           : cell array of times corresponding to the dataset (each set must be a column)
%   plotTitle       : name of the plot
%   inOutFlag       : cell array of flags indicating inlier data points
%                     (value 1) and outlier data points (value 0)
%   aveInterval     : Interval for averaging. Default: 3 min
%   shiftTime       : Value by which time is shifted. Needed to shift time
%                     back temporarily for this analysis.
%   ignoreIsolatedPts: 1 to ignore isolated points when averaging (and then
%                     later when doign the spline fit; 0 otherwise.
%                     Default: 1.
%
%OUTPUT
%   aveData         :
%
% Khuloud Jaqaman, March 2017
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

%% Initialization
%assign default value
if nargin<4 || isempty(aveInterval)
    aveInterval = 3;
end
if nargin<5 || isempty(shiftTime)
    shiftTime = zeros(size(data));
end
if nargin<6 || isempty(ignoreIsolatedPts)
    ignoreIsolatedPts = 1;
end

shiftTimeCell = cell(size(data));
for iCell = 1 : length(shiftTimeCell)
    shiftTimeCell{iCell} = shiftTime(iCell);
end
if max(shiftTime>0)
    shiftTimeAll = mean(shiftTime(shiftTime~=0));
else
    shiftTimeAll = 0;
end

%creates figure and stores the figure handle
fxn = @(varargin) calcMultipleDataAvePerCondition(varargin{:}, aveInterval, ignoreIsolatedPts, shiftTimeAll);
[dataAve, inOutFlag] = cellfun(fxn, data, times, inOutFlag, shiftTimeCell, 'UniformOutput',false);

end


%% sub-function

function [dataAve, inOutFlag] = calcMultipleDataAvePerCondition(data, times, inOutFlag, shiftTime, aveInterval, ignoreIsolatedPts, shiftTimeAll)

try
    %     %remove nan and inf
    %     mask = isfinite(data);
    %     data = data(mask);
    %     times = times(mask);
    %     inOutFlag = inOutFlag(mask);
    %
    %     %KJ: discard outliers
    %     timesIn = times(inOutFlag);
    %     dataIn = data(inOutFlag);
    
    %divide time points into bins and take time and data average and
    %std
    binID = floor((times-shiftTimeAll)/aveInterval);
    binUnique = unique(binID);
    nBin = length(binUnique);
    [nDataPoints,dataInAve,dataInStd,timeInAve,timeInStd] = deal(NaN(nBin,1));
    for iBin = 1 : nBin
        indx = find(binID == binUnique(iBin));
        finiteFlagBin = isfinite(data(indx));
        inOutFlagBin = inOutFlag(indx);
        indxGood = indx(finiteFlagBin&inOutFlagBin==1);
        nDP = length(indxGood);
        nDataPoints(iBin) = nDP;
        if ignoreIsolatedPts && nDP < 5
            inOutFlag(indx) = -1;
        else
            dataInAve(iBin) = mean(data(indxGood));
            timeInAve(iBin) = mean(times(indxGood));
            dataInStd(iBin) = std(data(indxGood));
            timeInStd(iBin) = std(times(indxGood));
        end
    end
    
    %output
    timeBinEdges = [binUnique binUnique+1]*aveInterval + shiftTimeAll;
    dataAve = struct('timeBinEdges',timeBinEdges,'timeAve',timeInAve,'dataAve',dataInAve,...
        'timeStd',timeInStd,'dataStd',dataInStd,'nDataPts',nDataPoints);
    
catch
    dataAve = [];
end

end