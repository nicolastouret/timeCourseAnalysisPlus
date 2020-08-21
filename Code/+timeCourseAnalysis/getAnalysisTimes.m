function [ analysisTimes, timeLimit, timeLimitIndx ] = getAnalysisTimes( times, timeResolution, inOutFlag )
%getAnalysisTimes Get analysis times
%
% INPUT
% times - from commonInfo
% timeResolution - from params
%
% OUTPUT
% analysisTimes - for commonInfo
% timeLimit - 
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

times = cellfun(@(x,y) x(y==1), times, inOutFlag,'UniformOutput',false);
timeMax = cellfun(@(x) max(x), times);
timeMin = cellfun(@(x) min(x), times);
%determine the overall range of all conditions (This is useful later when comparing two curves)
timeMaxMax = max(timeMax);
timeMinMin = min(timeMin);
%convert to number divisible by timeResolution
timeMinMin = timeMinMin - mod(timeMinMin, timeResolution);
timeMaxMax = timeMaxMax - mod(timeMaxMax, -timeResolution);
analysisTime_Union = timeMinMin : timeResolution : timeMaxMax;
%Determine the indx and value of timemin and max for each conditions
nConditions = numel(times);
timeLimit = cell(size(times));
timeLimitIndx = cell(size(times));
for iCond = 1:nConditions
    timeLimitIndx{iCond} = [find(analysisTime_Union <= timeMin(iCond), 1, 'last'), find(analysisTime_Union >= timeMax(iCond), 1, 'first')];
    timeLimit{iCond} = [analysisTime_Union(timeLimitIndx{iCond}(1)), analysisTime_Union(timeLimitIndx{iCond}(2))];
end
% timeLimit = reshape(timeLimit,size(data));
%store this in commonInfo
analysisTimes = cellfun(@(x) x(1) : timeResolution : x(2), timeLimit, 'UniformOutput', false);

end

