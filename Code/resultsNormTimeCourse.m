function [resSummaryInd] = resultsNormTimeCourse(resSummaryInd)
%RESULTSNORMTIMECOURSE normalizes single-molecule properties by values in first 3 time points of each time course
%
%SYNOPSIS [resSummaryInd] = resultsNormTimeCourse(resSummaryInd)
%
%INPUT  resSummaryInd: Output of resultsCombTimeCourseMod.
%    
%OUTPUT resSummaryInd: Same as input but with added fields numNorm0Class,
%                      probNorm0Class and densityNorm0Class.
%
%Khuloud Jaqaman, February 2017
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

%get number of movieLists
nML = max(resSummaryInd.mlIndex);

%reserve memory
[numNorm0Class,densityNorm0Class] = deal(zeros(size(resSummaryInd.numAbsClass)));
probNorm0Class = zeros(size(resSummaryInd.probAbsClass));

for iML = 1 : nML
    
    %get indices of movies in this movieList
    mdIndex = find(resSummaryInd.mlIndex == iML);
    
    %get number of time points in this movielist
    nTP = length(mdIndex);
    
    %specify indices of movies used for normalization (these are the first
    %3 movies, unless movieList has less than 3 movies)
    normIndex = mdIndex(1:min(3,nTP));
    
    %normalize numbers in various motion classes
    normConst = nanmean(resSummaryInd.numAbsClass(normIndex,:));
    numNorm0Class(mdIndex,:) = ...
        resSummaryInd.numAbsClass(mdIndex,:) ./ repmat(normConst,nTP,1);
    
    %normalize probabilities of various motion classes
    normConst = nanmean(resSummaryInd.probAbsClass(normIndex,:));
    probNorm0Class(mdIndex,:) = ...
        resSummaryInd.probAbsClass(mdIndex,:) ./ repmat(normConst,nTP,1);
    
    %normalize densities in various motion classes
    normConst = nanmean(resSummaryInd.densityAbsClass(normIndex,:));
    densityNorm0Class(mdIndex,:) = ...
        resSummaryInd.densityAbsClass(mdIndex,:) ./ repmat(normConst,nTP,1);
    
end

%store in output structure
resSummaryInd.numNorm0Class = numNorm0Class;
resSummaryInd.probNorm0Class = probNorm0Class;
resSummaryInd.densityNorm0Class = densityNorm0Class;

%% ~~~ the end ~~~
