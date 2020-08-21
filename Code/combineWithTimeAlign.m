function [valCombMSN,valCombSample] = combineWithTimeAlign(valIndSets,timeIndSets,timeComb,timeIndxComb)
%COMBINEWITHTIMEALIGN calculates average properties from multiple datasets after aligning them based on their time
%
%SYNOPSIS [valCombMSN,valCombSample] = combineWithTimeAlign(valIndSets,timeIndSets,timeComb,timeIndxComb)
%
%INPUT  valIndSets  : Cell array with parameter values from individual
%                     datasets.
%       timeIndSets : Cell array with time points of individual datasets.
%       timeComb    : Column vector of times to use for combined dataset
%                     output.
%       timeIndxComb: Index of column in timeIndSets to use for aligning
%                     the different datasets
%                     Optional. Default: 1.
%    
%OUTPUT valCombMSN: 3D array storing the mean, standard deviation and
%                   number of datapoints used to calculate each property
%                   value. 
%                   Number of rows = number of combined timecourse
%                   time points. 
%                   Number of columns = number of column in input
%                   valIndSets, which is basically the number of properties
%                   handled together.
%       valCombSample: 3D array storing the individual property values from
%                   the individual datasets used to calculate the mean and
%                   standard deviation. 
%                   Number of rows and columns = same as above.
%                   Number of entries in 3rd dimension = maximum number of
%                   datapoints contributing to any individual property.
%
%Khuloud Jaqaman, March 2015
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

%% Input

%get number of individual datasets
numIndSet = length(valIndSets);

%get number of columns (i.e. number of different properties) - use same for output
numCol = size(valIndSets{1},2);

%get number of time point in combined data
numTP = length(timeComb);

%calculate dividers between time points in order to collect relevant data points
timeCombTmp = [2*timeComb(1)-timeComb(2); timeComb; 2*timeComb(end)-timeComb(end-1)];
tpLimit = mean([timeCombTmp(1:end-1) timeCombTmp(2:end)],2);

%% Data grouping and calculation

%reserve memory
valCombSample = NaN(numTP,numCol,2*numIndSet);
valCombMSN = NaN(numTP,numCol,3);
numVal = zeros(numTP,1);

%go over each combined set data point
for iTP = 1 : numTP
    
    %collect relevant individual values
    valSample = [];
    for iSet = 1 : numIndSet
        valTmp = valIndSets{iSet}( timeIndSets{iSet}(:,timeIndxComb) >= tpLimit(iTP) & ...
            timeIndSets{iSet}(:,timeIndxComb) < tpLimit(iTP+1) ,: );
        if ~isempty(valTmp)
            valSample = [valSample; valTmp];
        end
    end
    
    %store in output variable
    numVal(iTP) = size(valSample,1);
    if numVal(iTP)~=0
        valCombSample(iTP,:,1:numVal(iTP)) = reshape(valSample',[1 numCol numVal(iTP)]);
    end
    
end

%output also mean, standard deviation and number of data points
valCombSample = valCombSample(:,:,1:max(numVal));
valCombMSN(:,:,1) = nanmean(valCombSample,3);
valCombMSN(:,:,2) = nanstd(valCombSample,[],3);
valCombMSN(:,:,3) = repmat(numVal,1,numCol);

%% ~~~ the end ~~~
