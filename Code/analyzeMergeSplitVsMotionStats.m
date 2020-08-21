function analyzeMergeSplitVsMotionStats(prepTracksFullFileName,nBS)
%ANALYZEMERGESPLITVSMOTIONSTATS analyzes the merge/split/motion properties of combined tracks from time course analysis with bootstrapping for statistics
%
% SYNOPSIS:
%
%       analyzeMergeSplitVsMotionStats(prepTracksFullFileName,nBS)
%
% INPUT:
%
%   prepTracksFullFileName: Full name (including path) of .mat file saving
%   the tracks to be analyzed, as saved by prepareTracksForMSAnalysis.
%
%   nBS: Number of boostrap samples. If 0, no bootstrapping.
%
% OUTPUT:
%
%   No direct output, but analysis results are  stored in a .mat file in
%   the same directory as the input file. File name would be
%   combinedMSAnalysisRaw_nameCond_timeStart-timeEnd.mat, where nameCond,
%   timeStart and timeEnd are extracted from the input file name. See final
%   line of code for stored output.
%
% Khuloud Jaqaman, April 2019
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

%% Initalization and setup

%get prepared tracks etc.
load(prepTracksFullFileName,'tracksAltCell','tracksDefCell','diffTypeCell',...
    'diffModeCell','clustHistCell','cellMaskArea');

%read condition and time points to use for saving the function's output
%also read directory where to save output (same as input)
startIndex = regexp(prepTracksFullFileName,'tracksEtcForMSAnalysis_');
endIndex = regexp(prepTracksFullFileName,'.mat');
if isempty(endIndex)
    condName = prepTracksFullFileName(startIndex+length('tracksEtcForMSAnalysis_'):end);
else
    condName = prepTracksFullFileName(startIndex+length('tracksEtcForMSAnalysis_'):endIndex-1);
end
pathName = prepTracksFullFileName(1:startIndex-1);

%Initialize various variables
numInd = size(tracksAltCell,1);
[bsIndxSamp,tracksAltSamp,tracksDefSamp,diffTypeSamp,diffModeSamp,clustHistSamp] ...
    = deal(cell(numInd,1));

%get bootstrap samples for each movie in ensemble
for k = 1 : numInd
    
    %Get tracks
    numTracks = length(tracksAltCell{k});
    
    %determine bootstrap samples
    %first sample is original sample
    randSampleIndx = NaN(numTracks,nBS+1);
    randSampleIndx(:,1) = (1 : numTracks)';
    for iBS = 1 : nBS
        randSampleIndx(:,iBS+1) = (randsample(numTracks,numTracks,true))';
    end
    bsIndxSamp{k} = randSampleIndx;
    
end %(for k = 1 : numInd)

%% Analysis and bootstrap

% How many diffusion modes are being analyzed?
diffModeTotal = vertcat(diffModeCell{:});
dClass = max(vertcat(diffModeTotal.diffMode));

%go over original sample (1) and bootstrap samples (2 -> nBS+1)
for iBS = nBS+1 : -1 : 1
    
    %get current sample
    for k = 1 : numInd
        tracksAltSamp{k} = tracksAltCell{k}(bsIndxSamp{k}(:,iBS));
        tracksDefSamp{k} = tracksDefCell{k}(bsIndxSamp{k}(:,iBS));
        diffTypeSamp{k}  = diffTypeCell{k}(bsIndxSamp{k}(:,iBS));
        diffModeSamp{k}  = diffModeCell{k}(bsIndxSamp{k}(:,iBS));
        clustHistSamp{k} = clustHistCell{k}(bsIndxSamp{k}(:,iBS));
    end
    
    %put information together, forgetting about individual movies
    %         tracksAltTotal = vertcat(tracksAltSamp{:});
    tracksDefTotal = vertcat(tracksDefSamp{:});
    %         diffTypeTotal = vertcat(diffTypeSamp{:});
    %         diffModeTotal = vertcat(diffModeSamp{:});
    clustHistTotal = vertcat(clustHistSamp{:});
    clustHistTotalMerge = vertcat(clustHistTotal{:});
    
    %probability of diffusion modes
    %(total comes later from calcStatsMSAlt)
    mergeEventsAll = clustHistTotalMerge(clustHistTotalMerge(:,6)==2,:);
    mergeThenSplitEvents = mergeEventsAll(mergeEventsAll(:,7)==1,:);
    tmp = hist(mergeEventsAll(:,10),1:dClass);
    probModeAfterMerge(:,iBS) = tmp(end:-1:1)'/sum(tmp); %go backwards to be consistent with output of calcStatsMSAlt
    tmp = hist(mergeThenSplitEvents(:,10),1:dClass);
    probModeAfterMergeBefSplit(:,iBS) = tmp(end:-1:1)'/sum(tmp); %go backwards to be consistent with output of calcStatsMSAlt
    
    %probability that a split event is followed by a merge event
    %probability that the merge is between the same two molecules that split
    %(overall probability of splitting comes later from calcStatsMSAlt)
    msTimeInfo = calcMergeSplitTimesNew(tracksDefTotal,5,[],0);
    numSplit = length(msTimeInfo.all.timeStart2Split) + length(msTimeInfo.all.timeMerge2Split);
    numSplit2MergeSelf = length(msTimeInfo.all.timeSplit2MergeSelf);
    numSplit2MergeOther = length(msTimeInfo.all.timeSplit2MergeOther);
    numSplit2MergeTot = numSplit2MergeSelf + numSplit2MergeOther;
    probSplit2Merge(1,iBS) = numSplit2MergeTot / numSplit;
    probSplit2MergeSelf(1,iBS) = numSplit2MergeSelf / numSplit2MergeTot;
    
    %lots of merging/splitting vs. motion stats - ignore unclassified
    [statsGeneralNoUn(iBS),statsMotionNoUn(iBS),statsMergingNoUn(iBS),statsSplittingNoUn(iBS)] = ...
        calcStatsMSAlt(tracksAltSamp,diffTypeSamp,diffModeSamp,5,1);
    
    %lots of merging/splitting vs. motion stats - keep unclassified
    [statsGeneralAll(iBS),statsMotionAll(iBS),statsMergingAll(iBS),statsSplittingAll(iBS)] = ...
        calcStatsMSAlt(tracksAltSamp,diffTypeSamp,diffModeSamp,5,0);
    
    %labeled receptor on and off rates and oligomer densities
    infoSpaceTime = struct('probDim',2,'timeStep',1,'area',cellMaskArea);
    [rateOnPerClust,rateOffPerClust,densityPerClust,numClustForRateCalc,clustStats(iBS),...
        rateOnPerClustM,rateOffPerClustM,densityPerClustM,numClustForRateCalcM,clustStatsM(iBS)] = ...
        clusterOnOffRatesAndDensityFromClustHistory(clustHistTotalMerge,infoSpaceTime);
    rateStats(iBS).rateOnPerClust = rateOnPerClust;
    rateStats(iBS).rateOffPerClust = rateOffPerClust;
    rateStats(iBS).densityPerClust = densityPerClust;
    rateStats(iBS).numClustForRateCalc = numClustForRateCalc;
    rateStatsM(iBS).rateOnPerClust = rateOnPerClustM;
    rateStatsM(iBS).rateOffPerClust = rateOffPerClustM;
    rateStatsM(iBS).densityPerClust = densityPerClustM;
    rateStatsM(iBS).numClustForRateCalc = numClustForRateCalcM;
    
end %(for iBS = 1 : nBS+1)

%% Save all results

save([pathName 'combinedMSAnalysisRaw_' condName '.mat'],'bsIndxSamp',...
    'statsGeneralNoUn','statsMotionNoUn','statsMergingNoUn','statsSplittingNoUn',...
    'statsGeneralAll','statsMotionAll','statsMergingAll','statsSplittingAll',...
    'probModeAfterMerge','probModeAfterMergeBefSplit',...
    'probSplit2Merge','probSplit2MergeSelf',...
    'rateStats','rateStatsM','clustStats','clustStatsM');

%% ~~~ the end ~~~
