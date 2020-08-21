function summarizeMergeSplitVsMotionStats(analysisRawFullFileName)
%SUMMARIZEMERGESPLITVSMOTIONSTATS summarizes the merge/split/motion properties of combined tracks from time course analysis with bootstrapping for statistics
%
% SYNOPSIS:
%
%       summarizeMergeSplitVsMotionStats(analysisRawFullFileName)
%
% INPUT:
%
%   analysisRawFullFileName: Full name (including path) of .mat file saving
%   the raw analysis, as saved by analyzeMergeSplitVsMotionStats.
%
% OUTPUT:
%
%   No direct output, but analysis summary is stored in a .mat file in
%   the same directory as the input file. File name would be
%   combinedMSAnalysisSummary_nameCond_timeStart-timeEnd.mat, where nameCond,
%   timeStard and timeEnd are extracted from the input file name. See final
%   line of code for stored output.
%
% Khuloud Jaqaman, May 2019
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

%get analysis results
load(analysisRawFullFileName,'statsMotionNoUn',...
    'statsMergingNoUn','statsSplittingNoUn',...
    'statsMergingAll','statsSplittingAll',...
    'probModeAfterMerge','probModeAfterMergeBefSplit',...
    'probSplit2Merge','probSplit2MergeSelf',...
    'rateStats','rateStatsM','clustStats','clustStatsM'); %#ok<NASGU>

%read condition and time points to use in saving the function's output
%also read directory where to save output (same as input)
startIndex = regexp(analysisRawFullFileName,'combinedMSAnalysisRaw_');
endIndex = regexp(analysisRawFullFileName,'.mat');
if isempty(endIndex)
    condName = analysisRawFullFileName(startIndex+length('combinedMSAnalysisRaw_'):end);
else
    condName = analysisRawFullFileName(startIndex+length('combinedMSAnalysisRaw_'):endIndex-1);
end
pathName = analysisRawFullFileName(1:startIndex-1);

%determine total number of samples (original + bootstrap)
nSample = length(clustStats);

%% Statistics extraction and summary

%probabilities related to splitting and then merging
tmp = horzcat(statsSplittingAll.probOverall);
probSplitOverall = [tmp(1) nanstd(tmp)];
probSplit2Merge = [probSplit2Merge(1) nanstd(probSplit2Merge)];
probSplit2MergeSelf = [probSplit2MergeSelf(1) nanstd(probSplit2MergeSelf)];

%probabilties related to motion modes overall and upon merging
tmp = catStruct(2,'statsMotionNoUn.mode.prob');
probModeOverall = [tmp(:,1) nanstd(tmp,[],2)];
probModeAfterMerge = [probModeAfterMerge(:,1) nanstd(probModeAfterMerge,[],2)];
probModeAfterMergeBefSplit = [probModeAfterMergeBefSplit(:,1) nanstd(probModeAfterMergeBefSplit,[],2)];

%overall probability of merging
tmp = horzcat(statsMergingAll.probOverall);
probMergeOverall = [tmp(1) nanstd(tmp)];

%probability of merging vs motion mode before merging
tmp = catStruct(2,'statsMergingAll.mode.probVsSepSegments(:)');
probMergeVsModeSep = [tmp(:,1) nanstd(tmp,[],2)];

%probability of splitting vs motion mode before splitting
tmp = catStruct(2,'statsSplittingAll.mode.probVsTogSegment');
probSplitVsModeTog = [tmp(:,1) nanstd(tmp,[],2)];

%probability of motion after merging vs. before merging
tmp = catStruct(2,'statsMergingNoUn.mode.probMotionTogVsSep(:)');
probModeAftMergeVsModeBef = [tmp(:,1) nanstd(tmp,[],2)];

%probability of motion after splitting vs. before splitting
tmp = catStruct(2,'statsSplittingNoUn.mode.probMotionSepVsTog(:)');
probModeAftSplitVsModeBef = [tmp(:,1) nanstd(tmp,[],2)];

%mean diffusion coefficient before/after merging/splitting
tmp = catStruct(2,'statsMergingNoUn.mode.meanDiffCoefAftBef12(:,4)');
diffCoefMeanAftMerge = [tmp(:,1) nanstd(tmp,[],2)];
tmp = catStruct(2,'statsMergingNoUn.mode.meanDiffCoefAftBef12(:,5)');
diffCoefMeanBefMergeH = [tmp(:,1) nanstd(tmp,[],2)];
tmp = catStruct(2,'statsMergingNoUn.mode.meanDiffCoefAftBef12(:,6)');
diffCoefMeanBefMergeL = [tmp(:,1) nanstd(tmp,[],2)];
tmp = catStruct(2,'statsSplittingNoUn.mode.meanDiffCoefBefAft12(:,4)');
diffCoefMeanBefSplit = [tmp(:,1) nanstd(tmp,[],2)];
tmp = catStruct(2,'statsSplittingNoUn.mode.meanDiffCoefBefAft12(:,5)');
diffCoefMeanAftSplitH = [tmp(:,1) nanstd(tmp,[],2)];
tmp = catStruct(2,'statsSplittingNoUn.mode.meanDiffCoefBefAft12(:,6)');
diffCoefMeanAftSplitL = [tmp(:,1) nanstd(tmp,[],2)];

%mean ratio of diffusion coefficients before/after merging/splitting
tmp = catStruct(2,'statsMergingNoUn.mode.meanDiffCoefAftBef12(:,7)');
diffCoefRatioMergeH = [tmp(:,1) nanstd(tmp,[],2)];
tmp = catStruct(2,'statsMergingNoUn.mode.meanDiffCoefAftBef12(:,8)');
diffCoefRatioMergeL = [tmp(:,1) nanstd(tmp,[],2)];
tmp = catStruct(2,'statsSplittingNoUn.mode.meanDiffCoefBefAft12(:,7)');
diffCoefRatioSplitH = [tmp(:,1) nanstd(tmp,[],2)];
tmp = catStruct(2,'statsSplittingNoUn.mode.meanDiffCoefBefAft12(:,8)');
diffCoefRatioSplitL = [tmp(:,1) nanstd(tmp,[],2)];

%cluster fractions
%per motion mode
tmp1 = clustStatsM(1).clusterFrac;
[numClustGlobal,~,numMode] = size(tmp1); %number of cluster and motion modes, used in subsequent steps as well
tmp = [mean(tmp1,2) zeros(numClustGlobal,nSample-1,numMode)];
for iSamp = 2 : nSample
    tmp1 = clustStatsM(iSamp).clusterFrac;
    numClust = size(tmp1,1);
    tmp(1:numClust,iSamp,:) = mean(tmp1,2);
end
clustFracMotion = [tmp(:,1,:) nanstd(tmp,[],2)];
%overall
tmp = zeros(numClustGlobal,nSample);
for iSamp = 1 : nSample
    tmp1 = clustStats(iSamp).clusterFrac;
    numClust = size(tmp1,1);
    tmp(1:numClust,iSamp) = mean(tmp1,2);
end
clustFracOverall = [tmp(:,1) nanstd(tmp,[],2)];

%off rate
%overall
tmp = NaN(numClustGlobal,nSample);
for iSamp = 1 : nSample
    tmp1 = rateStats(iSamp).rateOffPerClust;
    numClust = size(tmp1,1);
    tmp(1:numClust,iSamp) = tmp1;
end
rateOffPerClustOverall = [tmp(:,1) nanstd(tmp,[],2)];
%per motion mode
tmp = NaN(numClustGlobal,nSample,numMode);
for iSamp = 1 : nSample
    tmp1 = rateStatsM(iSamp).rateOffPerClust;
    numClust = size(tmp1,1);
    tmp(1:numClust,iSamp,:) = tmp1;
end
rateOffPerClustMotion = [tmp(:,1,:) nanstd(tmp,[],2)];

%on rate, overall (no per motion mode)
tmp = NaN(numClustGlobal,nSample);
for iSamp = 1 : nSample
    tmp1 = rateStats(iSamp).rateOnPerClust;
    numClust = size(tmp1,1);
    tmp(1:numClust,iSamp) = tmp1;
end
rateOnPerClustOverall = [tmp(:,1) nanstd(tmp,[],2)];

%% Save all results

save([pathName 'combinedMSAnalysisSummary_' condName '.mat'],...
    'probSplitOverall','probSplit2Merge','probSplit2MergeSelf',...
    'probModeOverall','probModeAfterMerge','probModeAfterMergeBefSplit',...
    'probMergeOverall','probMergeVsModeSep','probSplitVsModeTog',...
    'probModeAftMergeVsModeBef','probModeAftSplitVsModeBef',...
    'diffCoefMeanAftMerge','diffCoefMeanBefMergeH','diffCoefMeanBefMergeL',...
    'diffCoefMeanBefSplit','diffCoefMeanAftSplitH','diffCoefMeanAftSplitL',...
    'diffCoefRatioMergeH','diffCoefRatioMergeL','diffCoefRatioSplitH','diffCoefRatioSplitL',...
    'clustFracOverall','clustFracMotion','rateOnPerClustOverall',...
    'rateOffPerClustOverall','rateOffPerClustMotion');

%% ~~~ the end ~~~
