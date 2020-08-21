function summarizeMergeSplitMotionVsFollowerStats(analysisRawFullFileName)
%SUMMARIZEMERGESPLITVSFOLLOWERSTATS summarizes the merge/split/motion properties of combined tracks from time course analysis with bootstrapping for statistics
%
% SYNOPSIS:
%
%       summarizeMergeSplitMotionVsFollowerStats(analysisRawFullFileName)
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
%   combinedFollowerAnalysisSummary_nameCond_timeStart-timeEnd.mat, where
%   nameCond, timeStart and timeEnd are extracted from the input file name.
%   See final line of code for stored output.
%
% Khuloud Jaqaman, July 2019
% May 2020 Updated by Zachariah Malik to treat individual movies separately
% depending on whether we bootstrapped input data.
% June 2020 Updated by Zachariah Malik; when we analyze output
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
%%ZMALIK 20200528 load nBS from analysis results as well

load(analysisRawFullFileName,'statsMergingMF','statsSplittingMF',...
    'tracks1SegMobilityMF','tracks2plusSegMobilityMF','tracksMergeSplitMF',...
    'tracks1SegFollowerPropMF','tracks2plusSegFollowerPropMF',...
    'rateStats','rateStatsF','clustStats','clustStatsF','nBS');  %#ok<NASGU>

%read condition and time points to use in saving the function's output
%also read directory where to save output (same as input)
startIndex = regexp(analysisRawFullFileName,'combinedFollowerAnalysisRaw_');
endIndex = regexp(analysisRawFullFileName,'.mat');
if isempty(endIndex)
    condName = analysisRawFullFileName(startIndex+length('combinedFollowerAnalysisRaw_'):end);
else
    condName = analysisRawFullFileName(startIndex+length('combinedFollowerAnalysisRaw_'):endIndex-1);
end
pathName = analysisRawFullFileName(1:startIndex-1);

%determine total number of samples (original + bootstrap)
nSample = length(clustStats);

% ZMALIK 20200518
% Note that if input does not have a value for nBS, then we assume
% bootstrapping case.
if exist('nBS','var')
    if nBS+1 == nSample
        BSFlag = 1; %bootstrapping
    elseif nBS == -1
        BSFlag = 0; %no bootstrapping
    else
        disp('nBS does not match nSample');
        return
    end
else %nBS not input, assume bootstrapping
    BSFlag = 1;
end
 

%% Statistics extraction and summary

if BSFlag == 1 %Bootstrapping
    %probability of merging vs follower presence before merging
    tmp = catStruct(2,'statsMergingMF.type.probVsSepSegments(:)');
    probMergeVsFollowerSep = [tmp(:,1) nanstd(tmp,[],2)];

    %probability of splitting vs follower presence before splitting
    tmp = catStruct(2,'statsSplittingMF.type.probVsTogSegment');
    probSplitVsFollowerTog = [tmp(:,1) nanstd(tmp,[],2)];

    %probability of follower presence after merging vs. before merging
    tmp = catStruct(2,'statsMergingMF.type.probFollowerTogVsSep(:)');
    probFollowerAftMergeVsFollowerBef = [tmp(:,1) nanstd(tmp,[],2)];

    %probability of follower presence after splitting vs. before splitting
    tmp = catStruct(2,'statsSplittingMF.type.probFollowerSepVsTog(:)');
    probFollowerAftSplitVsFollowerBef = [tmp(:,1) nanstd(tmp,[],2)];

    %cluster fractions
    %per follower presence
    tmp1 = clustStatsF(1).clusterFrac;
    numClustGlobal = size(tmp1,1); %number of cluster sizes, used in subsequent steps as well
    numMode = 2;
    tmp = [mean(tmp1,2) zeros(numClustGlobal,nSample-1,numMode)];
    for iSamp = 2 : nSample
    tmp1 = clustStatsF(iSamp).clusterFrac;
    numClust = size(tmp1,1);
    tmp(1:numClust,iSamp,:) = mean(tmp1,2);
    end
    clustFracFollower = [tmp(:,1,:) nanstd(tmp,[],2)];
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
    %per follower presence
    tmp = NaN(numClustGlobal,nSample,numMode);
    for iSamp = 1 : nSample
    tmp1 = rateStatsF(iSamp).rateOffPerClust;
    numClust = size(tmp1,1);
    tmp(1:numClust,iSamp,:) = tmp1;
    end
    rateOffPerClustFollower = [tmp(:,1,:) nanstd(tmp,[],2)];

    %movement vs. follower for tracks with 1 segment

    %fraction of tracks with follower
    %overall
    tmp = horzcat(tracks1SegMobilityMF.fracTracksWithF);
    fracTracksWithFollowerOverall1Seg = [tmp(1) nanstd(tmp(2:end))];
    %per motion mode
    tmp = horzcat(tracks1SegMobilityMF.fracTracksWithFperMode);
    fracTracksWithFollowerPerMode1Seg = [tmp(:,1) nanstd(tmp,[],2)];

    %fraction of diffusion modes among tracks
    %without follower
    tmp = horzcat(tracks1SegMobilityMF.fracDiffModeTracksNoF);
    fracDiffModeTracksNoF1Seg = [tmp(:,1) nanstd(tmp,[],2)];
    %with follower - during 
    tmp = cat(3,tracks1SegMobilityMF.fracDiffModeComboTracksWithF);
    tmp1 = squeeze(sum(tmp,1));
    fracDiffModeTracksWithFDur1Seg = [tmp1(:,1) nanstd(tmp1,[],2)];
    %with follower - before
    tmp1 = squeeze(sum(tmp,2));
    fracDiffModeTracksWithFBef1Seg = [tmp1(:,1) nanstd(tmp1,[],2)];

    %For each mode before follower, fraction of tracks that swtich to slowest
    %mode during follower
    tmp1 = squeeze(tmp(:,end,:)) ./ tmp1;
    fracMode1DurPerModeBefF1Seg = [tmp1(:,1) nanstd(tmp1,[],2)];

    %diffusion coefficient
    %without follower
    tmp = horzcat(tracks1SegMobilityMF.meanDiffCoefModeNoF);
    diffCoefModeNoF1Seg = [tmp(:,1) nanstd(tmp,[],2)];
    %with follower - before
    tmp = cat(3,tracks1SegMobilityMF.meanDiffCoefBeforeModeCombo);
    for i=1:size(tmp,3)
    tmp1(:,i) = diag(tmp(:,:,i));
    end
    diffCoefModeWithFBef1Seg = [tmp1(:,1) nanstd(tmp1,[],2)];
    %with follower - during 
    tmp = cat(3,tracks1SegMobilityMF.meanDiffCoefDuringModeCombo);
    for i=1:size(tmp,3)
    tmp1(:,i) = diag(tmp(:,:,i));
    end
    diffCoefModeWithFDur1Seg = [tmp1(:,1) nanstd(tmp1,[],2)];

    %movement vs. follower for tracks with 2+ segment

    %fraction of tracks with follower
    %overall
    tmp = horzcat(tracks2plusSegMobilityMF.fracTracksWithF);
    fracTracksWithFollowerOverall2plusSeg = [tmp(1) nanstd(tmp(2:end))];
    %per motion mode
    tmp = horzcat(tracks2plusSegMobilityMF.fracTracksWithFperMode);
    fracTracksWithFollowerPerMode2plusSeg = [tmp(:,1) nanstd(tmp,[],2)];

    %fraction of diffusion modes among tracks
    %without follower
    tmp = horzcat(tracks2plusSegMobilityMF.fracDiffModeTracksNoF);
    fracDiffModeTracksNoF2plusSeg = [tmp(:,1) nanstd(tmp,[],2)];
    %with follower - during 
    tmp = cat(3,tracks2plusSegMobilityMF.fracDiffModeComboTracksWithF);
    tmp1 = squeeze(sum(tmp,1));
    fracDiffModeTracksWithFDur2plusSeg = [tmp1(:,1) nanstd(tmp1,[],2)];
    %with follower - before
    tmp1 = squeeze(sum(tmp,2));
    fracDiffModeTracksWithFBef2plusSeg = [tmp1(:,1) nanstd(tmp1,[],2)];

    %For each mode before follower, fraction of tracks that swtich to slowest
    %mode during follower
    tmp1 = squeeze(tmp(:,end,:)) ./ tmp1;
    fracMode1DurPerModeBefF2plusSeg = [tmp1(:,1) nanstd(tmp1,[],2)];

    %diffusion coefficient
    %without follower
    tmp = horzcat(tracks2plusSegMobilityMF.meanDiffCoefModeNoF);
    diffCoefModeNoF2plusSeg = [tmp(:,1) nanstd(tmp,[],2)];
    %with follower - before
    tmp = cat(3,tracks2plusSegMobilityMF.meanDiffCoefBeforeModeCombo);
    for i=1:size(tmp,3)
    tmp1(:,i) = diag(tmp(:,:,i));
    end
    diffCoefModeWithFBef2plusSeg = [tmp1(:,1) nanstd(tmp1,[],2)];
    %with follower - during
    tmp = cat(3,tracks2plusSegMobilityMF.meanDiffCoefDuringModeCombo);
    for i=1:size(tmp,3)
    tmp1(:,i) = diag(tmp(:,:,i));
    end
    diffCoefModeWithFDur2plusSeg = [tmp1(:,1) nanstd(tmp1,[],2)];

    %merging vs follower
    %all of the below are for merging events with follower after merge

    %fraction with follower before merge
    tmp = horzcat(tracksMergeSplitMF.fracFollowerBefAmongMergesWithFollowerAft);
    fracFBefAmongMergewithFAft = [tmp(1) nanstd(tmp(2:end))];

    %for those with follower before merge
    %fraction of monomers
    tmp = vertcat(tracksMergeSplitMF.fracMonoWithFollowerBefMerge)';
    fracMonoWithFBefAmongMergeWithFAft = [tmp(:,1) nanstd(tmp(:,2:end),[],2)];
    %mean oligomeric state
    tmp = vertcat(tracksMergeSplitMF.meanOligoWithFollowerBefMerge)';
    meanOligoWithFBefAmongMergeWithFAft = [tmp(:,1) nanstd(tmp(:,2:end),[],2)];

    %for those without follower before merge
    %fraction of monomers
    tmp = vertcat(tracksMergeSplitMF.fracMonoNoFollowerBefMerge)';
    fracMonoNoFBefAmongMergeWithFAft = [tmp(:,1) nanstd(tmp(:,2:end),[],2)];
    %mean oligomeric state
    tmp = vertcat(tracksMergeSplitMF.meanOligoNoFollowerBefMerge)';
    meanOligoNoFBefAmongMergeWithFAft = [tmp(:,1) nanstd(tmp(:,2:end),[],2)];

    %follower properties vs. oligomeric state for tracks with 1 segment and
    %with follower

    %fraction of tracks that have follower from their start
    tmp = horzcat(tracks1SegFollowerPropMF.fracTime0Start2Follower);
    fracTime0Start2Follower1Seg = [tmp(:,1) nanstd(tmp,[],2)];

    %mean time from track start to follower appearance
    tmp = horzcat(tracks1SegFollowerPropMF.meanTimeStart2Follower);
    meanTimeStart2Follower1Seg = [tmp(:,1) nanstd(tmp,[],2)];

    %fraction of time track has follower
    tmp = horzcat(tracks1SegFollowerPropMF.meanFracTimeWithFollower);
    fracTimeWithFollower1Seg = [tmp(:,1) nanstd(tmp,[],2)];

    %follower consistency (i.e. fraction of time it is present between its
    %appearance and disappearance)
    tmp = horzcat(tracks1SegFollowerPropMF.meanFollowerConsistency);
    followerConsistency1Seg = [tmp(:,1) nanstd(tmp,[],2)];

    %mean track duration
    tmp = horzcat(tracks1SegFollowerPropMF.meanTrackDuration);
    meanTrackDuration1Seg = [tmp(:,1) nanstd(tmp,[],2)];

    %follower properties vs. oligomeric state for tracks with 2+ segments and
    %with follower

    %fraction of tracks that have follower from their start
    tmp = horzcat(tracks2plusSegFollowerPropMF.fracTime0Start2Follower);
    fracTime0Start2Follower2plusSeg = [tmp(:,1) nanstd(tmp,[],2)];

    %mean time from track start to follower appearance
    tmp = horzcat(tracks2plusSegFollowerPropMF.meanTimeStart2Follower);
    meanTimeStart2Follower2plusSeg = [tmp(:,1) nanstd(tmp,[],2)];

    %fraction of time track has follower
    tmp = horzcat(tracks2plusSegFollowerPropMF.meanFracTimeWithFollower);
    fracTimeWithFollower2plusSeg = [tmp(:,1) nanstd(tmp,[],2)];

    %follower consistency (i.e. fraction of time it is present between its
    %appearance and disappearance)
    tmp = horzcat(tracks2plusSegFollowerPropMF.meanFollowerConsistency);
    followerConsistency2plusSeg = [tmp(:,1) nanstd(tmp,[],2)];

    %mean track duration
    tmp = horzcat(tracks2plusSegFollowerPropMF.meanTrackDuration);
    meanTrackDuration2plusSeg = [tmp(:,1) nanstd(tmp,[],2)];

else %Treat movies/simulations separately
    %probability of merging vs follower presence before merging
    tmp = catStruct(2,'statsMergingMF.type.probVsSepSegments(:)');
    probMergeVsFollowerSep = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(probMergeVsFollowerSep(1,:))
        probMergeVsFollowerSep(rowInd,2) = ...
            probMergeVsFollowerSep(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %probability of splitting vs follower presence before splitting
    tmp = catStruct(2,'statsSplittingMF.type.probVsTogSegment');
    probSplitVsFollowerTog = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(probSplitVsFollowerTog(:,1))
        probSplitVsFollowerTog(rowInd,2) = ...
            probSplitVsFollowerTog(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %probability of follower presence after merging vs. before merging
    tmp = catStruct(2,'statsMergingMF.type.probFollowerTogVsSep(:)');
    probFollowerAftMergeVsFollowerBef = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(probFollowerAftMergeVsFollowerBef(:,1))
        probFollowerAftMergeVsFollowerBef(rowInd,2) = ...
            probFollowerAftMergeVsFollowerBef(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %probability of follower presence after splitting vs. before splitting
    tmp = catStruct(2,'statsSplittingMF.type.probFollowerSepVsTog(:)');
    probFollowerAftSplitVsFollowerBef = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(probFollowerAftSplitVsFollowerBef(:,1))
        probFollowerAftSplitVsFollowerBef(rowInd,2) = ...
            probFollowerAftSplitVsFollowerBef(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %cluster fractions
    %per follower presence
    tmp1 = clustStatsF(1).clusterFrac;
    numClustGlobal = size(tmp1,1); %number of cluster sizes, used in subsequent steps as well
    numMode = 2;
    tmp = [mean(tmp1,2) zeros(numClustGlobal,nSample-1,numMode)];
    for iSamp = 2 : nSample
    tmp1 = clustStatsF(iSamp).clusterFrac;
    numClust = size(tmp1,1);
    tmp(1:numClust,iSamp,:) = mean(tmp1,2);
    end
    clustFracFollower = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(clustFracFollower(:,1))
        clustFracFollower(rowInd,2) = ...
            clustFracFollower(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %overall
    tmp = zeros(numClustGlobal,nSample);
    for iSamp = 1 : nSample
    tmp1 = clustStats(iSamp).clusterFrac;
    numClust = size(tmp1,1);
    tmp(1:numClust,iSamp) = mean(tmp1,2);
    end
    clustFracOverall = [nanmean(tmp,2) nanstd(tmp,[],2)];

    for rowInd = 1:length(clustFracOverall(:,1))
        clustFracOverall(rowInd,2) = ...
            clustFracOverall(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    
    %off rate
    %overall
    tmp = NaN(numClustGlobal,nSample);
    for iSamp = 1 : nSample
    tmp1 = rateStats(iSamp).rateOffPerClust;
    numClust = size(tmp1,1);
    tmp(1:numClust,iSamp) = tmp1;
    end
    rateOffPerClustOverall = [nanmean(tmp,2) nanstd(tmp,[],2)];
    
    for rowInd = 1:length(rateOffPerClustOverall(:,1))
        rateOffPerClustOverall(rowInd,2) = ...
            rateOffPerClustOverall(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %per follower presence
    tmp = NaN(numClustGlobal,nSample,numMode);
    for iSamp = 1 : nSample
    tmp1 = rateStatsF(iSamp).rateOffPerClust;
    numClust = size(tmp1,1);
    tmp(1:numClust,iSamp,:) = tmp1;
    end
    rateOffPerClustFollower = [nanmean(tmp,2) nanstd(tmp,[],2)];

    for rowInd = 1:length(rateOffPerClustFollower(:,1))
        rateOffPerClustFollower(rowInd,2) = ...
            rateOffPerClustFollower(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    
    %movement vs. follower for tracks with 1 segment

    %fraction of tracks with follower
    %overall
    tmp = horzcat(tracks1SegMobilityMF.fracTracksWithF);
    fracTracksWithFollowerOverall1Seg = [nanmean(tmp) nanstd(tmp)];
    for rowInd = 1:length(fracTracksWithFollowerOverall1Seg(:,1))
        fracTracksWithFollowerOverall1Seg(rowInd,2) = ...
            fracTracksWithFollowerOverall1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %per motion mode
    tmp = horzcat(tracks1SegMobilityMF.fracTracksWithFperMode);
    fracTracksWithFollowerPerMode1Seg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(fracTracksWithFollowerPerMode1Seg(:,1))
        fracTracksWithFollowerPerMode1Seg(rowInd,2) = ...
            fracTracksWithFollowerPerMode1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %fraction of diffusion modes among tracks
    %without follower
    tmp = horzcat(tracks1SegMobilityMF.fracDiffModeTracksNoF);
    fracDiffModeTracksNoF1Seg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(fracDiffModeTracksNoF1Seg(:,1))
        fracDiffModeTracksNoF1Seg(rowInd,2) = ...
            fracDiffModeTracksNoF1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %with follower - during 
    tmp = cat(3,tracks1SegMobilityMF.fracDiffModeComboTracksWithF);
    tmp1 = squeeze(sum(tmp,1));
    fracDiffModeTracksWithFDur1Seg = [nanmean(tmp1,2) nanstd(tmp1,[],2)];
    for rowInd = 1:length(fracDiffModeTracksWithFDur1Seg(:,1))
        fracDiffModeTracksWithFDur1Seg(rowInd,2) = ...
            fracDiffModeTracksWithFDur1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp1(rowInd,:))));
    end
    %with follower - before
    tmp1 = squeeze(sum(tmp,2));
    fracDiffModeTracksWithFBef1Seg = [nanmean(tmp1,2) nanstd(tmp1,[],2)];
    for rowInd = 1:length(fracDiffModeTracksWithFBef1Seg(:,1))
         fracDiffModeTracksWithFBef1Seg(rowInd,2) = ...
             fracDiffModeTracksWithFBef1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp1(rowInd,:))));
    end
    %For each mode before follower, fraction of tracks that swtich to slowest
    %mode during follower
    tmp1 = squeeze(tmp(:,end,:)) ./ tmp1;
    fracMode1DurPerModeBefF1Seg = [nanmean(tmp1,2) nanstd(tmp1,[],2)];
    for rowInd = 1:length(fracMode1DurPerModeBefF1Seg(:,1))
         fracMode1DurPerModeBefF1Seg(rowInd,2) = ...
             fracMode1DurPerModeBefF1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp1(rowInd,:))));
    end
    %diffusion coefficient
    %without follower
    tmp = horzcat(tracks1SegMobilityMF.meanDiffCoefModeNoF);
    diffCoefModeNoF1Seg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(diffCoefModeNoF1Seg(:,1))
         diffCoefModeNoF1Seg(rowInd,2) = ...
             diffCoefModeNoF1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %with follower - before
    tmp = cat(3,tracks1SegMobilityMF.meanDiffCoefBeforeModeCombo);
    for i=1:size(tmp,3)
    tmp1(:,i) = diag(tmp(:,:,i));
    end
    diffCoefModeWithFBef1Seg = [nanmean(tmp1,2) nanstd(tmp1,[],2)];
    for rowInd = 1:length(diffCoefModeWithFBef1Seg(:,1))
         diffCoefModeWithFBef1Seg(rowInd,2) = ...
             diffCoefModeWithFBef1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp1(rowInd,:))));
    end
    %with follower - during 
    tmp = cat(3,tracks1SegMobilityMF.meanDiffCoefDuringModeCombo);
    for i=1:size(tmp,3)
    tmp1(:,i) = diag(tmp(:,:,i));
    end
    diffCoefModeWithFDur1Seg = [nanmean(tmp1,2) nanstd(tmp1,[],2)];
    for rowInd = 1:length(diffCoefModeWithFDur1Seg(:,1))
         diffCoefModeWithFDur1Seg(rowInd,2) = ...
             diffCoefModeWithFDur1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp1(rowInd,:))));
    end
    %movement vs. follower for tracks with 2+ segment

    %fraction of tracks with follower
    %overall
    tmp = horzcat(tracks2plusSegMobilityMF.fracTracksWithF);
    fracTracksWithFollowerOverall2plusSeg = [nanmean(tmp) nanstd(tmp)];
    for rowInd = 1:length(fracTracksWithFollowerOverall2plusSeg(:,1))
         fracTracksWithFollowerOverall2plusSeg(rowInd,2) = ...
             fracTracksWithFollowerOverall2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %per motion mode
    tmp = horzcat(tracks2plusSegMobilityMF.fracTracksWithFperMode);
    fracTracksWithFollowerPerMode2plusSeg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(fracTracksWithFollowerPerMode2plusSeg(:,1))
         fracTracksWithFollowerPerMode2plusSeg(rowInd,2) = ...
             fracTracksWithFollowerPerMode2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %fraction of diffusion modes among tracks
    %without follower
    tmp = horzcat(tracks2plusSegMobilityMF.fracDiffModeTracksNoF);
    fracDiffModeTracksNoF2plusSeg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(fracDiffModeTracksNoF2plusSeg(:,1))
         fracDiffModeTracksNoF2plusSeg(rowInd,2) = ...
             fracDiffModeTracksNoF2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %with follower - during 
    tmp = cat(3,tracks2plusSegMobilityMF.fracDiffModeComboTracksWithF);
    tmp1 = squeeze(sum(tmp,1));
    fracDiffModeTracksWithFDur2plusSeg = [nanmean(tmp1,2) nanstd(tmp1,[],2)];
    for rowInd = 1:length(fracDiffModeTracksWithFDur2plusSeg(:,1))
         fracDiffModeTracksWithFDur2plusSeg(rowInd,2) = ...
             fracDiffModeTracksWithFDur2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp1(rowInd,:))));
    end
    %with follower - before
    tmp1 = squeeze(sum(tmp,2));
    fracDiffModeTracksWithFBef2plusSeg = [nanmean(tmp1,2) nanstd(tmp1,[],2)];
    for rowInd = 1:length(fracDiffModeTracksWithFBef2plusSeg(:,1))
         fracDiffModeTracksWithFBef2plusSeg(rowInd,2) = ...
             fracDiffModeTracksWithFBef2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp1(rowInd,:))));
    end
    %For each mode before follower, fraction of tracks that swtich to slowest
    %mode during follower
    tmp1 = squeeze(tmp(:,end,:)) ./ tmp1;
    fracMode1DurPerModeBefF2plusSeg = [nanmean(tmp1,2) nanstd(tmp1,[],2)];
    for rowInd = 1:length(fracMode1DurPerModeBefF2plusSeg(:,1))
         fracMode1DurPerModeBefF2plusSeg(rowInd,2) = ...
             fracMode1DurPerModeBefF2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp1(rowInd,:))));
    end
    %diffusion coefficient
    %without follower
    tmp = horzcat(tracks2plusSegMobilityMF.meanDiffCoefModeNoF);
    diffCoefModeNoF2plusSeg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(diffCoefModeNoF2plusSeg(:,1))
         diffCoefModeNoF2plusSeg(rowInd,2) = ...
             diffCoefModeNoF2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %with follower - before
    tmp = cat(3,tracks2plusSegMobilityMF.meanDiffCoefBeforeModeCombo);
    for i=1:size(tmp,3)
    tmp1(:,i) = diag(tmp(:,:,i));
    end
    diffCoefModeWithFBef2plusSeg = [nanmean(tmp1,2) nanstd(tmp1,[],2)];
    for rowInd = 1:length(diffCoefModeWithFBef2plusSeg(:,1))
         diffCoefModeWithFBef2plusSeg(rowInd,2) = ...
             diffCoefModeWithFBef2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp1(rowInd,:))));
    end
    %with follower - during
    tmp = cat(3,tracks2plusSegMobilityMF.meanDiffCoefDuringModeCombo);
    for i=1:size(tmp,3)
    tmp1(:,i) = diag(tmp(:,:,i));
    end
    diffCoefModeWithFDur2plusSeg = [nanmean(tmp1,2) nanstd(tmp1,[],2)];
    for rowInd = 1:length(diffCoefModeWithFBef2plusSeg(:,1))
         diffCoefModeWithFDur2plusSeg(rowInd,2) = ...
             diffCoefModeWithFDur2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp1(rowInd,:))));
    end
    %merging vs follower
    %all of the below are for merging events with follower after merge

    %fraction with follower before merge
    tmp = horzcat(tracksMergeSplitMF.fracFollowerBefAmongMergesWithFollowerAft);
    fracFBefAmongMergewithFAft = [nanmean(tmp) nanstd(tmp)];
    for rowInd = 1:length(fracFBefAmongMergewithFAft(:,1))
         fracFBefAmongMergewithFAft(rowInd,2) = ...
             fracFBefAmongMergewithFAft(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %for those with follower before merge
    %fraction of monomers
    tmp = vertcat(tracksMergeSplitMF.fracMonoWithFollowerBefMerge)';
    fracMonoWithFBefAmongMergeWithFAft = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(fracMonoWithFBefAmongMergeWithFAft(:,1))
         fracMonoWithFBefAmongMergeWithFAft(rowInd,2) = ...
             fracMonoWithFBefAmongMergeWithFAft(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %mean oligomeric state
    tmp = vertcat(tracksMergeSplitMF.meanOligoWithFollowerBefMerge)';
    meanOligoWithFBefAmongMergeWithFAft = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(meanOligoWithFBefAmongMergeWithFAft(:,1))
         meanOligoWithFBefAmongMergeWithFAft(rowInd,2) = ...
             meanOligoWithFBefAmongMergeWithFAft(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %for those without follower before merge
    %fraction of monomers
    tmp = vertcat(tracksMergeSplitMF.fracMonoNoFollowerBefMerge)';
    fracMonoNoFBefAmongMergeWithFAft = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(fracMonoNoFBefAmongMergeWithFAft(:,1))
         fracMonoNoFBefAmongMergeWithFAft(rowInd,2) = ...
             fracMonoNoFBefAmongMergeWithFAft(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %mean oligomeric state
    tmp = vertcat(tracksMergeSplitMF.meanOligoNoFollowerBefMerge)';
    meanOligoNoFBefAmongMergeWithFAft = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(meanOligoNoFBefAmongMergeWithFAft(:,1))
         meanOligoNoFBefAmongMergeWithFAft(rowInd,2) = ...
             meanOligoNoFBefAmongMergeWithFAft(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end

    %follower properties vs. oligomeric state for tracks with 1 segment and
    %with follower

    %fraction of tracks that have follower from their start
    tmp = horzcat(tracks1SegFollowerPropMF.fracTime0Start2Follower);
    fracTime0Start2Follower1Seg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(fracTime0Start2Follower1Seg(:,1))
         fracTime0Start2Follower1Seg(rowInd,2) = ...
             fracTime0Start2Follower1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %mean time from track start to follower appearance
    tmp = horzcat(tracks1SegFollowerPropMF.meanTimeStart2Follower);
    meanTimeStart2Follower1Seg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(meanTimeStart2Follower1Seg(:,1))
         meanTimeStart2Follower1Seg(rowInd,2) = ...
             meanTimeStart2Follower1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %fraction of time track has follower
    tmp = horzcat(tracks1SegFollowerPropMF.meanFracTimeWithFollower);
    fracTimeWithFollower1Seg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(fracTimeWithFollower1Seg(:,1))
         fracTimeWithFollower1Seg(rowInd,2) = ...
             fracTimeWithFollower1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %follower consistency (i.e. fraction of time it is present between its
    %appearance and disappearance)
    tmp = horzcat(tracks1SegFollowerPropMF.meanFollowerConsistency);
    followerConsistency1Seg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(followerConsistency1Seg(:,1))
         followerConsistency1Seg(rowInd,2) = ...
             followerConsistency1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %mean track duration
    tmp = horzcat(tracks1SegFollowerPropMF.meanTrackDuration);
    meanTrackDuration1Seg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(meanTrackDuration1Seg(:,1))
         meanTrackDuration1Seg(rowInd,2) = ...
             meanTrackDuration1Seg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %follower properties vs. oligomeric state for tracks with 2+ segments and
    %with follower

    %fraction of tracks that have follower from their start
    tmp = horzcat(tracks2plusSegFollowerPropMF.fracTime0Start2Follower);
    fracTime0Start2Follower2plusSeg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(fracTime0Start2Follower2plusSeg(:,1))
         fracTime0Start2Follower2plusSeg(rowInd,2) = ...
             fracTime0Start2Follower2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %mean time from track start to follower appearance
    tmp = horzcat(tracks2plusSegFollowerPropMF.meanTimeStart2Follower);
    meanTimeStart2Follower2plusSeg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(meanTimeStart2Follower2plusSeg(:,1))
         meanTimeStart2Follower2plusSeg(rowInd,2) = ...
             meanTimeStart2Follower2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %fraction of time track has follower
    tmp = horzcat(tracks2plusSegFollowerPropMF.meanFracTimeWithFollower);
    fracTimeWithFollower2plusSeg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(fracTimeWithFollower2plusSeg(:,1))
         fracTimeWithFollower2plusSeg(rowInd,2) = ...
             fracTimeWithFollower2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %follower consistency (i.e. fraction of time it is present between its
    %appearance and disappearance)
    tmp = horzcat(tracks2plusSegFollowerPropMF.meanFollowerConsistency);
    followerConsistency2plusSeg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(followerConsistency2plusSeg(:,1))
         followerConsistency2plusSeg(rowInd,2) = ...
             followerConsistency2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
    %mean track duration
    tmp = horzcat(tracks2plusSegFollowerPropMF.meanTrackDuration);
    meanTrackDuration2plusSeg = [nanmean(tmp,2) nanstd(tmp,[],2)];
    for rowInd = 1:length(meanTrackDuration2plusSeg(:,1))
         meanTrackDuration2plusSeg(rowInd,2) = ...
             meanTrackDuration2plusSeg(rowInd,2)/sqrt(nnz(~isnan(tmp(rowInd,:))));
    end
end
%% Save all results

save([pathName 'combinedFollowerAnalysisSummary_' condName '.mat'],...
    'probMergeVsFollowerSep','probSplitVsFollowerTog',...
    'probFollowerAftMergeVsFollowerBef','probFollowerAftSplitVsFollowerBef',...
    'clustFracFollower','clustFracOverall','rateOffPerClustOverall','rateOffPerClustFollower',...
    'fracTracksWithFollowerOverall1Seg','fracTracksWithFollowerPerMode1Seg',...
    'fracDiffModeTracksNoF1Seg','fracDiffModeTracksWithFDur1Seg',...
    'fracDiffModeTracksWithFBef1Seg','fracMode1DurPerModeBefF1Seg',...
    'diffCoefModeNoF1Seg','diffCoefModeWithFBef1Seg','diffCoefModeWithFDur1Seg',...
    'fracTracksWithFollowerOverall2plusSeg','fracTracksWithFollowerPerMode2plusSeg',...
    'fracDiffModeTracksNoF2plusSeg','fracDiffModeTracksWithFDur2plusSeg',...
    'fracDiffModeTracksWithFBef2plusSeg','fracMode1DurPerModeBefF2plusSeg',...
    'diffCoefModeNoF2plusSeg','diffCoefModeWithFBef2plusSeg',...
    'diffCoefModeWithFDur2plusSeg','fracFBefAmongMergewithFAft',...
    'fracMonoWithFBefAmongMergeWithFAft','meanOligoWithFBefAmongMergeWithFAft',...
    'fracMonoNoFBefAmongMergeWithFAft','meanOligoNoFBefAmongMergeWithFAft',...
    'fracTime0Start2Follower1Seg','meanTimeStart2Follower1Seg','fracTimeWithFollower1Seg',...
    'followerConsistency1Seg','meanTrackDuration1Seg',...
    'fracTime0Start2Follower2plusSeg','meanTimeStart2Follower2plusSeg','fracTimeWithFollower2plusSeg',...
    'followerConsistency2plusSeg','meanTrackDuration2plusSeg');

%% ~~~ the end ~~~
