function analyzeMergeSplitMotionVsFollowerStats(prepTracksFullFileName,diffModeDividerStruct,nBS)
%ANALYZEMERGESPLITMOTIONVSFOLLOWERSTATS analyzes the merge/split/motion/follower properties of combined tracks from time course analysis with bootstrapping for statistics
%
% SYNOPSIS:
%
%       analyzeMergeSplitMotionVsFollower(prepTracksFullFileName,diffModeDividerStruct,nBS)
%
% INPUT:
%
%   prepTracksFullFileName: Full name (including path) of .mat file saving
%   the tracks to be analyzed, as saved by prepareTracksForMSAnalysis.
%
%   diffModeDividerStruct: Diffusion mode divider structure, as input to
%   trackDiffModeAnalysis.
%
%   nBS: Number of boostrap samples. If 0, no bootstrapping.
%       If -1, No longer bootstrapping, but we use each simulation/movie
%       independently.
%
% OUTPUT:
%
%   No direct output, but analysis results are  stored in a .mat file in
%   the same directory as the input file. File name would be
%   combinedFollowerAnalysisRaw_nameCond_timeStart-timeEnd.mat, where nameCond,
%   timeStart and timeEnd are extracted from the input file name. See final
%   line of code for stored output.
%
% Khuloud Jaqaman, July 2019
% May 2020 Updated by Zachariah Malik to treat individual movies separately
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
    'diffModeCell','clustHistCell','cellMaskArea','tracksFollowerCell','fracFollowerCell');

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

if nBS ~= -1
    %Initialize variables if bootstrapping
    [bsIndxSamp,tracksAltTrial,tracksDefTrial,diffTypeTrial,diffModeTrial,...
    clustHistTrial,tracksFollowerTrial,fracFollowerTrial] ...
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
end %(if nBS ~= -1)

%% Analysis and bootstrap

%go over original sample (1) and bootstrap samples (2 -> nBS+1)
if nBS ~= -1
    for iBS = nBS+1 : -1 : 1

        %get current sample
        for k = 1 : numInd
            tracksAltTrial{k} = tracksAltCell{k}(bsIndxSamp{k}(:,iBS));
            tracksDefTrial{k} = tracksDefCell{k}(bsIndxSamp{k}(:,iBS));
            diffTypeTrial{k}  = diffTypeCell{k}(bsIndxSamp{k}(:,iBS));
            diffModeTrial{k}  = diffModeCell{k}(bsIndxSamp{k}(:,iBS));
            clustHistTrial{k} = clustHistCell{k}(bsIndxSamp{k}(:,iBS));
            tracksFollowerTrial{k} = tracksFollowerCell{k}(bsIndxSamp{k}(:,iBS));
            fracFollowerTrial{k} = fracFollowerCell{k}(bsIndxSamp{k}(:,iBS));
        end

        %put information together, forgetting about individual movies
        tracksAltTotal = vertcat(tracksAltTrial{:});
        %         tracksDefTotal = vertcat(tracksDefSamp{:});
        %         diffTypeTotal = vertcat(diffTypeSamp{:});
        %         diffModeTotal = vertcat(diffModeSamp{:});
        clustHistTotal = vertcat(clustHistTrial{:});
        clustHistTotalMerge = vertcat(clustHistTotal{:});
        tracksFollowerTotal = vertcat(tracksFollowerTrial{:});
        fracFollowerTotal = vertcat(fracFollowerTrial{:});
        fracFollowerTotalMerge = vertcat(fracFollowerTotal{:});

        %merging, splitting and motion analysis relative to follower presence
        [tracks1SegMobilityMF(iBS),tracks2plusSegMobilityMF(iBS),tracksMergeSplitMF(iBS),...
            tracks1SegFollowerPropMF(iBS),tracks2plusSegFollowerPropMF(iBS)] = ...
            analyzeMasterMoveMergeSplitRelFollower(tracksAltTotal,tracksFollowerTotal,...
            diffModeDividerStruct,10);

        %more merging/splitting stats vs. follower presence
        [statsGeneralMF(iBS),statsFollowerMF(iBS),statsMergingMF(iBS),statsSplittingMF(iBS)] = ...
            calcStatsMSperF(tracksAltTrial,5);

        %labeled receptor on and off rates and oligomer densities
        infoSpaceTime = struct('probDim',2,'timeStep',1,'area',cellMaskArea);
        [rateOnPerClust,rateOffPerClust,densityPerClust,numClustForRateCalc,clustStats(iBS),...
            rateOnPerClustM,rateOffPerClustM,densityPerClustM,numClustForRateCalcM,clustStatsM(iBS), ...
            rateOnPerClustF,rateOffPerClustF,densityPerClustF,numClustForRateCalcF,clustStatsF(iBS)] = ...
            clusterOnOffRatesAndDensityFromClustHistory(clustHistTotalMerge,infoSpaceTime,fracFollowerTotalMerge);

        rateStats(iBS).rateOnPerClust = rateOnPerClust;
        rateStats(iBS).rateOffPerClust = rateOffPerClust;
        rateStats(iBS).densityPerClust = densityPerClust;
        rateStats(iBS).numClustForRateCalc = numClustForRateCalc;

        rateStatsM(iBS).rateOnPerClust = rateOnPerClustM;
        rateStatsM(iBS).rateOffPerClust = rateOffPerClustM;
        rateStatsM(iBS).densityPerClust = densityPerClustM;
        rateStatsM(iBS).numClustForRateCalc = numClustForRateCalcM;

        rateStatsF(iBS).rateOnPerClust = rateOnPerClustF;
        rateStatsF(iBS).rateOffPerClust = rateOffPerClustF;
        rateStatsF(iBS).densityPerClust = densityPerClustF;
        rateStatsF(iBS).numClustForRateCalc = numClustForRateCalcF;

    end %(for iBS = 1 : nBS+1)

    %take care of size inconsistencies
    %tracks with 1 segment
    numRowGlob = length(tracks1SegFollowerPropMF(1).fracTime0Start2Follower);
    for iBS = 1 : nBS+1
        numRowLoc = length(tracks1SegFollowerPropMF(iBS).fracTime0Start2Follower);
        if numRowLoc >= numRowGlob
            tracks1SegFollowerPropMF(iBS).fracTime0Start2Follower = ...
                tracks1SegFollowerPropMF(iBS).fracTime0Start2Follower(1:numRowGlob);
            tracks1SegFollowerPropMF(iBS).meanTimeStart2Follower = ...
                tracks1SegFollowerPropMF(iBS).meanTimeStart2Follower(1:numRowGlob);
            tracks1SegFollowerPropMF(iBS).meanFracTimeWithFollower = ...
                tracks1SegFollowerPropMF(iBS).meanFracTimeWithFollower(1:numRowGlob);
            tracks1SegFollowerPropMF(iBS).meanFollowerConsistency = ...
                tracks1SegFollowerPropMF(iBS).meanFollowerConsistency(1:numRowGlob);
            tracks1SegFollowerPropMF(iBS).meanTrackDuration = ...
                tracks1SegFollowerPropMF(iBS).meanTrackDuration(1:numRowGlob);
        else
            tracks1SegFollowerPropMF(iBS).fracTime0Start2Follower = ...
                [tracks1SegFollowerPropMF(iBS).fracTime0Start2Follower; NaN(numRowGlob-numRowLoc,1)];
            tracks1SegFollowerPropMF(iBS).meanTimeStart2Follower = ...
                [tracks1SegFollowerPropMF(iBS).meanTimeStart2Follower; NaN(numRowGlob-numRowLoc,1)];
            tracks1SegFollowerPropMF(iBS).meanFracTimeWithFollower = ...
                [tracks1SegFollowerPropMF(iBS).meanFracTimeWithFollower; NaN(numRowGlob-numRowLoc,1)];
            tracks1SegFollowerPropMF(iBS).meanFollowerConsistency = ...
                [tracks1SegFollowerPropMF(iBS).meanFollowerConsistency; NaN(numRowGlob-numRowLoc,1)];
            tracks1SegFollowerPropMF(iBS).meanTrackDuration = ...
                [tracks1SegFollowerPropMF(iBS).meanTrackDuration; NaN(numRowGlob-numRowLoc,1)];
        end
    end
    %tracks with 2+ segments
    numRowGlob = length(tracks2plusSegFollowerPropMF(1).fracTime0Start2Follower);
    for iBS = 1 : nBS+1
        numRowLoc = length(tracks2plusSegFollowerPropMF(iBS).fracTime0Start2Follower);
        if numRowLoc >= numRowGlob
            tracks2plusSegFollowerPropMF(iBS).fracTime0Start2Follower = ...
                tracks2plusSegFollowerPropMF(iBS).fracTime0Start2Follower(1:numRowGlob);
            tracks2plusSegFollowerPropMF(iBS).meanTimeStart2Follower = ...
                tracks2plusSegFollowerPropMF(iBS).meanTimeStart2Follower(1:numRowGlob);
            tracks2plusSegFollowerPropMF(iBS).meanFracTimeWithFollower = ...
                tracks2plusSegFollowerPropMF(iBS).meanFracTimeWithFollower(1:numRowGlob);
            tracks2plusSegFollowerPropMF(iBS).meanFollowerConsistency = ...
                tracks2plusSegFollowerPropMF(iBS).meanFollowerConsistency(1:numRowGlob);
            tracks2plusSegFollowerPropMF(iBS).meanTrackDuration = ...
                tracks2plusSegFollowerPropMF(iBS).meanTrackDuration(1:numRowGlob);
        else
            tracks2plusSegFollowerPropMF(iBS).fracTime0Start2Follower = ...
                [tracks2plusSegFollowerPropMF(iBS).fracTime0Start2Follower; NaN(numRowGlob-numRowLoc,1)];
            tracks2plusSegFollowerPropMF(iBS).meanTimeStart2Follower = ...
                [tracks2plusSegFollowerPropMF(iBS).meanTimeStart2Follower; NaN(numRowGlob-numRowLoc,1)];
            tracks2plusSegFollowerPropMF(iBS).meanFracTimeWithFollower = ...
                [tracks2plusSegFollowerPropMF(iBS).meanFracTimeWithFollower; NaN(numRowGlob-numRowLoc,1)];
            tracks2plusSegFollowerPropMF(iBS).meanFollowerConsistency = ...
                [tracks2plusSegFollowerPropMF(iBS).meanFollowerConsistency; NaN(numRowGlob-numRowLoc,1)];
            tracks2plusSegFollowerPropMF(iBS).meanTrackDuration = ...
                [tracks2plusSegFollowerPropMF(iBS).meanTrackDuration; NaN(numRowGlob-numRowLoc,1)];
        end
    end
    %%%%ZMALIK 20200604 Add these variables to analyze individual movies
    %with follower - during 
    tmp = cat(3,tracks1SegMobilityMF.fracDiffModeComboTracksWithF);
    tmp1 = squeeze(sum(tmp,1));
    fracDiffModeTracksWithFDur1Seg = tmp1;
    %with follower - before
    tmp1 = squeeze(sum(tmp,2));
    fracDiffModeTracksWithFBef1Seg = tmp1;

    %For each mode before follower, fraction of tracks that swtich to slowest
    %mode during follower
    tmp1 = squeeze(tmp(:,end,:)) ./ tmp1;
    fracMode1DurPerModeBefF1Seg = tmp1;
    
    for indStruct = 1:length(tracks1SegMobilityMF)
        tracks1SegMobilityMF(indStruct).fracDiffModeTracksWithFDur1Seg = ...
            fracDiffModeTracksWithFDur1Seg(:,indStruct);
        tracks1SegMobilityMF(indStruct).fracDiffModeTracksWithFBef1Seg = ...
            fracDiffModeTracksWithFBef1Seg(:,indStruct);
        tracks1SegMobilityMF(indStruct).fracMode1DurPerModeBefF1Seg = ...
            fracMode1DurPerModeBefF1Seg(:,indStruct);
    end
else %No bootstrapping and treat each movie individually
%%% 2020/May/18 ZMalik Now we treat each movie separately
%Assume minTrackLen = 5
minTrackLen = 5;
    for trialNum = numInd : -1 : 1
        %get current sample
        tracksAltTrial{1} = tracksAltCell{trialNum}(:,1);
        tracksDefTrial{1} = tracksDefCell{trialNum}(:,1);
        diffTypeTrial{1}  = diffTypeCell{trialNum}(:,1);
        diffModeTrial{1}  = diffModeCell{trialNum}(:,1);
        clustHistTrial{1} = clustHistCell{trialNum}(:,1);
        tracksFollowerTrial{1} = tracksFollowerCell{trialNum}(:,1);
        fracFollowerTrial{1} = fracFollowerCell{trialNum}(:,1);


        %put information together
        tracksAltTotal = vertcat(tracksAltTrial{:});
        %         tracksDefTotal = vertcat(tracksDefSamp{:});
        %         diffTypeTotal = vertcat(diffTypeSamp{:});
        %         diffModeTotal = vertcat(diffModeSamp{:});
        clustHistTotal = vertcat(clustHistTrial{:});
        clustHistTotalMerge = vertcat(clustHistTotal{:});
        tracksFollowerTotal = vertcat(tracksFollowerTrial{:});
        fracFollowerTotal = vertcat(fracFollowerTrial{:});
        fracFollowerTotalMerge = vertcat(fracFollowerTotal{:});

        %keep only tracks with length >= minTrackLen
        criteria.lifeTime.min = minTrackLen;
        indx = chooseTracks(tracksAltTotal,criteria);
        clear criteria
        if isempty(indx)
            continue %skip this loop iteration
        end
        
        %merging, splitting and motion analysis relative to follower presence
        [tracks1SegMobilityMF(trialNum),tracks2plusSegMobilityMF(trialNum),...
            tracksMergeSplitMF(trialNum),tracks1SegFollowerPropMF(trialNum),...
            tracks2plusSegFollowerPropMF(trialNum)] = ...
            analyzeMasterMoveMergeSplitRelFollower(tracksAltTotal,tracksFollowerTotal,...
            diffModeDividerStruct,10);

        %more merging/splitting stats vs. follower presence
        [statsGeneralMF(trialNum),statsFollowerMF(trialNum),...
            statsMergingMF(trialNum),statsSplittingMF(trialNum)] = ...
            calcStatsMSperF(tracksAltTrial,5);

        %labeled receptor on and off rates and oligomer densities
        infoSpaceTime = struct('probDim',2,'timeStep',1,'area',cellMaskArea);
        [rateOnPerClust,rateOffPerClust,densityPerClust,numClustForRateCalc,...
            clustStats(trialNum),rateOnPerClustM,rateOffPerClustM,...
            densityPerClustM,numClustForRateCalcM,clustStatsM(trialNum), ...
            rateOnPerClustF,rateOffPerClustF,densityPerClustF,...
            numClustForRateCalcF,clustStatsF(trialNum)] = ...
            clusterOnOffRatesAndDensityFromClustHistory(clustHistTotalMerge,...
            infoSpaceTime,fracFollowerTotalMerge);

        rateStats(trialNum).rateOnPerClust = rateOnPerClust;
        rateStats(trialNum).rateOffPerClust = rateOffPerClust;
        rateStats(trialNum).densityPerClust = densityPerClust;
        rateStats(trialNum).numClustForRateCalc = numClustForRateCalc;

        rateStatsM(trialNum).rateOnPerClust = rateOnPerClustM;
        rateStatsM(trialNum).rateOffPerClust = rateOffPerClustM;
        rateStatsM(trialNum).densityPerClust = densityPerClustM;
        rateStatsM(trialNum).numClustForRateCalc = numClustForRateCalcM;

        rateStatsF(trialNum).rateOnPerClust = rateOnPerClustF;
        rateStatsF(trialNum).rateOffPerClust = rateOffPerClustF;
        rateStatsF(trialNum).densityPerClust = densityPerClustF;
        rateStatsF(trialNum).numClustForRateCalc = numClustForRateCalcF;

    end %(for numInd : -1 : 1)

    %take care of size inconsistencies
    %tracks with 1 segment
    numRowGlob = length(tracks1SegFollowerPropMF(1).fracTime0Start2Follower);
    for trialNum = 1 : numInd
        numRowLoc = length(tracks1SegFollowerPropMF(trialNum).fracTime0Start2Follower);
        if numRowLoc >= numRowGlob
            tracks1SegFollowerPropMF(trialNum).fracTime0Start2Follower =...
                tracks1SegFollowerPropMF(trialNum).fracTime0Start2Follower(1:numRowGlob);
            tracks1SegFollowerPropMF(trialNum).meanTimeStart2Follower = ...
                tracks1SegFollowerPropMF(trialNum).meanTimeStart2Follower(1:numRowGlob);
            tracks1SegFollowerPropMF(trialNum).meanFracTimeWithFollower =...
                tracks1SegFollowerPropMF(trialNum).meanFracTimeWithFollower(1:numRowGlob);
            tracks1SegFollowerPropMF(trialNum).meanFollowerConsistency = ...
                tracks1SegFollowerPropMF(trialNum).meanFollowerConsistency(1:numRowGlob);
            tracks1SegFollowerPropMF(trialNum).meanTrackDuration = ...
                tracks1SegFollowerPropMF(trialNum).meanTrackDuration(1:numRowGlob);
        else
            tracks1SegFollowerPropMF(trialNum).fracTime0Start2Follower = ...
                [tracks1SegFollowerPropMF(trialNum).fracTime0Start2Follower; NaN(numRowGlob-numRowLoc,1)];
            tracks1SegFollowerPropMF(trialNum).meanTimeStart2Follower = ...
                [tracks1SegFollowerPropMF(trialNum).meanTimeStart2Follower; NaN(numRowGlob-numRowLoc,1)];
            tracks1SegFollowerPropMF(trialNum).meanFracTimeWithFollower = ...
                [tracks1SegFollowerPropMF(trialNum).meanFracTimeWithFollower; NaN(numRowGlob-numRowLoc,1)];
            tracks1SegFollowerPropMF(trialNum).meanFollowerConsistency =...
                [tracks1SegFollowerPropMF(trialNum).meanFollowerConsistency; NaN(numRowGlob-numRowLoc,1)];
            tracks1SegFollowerPropMF(trialNum).meanTrackDuration = ...
                [tracks1SegFollowerPropMF(trialNum).meanTrackDuration; NaN(numRowGlob-numRowLoc,1)];
        end
    end
    %tracks with 2+ segments
    numRowGlob = length(tracks2plusSegFollowerPropMF(1).fracTime0Start2Follower);
    for trialNum = 1 : numInd
        numRowLoc = length(tracks2plusSegFollowerPropMF(trialNum).fracTime0Start2Follower);
        if numRowLoc >= numRowGlob
            tracks2plusSegFollowerPropMF(trialNum).fracTime0Start2Follower = ...
                tracks2plusSegFollowerPropMF(trialNum).fracTime0Start2Follower(1:numRowGlob);
            tracks2plusSegFollowerPropMF(trialNum).meanTimeStart2Follower = ...
                tracks2plusSegFollowerPropMF(trialNum).meanTimeStart2Follower(1:numRowGlob);
            tracks2plusSegFollowerPropMF(trialNum).meanFracTimeWithFollower =...
                tracks2plusSegFollowerPropMF(trialNum).meanFracTimeWithFollower(1:numRowGlob);
            tracks2plusSegFollowerPropMF(trialNum).meanFollowerConsistency = ...
                tracks2plusSegFollowerPropMF(trialNum).meanFollowerConsistency(1:numRowGlob);
            tracks2plusSegFollowerPropMF(trialNum).meanTrackDuration = ...
                tracks2plusSegFollowerPropMF(trialNum).meanTrackDuration(1:numRowGlob);
        else
            tracks2plusSegFollowerPropMF(trialNum).fracTime0Start2Follower = ...
                [tracks2plusSegFollowerPropMF(trialNum).fracTime0Start2Follower; NaN(numRowGlob-numRowLoc,1)];
            tracks2plusSegFollowerPropMF(trialNum).meanTimeStart2Follower = ...
                [tracks2plusSegFollowerPropMF(trialNum).meanTimeStart2Follower; NaN(numRowGlob-numRowLoc,1)];
            tracks2plusSegFollowerPropMF(trialNum).meanFracTimeWithFollower = ...
                [tracks2plusSegFollowerPropMF(trialNum).meanFracTimeWithFollower; NaN(numRowGlob-numRowLoc,1)];
            tracks2plusSegFollowerPropMF(trialNum).meanFollowerConsistency = ...
                [tracks2plusSegFollowerPropMF(trialNum).meanFollowerConsistency; NaN(numRowGlob-numRowLoc,1)];
            tracks2plusSegFollowerPropMF(trialNum).meanTrackDuration = ...
                [tracks2plusSegFollowerPropMF(trialNum).meanTrackDuration; NaN(numRowGlob-numRowLoc,1)];
        end
    end
    
    %%%%ZMALIK 20200604 Add these variables for individual cell analysis
    %with follower - during 
    tmp = cat(3,tracks1SegMobilityMF.fracDiffModeComboTracksWithF);
    tmp1 = squeeze(sum(tmp,1));
    fracDiffModeTracksWithFDur1Seg = tmp1;
    %with follower - before
    tmp1 = squeeze(sum(tmp,2));
    fracDiffModeTracksWithFBef1Seg = tmp1;

    %For each mode before follower, fraction of tracks that swtich to slowest
    %mode during follower
    tmp1 = squeeze(tmp(:,end,:)) ./ tmp1;
    fracMode1DurPerModeBefF1Seg = tmp1;
    
    for indStruct = 1:length(tracks1SegMobilityMF)
        tracks1SegMobilityMF(indStruct).fracDiffModeTracksWithFDur1Seg = ...
            fracDiffModeTracksWithFDur1Seg(:,indStruct);
        tracks1SegMobilityMF(indStruct).fracDiffModeTracksWithFBef1Seg = ...
            fracDiffModeTracksWithFBef1Seg(:,indStruct);
        tracks1SegMobilityMF(indStruct).fracMode1DurPerModeBefF1Seg = ...
            fracMode1DurPerModeBefF1Seg(:,indStruct);
    end
    
    %diffusion coefficient
    %with follower - before
    tmp = cat(3,tracks1SegMobilityMF.meanDiffCoefBeforeModeCombo);
    for i=1:size(tmp,3)
    tmp1(:,i) = diag(tmp(:,:,i));
    end
    diffCoefModeWithFBef1Seg = tmp1;
    %with follower - during 
    tmp = cat(3,tracks1SegMobilityMF.meanDiffCoefDuringModeCombo);
    for i=1:size(tmp,3)
    tmp1(:,i) = diag(tmp(:,:,i));
    end
    diffCoefModeWithFDur1Seg = tmp1;
    
    for indStruct = 1:length(tracks1SegMobilityMF)
    tracks1SegMobilityMF(indStruct).diffCoefModeWithFBef1Seg = ...
        diffCoefModeWithFBef1Seg(:,indStruct);
    tracks1SegMobilityMF(indStruct).diffCoefModeWithFDur1Seg = ...
        diffCoefModeWithFDur1Seg(:,indStruct);
    end
    
end %(if nBS ~= -1)

%% Save all results
if nBS ~= -1
    save([pathName 'combinedFollowerAnalysisRaw_' condName '.mat'],'bsIndxSamp',...
        'statsGeneralMF','statsFollowerMF','statsMergingMF','statsSplittingMF',...
        'tracks1SegMobilityMF','tracks2plusSegMobilityMF','tracksMergeSplitMF',...
        'tracks1SegFollowerPropMF','tracks2plusSegFollowerPropMF','rateStats',...
        'rateStatsM','rateStatsF','clustStats','clustStatsM','clustStatsF','nBS');
else
    save([pathName 'combinedFollowerAnalysisRaw_' condName '.mat'],...
        'statsGeneralMF','statsFollowerMF','statsMergingMF','statsSplittingMF',...
        'tracks1SegMobilityMF','tracks2plusSegMobilityMF','tracksMergeSplitMF',...
        'tracks1SegFollowerPropMF','tracks2plusSegFollowerPropMF','rateStats',...
        'rateStatsM','rateStatsF','clustStats','clustStatsM','clustStatsF','nBS');
end
%% ~~~ the end ~~~
