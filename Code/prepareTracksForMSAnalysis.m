function prepareTracksForMSAnalysis(nameCond,summary,diffModeDividerStruct,...
    timeStart,timeEnd,followerFlag,CML)
%PREPARETRACKSFORMSANALYSIS assembles tracks from time course analysis and prepares them for detailed merge/split/motion analysis
%
% SYNOPSIS: prepareTracksForMSAnalysis(nameCond,summary,diffModeDividerStruct,...
%    timeStart,timeEnd,followerFlag,CML)
%
% INPUT:
%
%   nameCond: Name of condition being processed, used for saving results.
%
%   summary: Cell from analysisData.mat with time course analysis data
%   for condition being processed.
%
%   diffModeDividerStruct: Structure containing threshold values to classify
%   diffusion mode of track.
%
%   timeStart: Single value for the earliest time point to analyse.
%
%   timeEnd: Single value for the latest time point to analyse.
%
%   followerFlag: 1 to extract follower info (for master-follower
%   experiments), 0 otherwise. Optional. Default: 0.
%
%   CML: Combined movie list. 
%
% OUTPUT: No output arguments, but the following variables are saved in a
% file named tracksEtcForMSAnalysis_nameCond_timeStart-timeEnd.mat, in the
% directory where the function is called. Saved variables:
%
%   tracksAltCell: Cell array of the tracks, in alternative format, for the
%   different movies assmebled together.
%
%   tracksDefCell: Cell array of the tracks, in default format, for the
%   different movies assembled together.
%
%   diffTypeCell: Cell array of the tracks' MSS diffusion analysis for the
%   different movies assembled together. Tracks analyzed are in alternative
%   format.
%
%   diffModeCell: Cell array of the tracks' diffusion mode analysis for the
%   different movies assembled together. Tracks analyzed are in alternative
%   format.
%
%   clustHistCell: Cluster history from the tracks for the different movies
%   assembled together.
%
%   cellMaskArea: 1D array of cell mask area per movie.
%
%   tracksFollowerCell: Cell array of follower tracks, in alternative
%   format, for the different movies assembled together. Only relevant for
%   master-follower analysis.
%
%   fracFollowerCell: Cell array indicating fraction of follower presence
%   per master track/cluster, for the different movies assembled together.
%   Only relevant for master-follower analysis.
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

%% Input

%check optional followerFlag input
if nargin < 6 || isempty(followerFlag)
    followerFlag = 0;
end

%Get all time points
timePts = summary{1}.time;

%Select min time point to analyze. First time point if empty
if nargin < 4 || isempty(timeStart)
    timeStart = min(timePts);
end
minTime = timeStart;

%Select max time point to analyze. Last time point + 1 if empty (the +1
%is to make sure everything is taken in the end, because of the "less than"
%statement)
if nargin < 5 || isempty(timeEnd)
    timeEnd = max(timePts) + 1;
end
maxTime = timeEnd;

%Find movies within time range
ind = find(timePts < maxTime & timePts >= minTime );
numInd = length(ind);

%get cell mask area for these movies
cellMaskArea = summary{1}.numAbsMode(:,1)./summary{1}.densityAbsMode(:,1);
cellMaskArea = cellMaskArea(ind);

%get the list of sorted movies, and append movie number within movie
%list to movie list number
mlIndex = summary{1}.mlIndex;
mlMax = max(mlIndex);
mdIndex = mlIndex;
for iML = 1 : mlMax
    mlEntries = find(mlIndex==iML);
    mdIndex(mlEntries) = (1:length(mlEntries))';
end
mlIndex = [mlIndex mdIndex];
mlIndex = mlIndex(ind,:);

%Initialize various variables
[tracksAltCell,tracksDefCell,diffTypeCell,diffModeCell,clustHistCell,...
    tracksFollowerCell,fracFollowerCell] = deal(cell(numInd,1));

%% Collect tracks and their basic analyses

%Go through each movie in time interval
for k = 1 : numInd
    
    %Get tracks
    tracksFinal = summary{1}.tracks{ind(k)};
    
    % 1. Remove simultaneous merges and splits
    tracksReform1 = removeSimultaneousMergeSplit(tracksFinal);
    
    % 2. reshape the seqOfEvents so that shorter-lived segments merge
    % with or split from longer lived segments
    tracksReform2 = reformatSplitsMergesKeepLongerLivedSegment(tracksReform1);
    
    % 3.Reform tracks to remove  artifacts coming from start to merge, split to end or
    % split to merge.
    tracksReform3 = removeSplitMergeArtifactsChronological(tracksReform2,[]);
    
    %get intensity of individual fluorophore for this movie,
    %taken as the intensity of Mode 2 / 2
    intensityInfo = summary{1}.ampStatsF20(ind(k),[3 5]);
    intensityInfo(1) = intensityInfo(1)/2;
    
    %estimate aggregation state
    %this also gives tracks in alternative format
    [compTracksOut] = aggregStateFromCompTracksMIQP(tracksReform3,intensityInfo,10);
    compTracksDef = compTracksOut.defaultFormatTracks;
    compTracksAlt = compTracksOut.alternativeFormatTracks;
    
    %MSS Diffusion Classification for tracks in alternative format
    [diffAnalysisRes] = trackDiffusionAnalysis1(compTracksAlt,[],2,[],[],0,[]);
    
    %Diffusion Mode analysis for tracks in alternative format
    diffModeAnalysisRes = trackDiffModeAnalysis(compTracksAlt,diffModeDividerStruct);
    
    %master-follower analysis
    if followerFlag
        MD = CML.movieLists_(mlIndex(k,1)).movies_{mlIndex(k,2)};
        iLocProc = MD.getProcessIndex('DetectionProcess',1,1);
        locProc = MD.processes_{iLocProc};
        movieInfoF = locProc.loadChannelOutput(2);
        [compTracksAlt,tracksFollAlt] = getFollowerTracksFromMaster(compTracksAlt,movieInfoF,4);
    end
    
    %Get full cluster history of all events, document diffusion modes
    [clustHistoryAll,~,followerFracAll] = clusterHistoryFromCompTracks_aggregState_motion(compTracksDef,compTracksAlt,diffModeAnalysisRes);
    
    % Store information
    tracksAltCell{k} = compTracksAlt;
    tracksDefCell{k} = compTracksDef;
    diffTypeCell{k}  = diffAnalysisRes;
    diffModeCell{k}  = diffModeAnalysisRes;
    clustHistCell{k} = clustHistoryAll;
    if followerFlag
        tracksFollowerCell{k} = tracksFollAlt;
        fracFollowerCell{k} = followerFracAll;
    end
    
end %(for k = 1 : numInd)

% Save all results
save(['tracksEtcForMSAnalysis_' nameCond '_' num2str(timeStart) '-' num2str(timeEnd)  '.mat'],...
    'tracksAltCell','tracksDefCell','diffTypeCell','diffModeCell','clustHistCell',...
    'cellMaskArea','tracksFollowerCell','fracFollowerCell');


%% ~~~ the end ~~~



