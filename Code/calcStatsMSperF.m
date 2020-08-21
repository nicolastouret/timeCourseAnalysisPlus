function [statsGeneral,statsFollower,statsMerging,statsSplitting] = ...
    calcStatsMSperF(tracksCell,minTrackLen)
%CALCSTATSMSAlt calculates merge/split statistics for different motion types and modes
%
%SYNOPSIS [statsGeneral,statsMotion,statsMerging,statsSplitting] = ...
%    calcStatsMSAlt(tracksCell,diffAnalysisResCell,diffModeAnalysisResCell,minTrackLen)
%
%INPUT  tracksCell     : Cell array of particle tracks in alternative format, as output by
%                    convTrackFormatDefault2Alt.
%       minTrackLen: Minimum length of a compound track to be used in getting
%                    merge/split statistics.
%                    Optional. Default: 5.
%
%OUTPUT
%
%
%Khuloud Jaqaman, July 2019
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

%% input

if nargin < 1 || isempty(tracksCell)
    disp('calcStatsMSperF: Missing master and/or follower tracks!');
    return
end

if nargin < 2 || isempty(minTrackLen)
    minTrackLen = 5;
end

%% preamble

%take data out of cell array, but record size
cellNum = length(tracksCell);
tracks = vertcat(tracksCell{:});

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
clear criteria

tracks = tracks(indx);

%get total number of tracks and number of frames
numTracks = length(tracks);
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1))*cellNum;

%put tracks in matrix format
[tracksMat,tracksIndxMat,trackStartRow] = convStruct2MatIgnoreMS(tracks);

%get number of track segments
numTrackSegments = size(tracksMat,1);

%% features

%get average number of features per frame
numFeatTot = length(find(tracksIndxMat(:)));
aveFeatPerFrame = numFeatTot / numFrames;

%% follower information

%assess whether follower is present
followerInfo = vertcat(tracks.followerInfo);
followerPresent = followerInfo(:,5) >= 0.1;

%get indices of track segments with and without follower
indxWith    = find( followerPresent );
indxWithout = find( ~followerPresent );

%calculate number of track segments per follower status
numSegmentsType = [length(indxWith) length(indxWithout)]';

%calculate fraction of track segments per follower status
fracSegmentsType = numSegmentsType / sum(numSegmentsType);

%calculate number of features per follower status
numFeatType = [length(find(tracksIndxMat(indxWith,:))) ...
    length(find(tracksIndxMat(indxWithout,:)))]';

%get fraction of features per follower status - this is the probability
%of a feature having or not having follower
probFeatType = numFeatType / sum(numFeatType);

%also calculate average number of features per follower status per frame
aveFeatTypePerFrame = numFeatType / numFrames;

%% merges/splits vs. features

%initialize lists of merges/splits
listOfMerges = [];
listOfSplits = [];

%go over all compound tracks
for iTrack = 1 : numTracks
    
    %get the sequence of events of this compound track
    seqOfEvents = tracks(iTrack).seqOfEvents;
    
    %find where this track has merges and splits
    indxSplit = find( seqOfEvents(:,2)==1 & ~isnan(seqOfEvents(:,4)) );
    indxMerge = find( seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)) );
    
    %store the merge/split time and the segment after merging/before
    %splitting (remember tracks are in alternative format)
    infoSplit = unique(seqOfEvents(indxSplit,[1 4]),'rows'); %#ok<FNDSB>
    infoMerge = unique(seqOfEvents(indxMerge,[1 4]),'rows'); %#ok<FNDSB>
    
    %go over the splits
    for iMS = 1 : size(infoSplit,1)
        
        %find the rows of this split
        rowSplit = find( seqOfEvents(:,1)==infoSplit(iMS,1) & seqOfEvents(:,4)==infoSplit(iMS,2) );
        
        %make a row vector of the segment before splitting and segments after
        %splitting
        segmentsMS = [infoSplit(iMS,2) seqOfEvents(rowSplit,3)']; %#ok<FNDSB>
        
        %get their indices in the global segment matrix
        segmentsMS = segmentsMS + trackStartRow(iTrack) - 1;
        
        %add this split to the list of splits
        listOfSplits = [listOfSplits; segmentsMS]; %#ok<AGROW>
        
    end
    
    %go over the merges
    for iMS = 1 : size(infoMerge,1)
        
        %find the rows of this merge
        rowMerge = find( seqOfEvents(:,1)==infoMerge(iMS,1) & seqOfEvents(:,4)==infoMerge(iMS,2) );
        
        %make a row vector of the segment after merging and segments before
        %merging
        segmentsMS = [infoMerge(iMS,2) seqOfEvents(rowMerge,3)']; %#ok<FNDSB>
        
        %get their indices in the global segment matrix
        segmentsMS = segmentsMS + trackStartRow(iTrack) - 1;
        
        %add this merge to the list of merges
        listOfMerges = [listOfMerges; segmentsMS]; %#ok<AGROW>
        
    end
    
end

%get total number of merges and splits
numMergesTot = size(listOfMerges,1);
numSplitsTot = size(listOfSplits,1);

%calculate average number of merges/splits per frame
aveMergePerFrame = numMergesTot / numFrames;
aveSplitPerFrame = numSplitsTot / numFrames;

%calculate overall probability of merging/splitting
%note that this considers that each merge/split involves two particles
probFeatMerge = 2 * aveMergePerFrame / ( aveFeatPerFrame * (aveFeatPerFrame-1) );
probFeatSplit = 2 * aveSplitPerFrame / ( aveFeatPerFrame * (aveFeatPerFrame-1) );

%% merges/splits vs. follower status

%get the follower status of segments participating in merges and splits
if numMergesTot > 0
    listOfMergeTypes = [followerPresent(listOfMerges(:,1)) ...
        followerPresent(listOfMerges(:,2)) followerPresent(listOfMerges(:,3))];
else
    listOfMergeTypes = zeros(0,3);
end
if numSplitsTot > 0
    listOfSplitTypes = [followerPresent(listOfSplits(:,1)) ...
        followerPresent(listOfSplits(:,2)) followerPresent(listOfSplits(:,3))];
else
    listOfSplitTypes = zeros(0,3);
end

%sort the lists so that the segments with follower before merging or after
%splitting comes first
listOfMergeTypes(:,2:3) = sort(listOfMergeTypes(:,2:3),2,'descend');
listOfSplitTypes(:,2:3) = sort(listOfSplitTypes(:,2:3),2,'descend');

%get number of merges/splits based on the follower status of the three
%track segments involved
typeList = [1 0];
numTypeIndx = length(typeList);
[numMergesTypeDetail,numSplitsTypeDetail] = deal(NaN(numTypeIndx*ones(1,3)));
for kType = 1 : numTypeIndx %segment after merging
    for iType = 1 : numTypeIndx %more dynamic segment before merging
        for jType = iType : numTypeIndx %less dynamic segment before merging
            numMergesTypeDetail(iType,jType,kType) = length(find( ...
                listOfMergeTypes(:,1) == typeList(kType) & ...
                listOfMergeTypes(:,2) == typeList(iType) & ...
                listOfMergeTypes(:,3) == typeList(jType) ));
        end
    end
end
for kType = 1 : numTypeIndx %segment before splitting
    for iType = 1 : numTypeIndx %more dynamic segment after splitting
        for jType = iType : numTypeIndx %less dynamic segment after splitting
            numSplitsTypeDetail(iType,jType,kType) = length(find( ...
                listOfSplitTypes(:,1) == typeList(kType) & ...
                listOfSplitTypes(:,2) == typeList(iType) & ...
                listOfSplitTypes(:,3) == typeList(jType) ));
        end
    end
end

%calculate average merges/splits per frame
aveMergesTypeDetail = numMergesTypeDetail / numFrames;
aveSplitsTypeDetail = numSplitsTypeDetail / numFrames;

%% probability of merging/splitting vs. follower status of segments before merging/after splitting

%collapse merge/split details to retain only folower status of
%segments before merging/after splitting
aveMergesTypeSeparate = sum(aveMergesTypeDetail,3);
aveSplitsTypeSeparate = sum(aveSplitsTypeDetail,3);

%calculate probability of merging/splitting per follower status pair,
%equivalent to calculation above for overall probability
[probMergeTypeSeparate,probSplitTypeSeparate] = deal(NaN(numTypeIndx));
for iType = 1 : numTypeIndx
    probMergeTypeSeparate(iType,iType) = 2 * aveMergesTypeSeparate(iType,iType) / ( aveFeatTypePerFrame(iType) * (aveFeatTypePerFrame(iType)-1) );
    for jType = iType + 1 : numTypeIndx
        probMergeTypeSeparate(iType,jType) = aveMergesTypeSeparate(iType,jType) / ( aveFeatTypePerFrame(iType) * aveFeatTypePerFrame(jType) );
    end
end
for iType = 1 : numTypeIndx
    probSplitTypeSeparate(iType,iType) = 2 * aveSplitsTypeSeparate(iType,iType) / ( aveFeatTypePerFrame(iType) * (aveFeatTypePerFrame(iType)-1) );
    for jType = iType + 1 : numTypeIndx
        probSplitTypeSeparate(iType,jType) = aveSplitsTypeSeparate(iType,jType) / ( aveFeatTypePerFrame(iType) * aveFeatTypePerFrame(jType) );
    end
end

%% probability of merging/splitting vs. follower status of segment after merging/before splitting

%collapse merge/split details to retain only follower status of
%segment after merging/before splitting
aveMergesTypeTogether = squeeze(nansum(nansum(aveMergesTypeDetail,1),2));
aveSplitsTypeTogether = squeeze(nansum(nansum(aveSplitsTypeDetail,1),2));

%calculate probability of merging/splitting per follower status after
%merging/before splitting
probMergeTypeTogether = aveMergesTypeTogether ./ aveFeatTypePerFrame;
probSplitTypeTogether = aveSplitsTypeTogether ./ aveFeatTypePerFrame;

%% probability of follower status after merging/before splitting vs. follower status before merging/after splitting and vice versa

probMergeTypeTogVsSep = aveMergesTypeDetail ./ repmat(aveMergesTypeSeparate,1,1,numTypeIndx);
probSplitTypeTogVsSep = aveSplitsTypeDetail ./ repmat(aveSplitsTypeSeparate,1,1,numTypeIndx);

probMergeTypeSepVsTog = aveMergesTypeDetail ./ repmat(reshape(aveMergesTypeTogether,1,1,numTypeIndx),numTypeIndx,numTypeIndx,1);
probSplitTypeSepVsTog = aveSplitsTypeDetail ./ repmat(reshape(aveSplitsTypeTogether,1,1,numTypeIndx),numTypeIndx,numTypeIndx,1);

%% output

%general statistics
statsGeneral = struct('numTracks',numTracks,'numTrackSegments',numTrackSegments,...
    'numFeatPerFrame',aveFeatPerFrame);

%follower statistics
statsType = struct('fracSegments',fracSegmentsType,...
    'prob',probFeatType);
statsFollower = struct('type',statsType);

%merging statistics
statsType = struct('numTotDetail',numMergesTypeDetail,...
    'numPerFrameDetail',aveMergesTypeDetail,...
    'probVsSepSegments',probMergeTypeSeparate,...
    'probVsTogSegment',probMergeTypeTogether,...
    'probFollowerTogVsSep',probMergeTypeTogVsSep,...
    'probFollowerSepVsTog',probMergeTypeSepVsTog);
statsMerging = struct('numTot',numMergesTot,'numPerFrame',aveMergePerFrame,...
    'probOverall',probFeatMerge,'type',statsType);

%splitting statistics
statsType = struct('numTotDetail',numSplitsTypeDetail,...
    'numPerFrameDetail',aveSplitsTypeDetail,...
    'probVsSepSegments',probSplitTypeSeparate,...
    'probVsTogSegment',probSplitTypeTogether,...
    'probFollowerTogVsSep',probSplitTypeTogVsSep,...
    'probFollowerSepVsTog',probSplitTypeSepVsTog);
statsSplitting = struct('numTot',numSplitsTot,'numPerFrame',aveSplitPerFrame,...
    'probOverall',probFeatSplit,'type',statsType);



%% ~~~ the end ~~~

