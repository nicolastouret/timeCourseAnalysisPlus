function statsGeneral = calcStatsMS_noMotionInfo(tracks,minTrackLen,...
    timeBetweenFrames,removePotArtifacts)
%CALCSTATSMS_NOMOTIONINFO calculates merge/split statistics without distinguishing between different motion types
%
%SYNOPSIS statsGeneral = calcStatsMS_noMotionInfo(tracks,minTrackLen,...
%    timeBetweenFrames,removePotArtifacts)
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       minTrackLen: Minimum length of a track to be used in getting
%                    merge/split statistics.
%                    Optional. Default: 5.
%       timeBetweenFrames: Time between frames (s).
%                    Optional. Devault: 1.
%       removePotArtifacts: 1 to remove potentially artifactual merges and
%                    splits, resulting for instance from detection
%                    artifacts, 0 otherwise. 
%                    Optional. Default: 1.
%
%OUTPUT statsGeneral: Row vector with entries: 
%                     (1) Number of features per frame.
%                     (2) Number of tracks with length >= minTrackLen.
%                     (3) Number of tracks segments belonging to the tracks
%                         with length >= minTrackLen.
%                     (4) Probability of a feature merging.
%                     (5) Probability of a feature splitting.
%                     (6) Rate of a feature merging (per unit time).
%                     (7) Rate of a feature splitting (per unit time).
%                     (8) Rate of a feature merging (per unit time per feature).
%                     (9) Rate of a feature splitting (per unit time per feature).
%                     (10) Probability of a feature merging (per feature).
%                     (11) Probability of a feature splitting (per feature).
%
%Khuloud Jaqaman, March 2013
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

if nargin < 1 || isempty(tracks)
    disp('calcStatsMS: Missing input argument "tracks"!');
    return
end

if nargin < 2 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 3 || isempty(timeBetweenFrames)
    timeBetweenFrames = 1;
end

if nargin < 4 || isempty(removePotArtifacts)
    removePotArtifacts = 1;
end

%% preamble

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
clear criteria
tracks = tracks(indx);

if isempty(indx)
    statsGeneral = [zeros(1,5) NaN(1,4)];
    return
end

%get number of tracks and number of frames
numTracks = length(tracks);
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%put tracks in matrix format
[tracksMat,tracksIndxMat,trackStartRow] = convStruct2MatIgnoreMS(tracks);

%get number of track segments
numTrackSegments = size(tracksMat,1);

%% features

%get average number of features per frame
numFeatTot = length(find(tracksIndxMat(:)));
aveFeatPerFrame = numFeatTot / numFrames;

%% merges/splits vs. features

%initialize lists of merges/splits
listOfMerges = [];
listOfSplits = [];

%go over all compound tracks
for iTrack = 1 : numTracks
    
    %get the sequence of events of this compound track
    seqOfEvents = tracks(iTrack).seqOfEvents;

    %if requested, remove splits and merges that are most likely artifacts
    if removePotArtifacts
        seqOfEvents = removeSplitMergeArtifacts(seqOfEvents,0);
    end

    %find where this track has merges/splits
    indxMS = find(~isnan(seqOfEvents(:,4)));

    %go over these merges/splits
    for iMS = indxMS'

        %determine whether it's a merge (2) or a split (1)
        msType = seqOfEvents(iMS,2);

        %get the indices of the participating segments within the compound
        %track
        segmentsMS = seqOfEvents(iMS,3:4);

        %get their indices in the global segment matrix
        segmentsMS = segmentsMS + trackStartRow(iTrack) - 1;

        %add this merge/split to the list of merges/splits
        if msType == 2 %merge
            listOfMerges = [listOfMerges; segmentsMS]; %#ok<AGROW>
        else %split
            listOfSplits = [listOfSplits; segmentsMS]; %#ok<AGROW>
        end

    end

end

%get total number of merges and splits
numMergesTot = size(listOfMerges,1);
numSplitsTot = size(listOfSplits,1);

%calculate average number of merges/splits per frame
aveMergePerFrame = numMergesTot / numFrames;
aveSplitPerFrame = numSplitsTot / numFrames;

%calculate average number of merges/splits per feature - this is the
%probability of a feature merging/splitting
probFeatMerge = aveMergePerFrame / aveFeatPerFrame;
probFeatSplit = aveSplitPerFrame / aveFeatPerFrame;

%calculate the rate of a feature merging/splitting per unit time
rateFeatMerge = probFeatMerge / timeBetweenFrames;
rateFeatSplit = probFeatSplit / timeBetweenFrames;

%calculate the rate per unit time per feature
%KJ 200716: updated rate calculations to make them consistent with probability calculations in calcStatsMSAlt
rateFeatMergeNorm = rateFeatMerge / ((aveFeatPerFrame-1)/2); 
rateFeatSplitNorm = rateFeatSplit / ((aveFeatPerFrame-1)/2);

%calculate the probability per feature
%KJ 200716: added this probability calculation, which is consitent with the
%probability calculations in calcStatsMSAlt
probFeatMergeNorm = probFeatMerge / ((aveFeatPerFrame-1)/2); 
probFeatSplitNorm = probFeatSplit / ((aveFeatPerFrame-1)/2);

%% output

%general statistics
statsGeneral = [aveFeatPerFrame numTracks numTrackSegments ...
    probFeatMerge probFeatSplit rateFeatMerge rateFeatSplit ...
    rateFeatMergeNorm rateFeatSplitNorm probFeatMergeNorm probFeatSplitNorm];

%% ~~~ the end ~~~

