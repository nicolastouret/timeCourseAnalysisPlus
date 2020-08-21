function [tracksM,tracksF] = getFollowerTracksFromMaster(tracksM,movieInfoF,...
    matchRad,cellMask)
%GETFOLLOWERTRACKSFROMMASTER extracts tracks of follower channel by association with master channel
%
%SYNPOSIS [tracksM,tracksF] = getFollowerTracksFromMaster(tracksM,movieInfoF,...
%    matchRad,cellMask)
%
%INPUT  tracksM     : Master channel tracks, in format of output of
%                     trackCloseGapsKalmanSparse (so-called default format)
%                     or in alternative format after running
%                     convTrackFormatDefault2Alt.
%       movieInfoF  : Follower channel detections, in format of output of
%                     detectSubResFeatures2D_StandAlone.
%       matchRad    : Matching radius between master and follower channel
%                     positions.
%       cellMask    : Cell mask, to analyze only tracks within it. 
%                     Optional. Default: no mask.
%
%OUTPUT tracksM     : Same as input, and in same format (default or
%                     alternative), but only for tracks within cell mask,
%                     and with the following added fields:
%           .followerFeatIndx: Equivalent to "tracksFeatIndxCG," but storing
%                              the associated follower feature indices. 0
%                              indicates no follower.
%           .followerCoordAmp: Equivalent to "tracksCoordAmpCG," but
%                              storing coordinates and amplitude
%                              information of associated follower. NaN
%                              indicates no follower.
%           .followerInfo    : Follower statistics for each segment in a
%                              track: Number of frames with follower, first
%                              frame with follower, last frame with
%                              follower, Lifetime with follower
%                              (last-first+1), and ratio of number of
%                              frames with follower to lifetime with
%                              follower.
%       tracksF     : Follower tracks, as simply deduced from association
%                     with master tracks. Same format as tracksM.
%
%Khuloud Jaqaman, June 2019
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

if nargin < 4 || isempty(cellMask)
    cellMask = [];
end

if ~isempty(cellMask)
    
    %get the mean position for each track
    numTracks = length(tracksM);
    coordMean = NaN(numTracks,2);
    for iTrack = 1 : length(tracksM)
        xCoord = tracksM(iTrack).tracksCoordAmpCG(:,1:8:end);
        yCoord = tracksM(iTrack).tracksCoordAmpCG(:,2:8:end);
        coordMean(iTrack) = [nanmean(xCoord(:)) nanmean(yCoord(:))];
    end
    coordLinIndx = sub2ind(size(cellMask),round(coordMean(:,1)),round(coordMean(:,2))); %check coordinates, image vs. matrix
    
    %keep only tracks within mask
    indxKeep = cellMask(coordLinIndx);
    tracksM = tracksM(indxKeep);

end
    
%% Analysis

%convert tracks to matrix format to fill in follower information
[trackInfoM,~,trackStartRow,numSeg] = convStruct2MatIgnoreMS(tracksM);
xCoordM = trackInfoM(:,1:8:end);
yCoordM = trackInfoM(:,2:8:end);
[numRow,numFrames] = size(xCoordM);
trackMatF = NaN(numRow,8*numFrames);
trackIdxF = zeros(numRow,numFrames);

%go over each frame, match follower to master, and store information
for iFrame = 1 : numFrames
    
    %get positions in master channel
    coordM = [xCoordM(:,iFrame) yCoordM(:,iFrame)];
    
    %proceed only if master and follower channels are not empty
    if ~isempty(~isnan(coordM)) && ~isempty(movieInfoF(iFrame).xCoord)
        
        %get positions in follower channel
        coordF = [movieInfoF(iFrame).xCoord(:,1) movieInfoF(iFrame).yCoord(:,1)];
        
        %match master and follower positions
        [idx1, idx2] = colocalizationLAP(coordM,coordF,matchRad);
        numMatch = length(idx1);
        
        %fill in follower information associated with master positions
        trackIdxF(idx1,iFrame) = idx2;
        trackMatF(idx1,(iFrame-1)*8+1:iFrame*8) = [coordF(idx2,:) ...
            zeros(numMatch,1) movieInfoF(iFrame).amp(idx2,1) ...
            movieInfoF(iFrame).xCoord(idx2,2) movieInfoF(iFrame).yCoord(idx2,2) ...
            zeros(numMatch,1) movieInfoF(iFrame).amp(idx2,2)];
        
    end
    
end

%add follower information to master tracks
%also make separate follower tracks in case needed later
tracksF = tracksM;
for iTrack = 1 : length(tracksM)
    
    startFrame = min(tracksM(iTrack).seqOfEvents(:,1));
    endFrame = max(tracksM(iTrack).seqOfEvents(:,1));
    startRowCurr = trackStartRow(iTrack);
    
    %this is equivalent to the usual track information
    tracksM(iTrack).followerFeatIndx = trackIdxF(startRowCurr:startRowCurr+numSeg(iTrack)-1,startFrame:endFrame);
    tracksM(iTrack).followerCoordAmp = trackMatF(startRowCurr:startRowCurr+numSeg(iTrack)-1,(startFrame-1)*8+1:endFrame*8);
    
    %information on when there is a follower, if at all
    [positiveFrames,firstFrame,lastFrame] = deal(NaN(numSeg(iTrack),1));
    for iSeg = 1 : numSeg(iTrack)
        tmp = find( tracksM(iTrack).followerFeatIndx(iSeg,:) ~= 0 );
        positiveFrames(iSeg) = length(tmp);
        if ~isempty(tmp)
            firstFrame(iSeg) = tmp(1);
            lastFrame(iSeg) = tmp(end);
        end
    end
    totalFrames = lastFrame - firstFrame + 1;
    tracksM(iTrack).followerInfo = [positiveFrames firstFrame lastFrame totalFrames positiveFrames./totalFrames];
    
    %copy follower tracks into own structure as well
    tracksF(iTrack).tracksFeatIndxCG = tracksM(iTrack).followerFeatIndx;
    tracksF(iTrack).tracksCoordAmpCG = tracksM(iTrack).followerCoordAmp;
    
end

