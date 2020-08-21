function [tracksMatNoF,tracksMatIntermF,tracksMatWithFBef,tracksMatWithFAft,...
    tracksMatWithFDuring,tracksMatWithFBetween,oligoInfo] = ...
    separateMasterTracksBasedOnFollower(tracksMat,followerIndxMat,oligoMat)
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

oligoFlag = 1;
if nargin < 3 || isempty(oligoMat)
    oligoFlag = 0;
    oligoMat = followerIndxMat;
end
oligoVec = max(oligoMat,[],2);

%% Analysis

%separate all master tracks that do not have any follower
sumFollowerMat = sum(followerIndxMat,2);
logicalFnoF = (sumFollowerMat ~= 0);
tracksMatNoF = tracksMat(~logicalFnoF,:);
oligoVecNoF = oligoVec(~logicalFnoF);

%keep the rest for further subdivision
tracksMat = tracksMat(logicalFnoF,:);
followerIndxMat = followerIndxMat(logicalFnoF,:);
sumFollowerMat = sumFollowerMat(logicalFnoF);
oligoVec = oligoVec(logicalFnoF);

%get the length of each track to compare to number of follower instances
%also get start time to compare to follower appearance time
trackLft = getTrackSEL(tracksMat);
trackStart = trackLft(:,1);
trackLft = trackLft(:,3);

%tracks that have <= 2 follower instances should be separated from the rest,
%unless the tracks are very short and have a follower for more than half of
%their lifetime
logicalFnoF = ( (sumFollowerMat > 2) | (sumFollowerMat./trackLft >= 0.5) );
tracksMatIntermF = tracksMat(~logicalFnoF,:);
oligoVecIntermF = oligoVec(~logicalFnoF);
    
%keep the rest for further subdivision
tracksMat = tracksMat(logicalFnoF,:);
followerIndxMat = followerIndxMat(logicalFnoF,:);
sumFollowerMat = sumFollowerMat(logicalFnoF);
oligoVec = oligoVec(logicalFnoF);
trackStart = trackStart(logicalFnoF);
trackLft = trackLft(logicalFnoF);

%now get the first and last follower instance for each track
numTracks = length(sumFollowerMat);
[firstInstance,lastInstance] = deal(NaN(numTracks,1));
for iTrack = 1 : numTracks
    firstInstance(iTrack) = find(followerIndxMat(iTrack,:),1,'first');
    lastInstance(iTrack) = find(followerIndxMat(iTrack,:),1,'last');
end

%then calculate time from first to last instance, and fraction of time spent
%with follower
first2lastInstance = lastInstance - firstInstance + 1;
fracFollower = sumFollowerMat ./ first2lastInstance;

%tracks with fracFollower < 0.1 should be added to the intermediate group
logicalFnoF = (fracFollower >= 0.1);
tracksMatIntermF = [tracksMatIntermF; tracksMat(~logicalFnoF,:)];
oligoVecIntermF = [oligoVecIntermF; oligoVec(~logicalFnoF)];

%keep the rest for further subdivision
tracksMat = tracksMat(logicalFnoF,:);
followerIndxMat = followerIndxMat(logicalFnoF,:);
firstInstance = firstInstance(logicalFnoF);
lastInstance = lastInstance(logicalFnoF);
fracFollower = fracFollower(logicalFnoF);
sumFollowerMat = sumFollowerMat(logicalFnoF);
oligoVec = oligoVec(logicalFnoF);
trackStart = trackStart(logicalFnoF);
trackLft = trackLft(logicalFnoF);
numTracks = length(firstInstance);

%subdivide to before and after follower
%WithFBef: with follower, before first instance
%WithFAft: with follower, after last instance
[tracksMatWithFBef,tracksMatWithFAft] = deal(tracksMat);
for iTrack = 1 : numTracks
    tracksMatWithFBef(iTrack,(firstInstance(iTrack)-1)*8+1:end) = NaN;
    tracksMatWithFAft(iTrack,1:lastInstance(iTrack)*8) = NaN;
end

%subdivide to during and between follower
%WithFDuring: with follower, at times when follower is there
%WithFBetween: with follower, but at times when follower is not there, between first instance and last instance

%expand followerIndxMat mat to have 8 columns per frame
followerIndxMat8 = repmat(followerIndxMat,1,8);
for iCol = 1 : 8
    followerIndxMat8(:,iCol:8:end) = followerIndxMat;
end
followerIndxMat8 = double(followerIndxMat8);
followerIndxMat8(followerIndxMat8==0) = NaN;

%multiply with tracksMat to retain relevant frames
tracksMatWithFDuring = tracksMat .* followerIndxMat8;

%use all of the above to obtain tracksWithFBetween
tracksMatWithFBetween = tracksMat;
tracksMatWithFBetween(~isnan(tracksMatWithFBef)) = NaN;
tracksMatWithFBetween(~isnan(tracksMatWithFAft)) = NaN;
tracksMatWithFBetween(~isnan(tracksMatWithFDuring)) = NaN;

%store oligomeric state of each track, and when follower appeared relative
%to the track's start time
if oligoFlag
    timeStart2FollowerAppear = firstInstance - trackStart;
    oligoInfo.noF = oligoVecNoF;
    oligoInfo.intermF = oligoVecIntermF;
    oligoInfo.withF = [oligoVec timeStart2FollowerAppear sumFollowerMat./trackLft fracFollower trackLft];
else
    oligoInfo = [];
end


