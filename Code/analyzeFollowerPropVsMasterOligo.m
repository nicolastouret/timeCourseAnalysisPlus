function [fracTime0Start2Follower,meanTimeStart2Follower,meanFracTimeWithFollower,...
    meanFollowerConsistency,meanTrackDuration,timeStart2Follower, ...
    fracTimeWithFollower,followerConsistency,trackDuration] = ...
    analyzeFollowerPropVsMasterOligo(tracksMF,tracksF)
%ANALYZEFOLLOWERPROPVSMASTEROLIGO analyzes follower properties vs. oligomerization state of master channel
%
%SYNPOSIS [fracTime0Start2Follower,meanTimeStart2Follower,meanFracTimeWithFollower,...
%    meanFollowerConsistency,meanTrackDuration,timeStart2Follower, ...
%    fracTimeWithFollower,followerConsistency,trackDuration] = ...
%    analyzeFollowerPropVsMasterOligo(tracksMF,tracksF)
%
%INPUT  tracksMF, tracksF      : Output of getFollowerTracksFromMaster.
%                                MUST BE IN ALTERNATIVE FORMAT (see
%                                getFollowerTracksFromMaster for format details).
%
%OUTPUT fracTime0Start2Follower: Fraction of cases where follower
%                                appearance time relative to track start
%                                time = 0, per oligomeric state.
%       meanTimeStart2Follower : Mean follower appearance time, per
%                                oligomeric state.
%       meanFracTimeWithFollower: Mean fraction of time with follower
%                                relative to each track's duration, per
%                                oligomeric state.
%       meanFollowerConsistency: Mean follower consistency, i.e. fraction
%                                of time that follower is present between
%                                its appearance and disappearance, per
%                                oligomeric state.
%       meanTrackDuration      : Mean length of time tracks last 
%       timeStart2Follower     : Cell array storing the distribution of
%                                follower appearance times, per
%                                oligomeric state. Follower appearance time
%                                is relative to track start time (so if a
%                                track has follower from its very
%                                beginning, the time will be 0).
%       fracTimeWithFollower   : Cell array storing the distribution of
%                                fraction of time with follower relative to
%                                each track's duration, per oligomeric
%                                state.
%       fracFollowerConsistency: Cell array storing the distribution of
%                                follower consistency, per oligomeric
%                                state.
%       trackDuration          : Cell array storing length of time tracks
%                                lasted for
%       
%Khuloud Jaqaman, July 2019
%Zachariah Malik, May 2020, account for case of no followers
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

%% Track separation and oligomeric state information

%convert to matrix format for ease of handling
[tracksMat,~,~,~,oligoMat] = convStruct2MatIgnoreMS(tracksMF);
[~,followerIndxMat] = convStruct2MatIgnoreMS(tracksF);
followerIndxMat = (followerIndxMat ~= 0); %logical indicating whether follower exists (1) or not (0)

%get oligomeric state information in relation to follower properties
[~,~,~,~,~,~,oligoInfo] = separateMasterTracksBasedOnFollower(tracksMat,...
    followerIndxMat,oligoMat);

% %oligomeric state information for tracks without follower
% tmp = oligoInfo.noF;
% if ~isempty(tmp)
%     tmp = tmp(:,1);
%     fracTracksOligoNoF = hist(tmp,1:max(tmp));
%     fracTracksOligoNoF = fracTracksOligoNoF / sum(fracTracksOligoNoF);
% else
%     fracTracksOligoNoF = NaN;
% end

%now for tracks with follower ...
oligoFollower = oligoInfo.withF;

if ~isempty(oligoFollower)
    
    %maximum oligomeric state
    maxOligo = max(oligoFollower(:,1));
    
    %allocate memory
    [timeStart2Follower,fracTimeWithFollower,followerConsistency,...
        trackDuration] = deal(cell(maxOligo,1));
    [fracTime0Start2Follower,meanTimeStart2Follower,meanFracTimeWithFollower,...
        meanFollowerConsistency,meanTrackDuration] = deal(NaN(maxOligo,1));
    
    %     %oligomeric state distribution
    %     fracTracksOligoWithF = hist(oligoFollower(:,1),1:maxOligo);
    %     fracTracksOligoWithF = fracTracksOligoWithF / sum(fracTracksOligoWithF);
    
    for iOligo = 1 : maxOligo
        
        %tracks with this oligomeric state
        indxOligo = oligoFollower(:,1)==iOligo;
        
        %their duration
        tmp = oligoFollower(indxOligo,5);
        trackDuration{iOligo} = tmp;
        meanTrackDuration(iOligo) = mean(tmp);
        
        %their follower appearance time
        tmp = oligoFollower(indxOligo,2);
        timeStart2Follower{iOligo} = tmp;
        fracTime0Start2Follower(iOligo) = length(find(tmp==0)) / length(tmp);
        meanTimeStart2Follower(iOligo) = mean(tmp);
        
        %their fraction with follower
        tmp = oligoFollower(indxOligo,3);
        fracTimeWithFollower{iOligo} = tmp;
        meanFracTimeWithFollower(iOligo) = mean(tmp);
        
        %their follower consistency (i.e. fraction of time follower is
        %there between its appearance and disappearance
        tmp = oligoFollower(indxOligo,4);
        followerConsistency{iOligo} = tmp;
        meanFollowerConsistency(iOligo) = mean(tmp);
        
    end

%%%%ZMALIK 20200528 Consider when no oligomers have a follower
else

    timeStart2Follower = {NaN};
    trackDuration = {NaN};
    fracTimeWithFollower = {NaN};
    followerConsistency = {NaN}; 
    meanTimeStart2Follower = NaN;
    meanTrackDuration = NaN;  
    fracTime0Start2Follower = NaN;
    meanFracTimeWithFollower = NaN;
    meanFollowerConsistency = NaN;

%%%%    
   
    
    % else
    %
    %     [fracTracksOligoWithF,fracTime0Start2Follower,meanTimeStart2Follower] = deal(NaN);
    %     timeStart2Follower = cell(1);
    
end

