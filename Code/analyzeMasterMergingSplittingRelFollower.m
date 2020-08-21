function [fracFollowerBefAmongMergesWithFollowerAft,fracMonoNoFollowerBefMerge,...
    fracMonoWithFollowerBefMerge,meanOligoNoFollowerBefMerge,...
    meanOligoWithFollowerBefMerge,infoMergeAftBef,infoSplitBefAft] = ...
    analyzeMasterMergingSplittingRelFollower(tracksMF,tracksF,timeTolFromMerge)
%ANALYZEMASTERMERGINGSPLITTINGRELFOLLOWER analyzes merging and splitting of master channel molecules relative to follower molecule presence
%
%SYNPOSIS [fracFollowerBefAmongMergesWithFollowerAft,fracMonoNoFollowerBefMerge,...
%    fracMonoWithFollowerBefMerge,meanOligoNoFollowerBefMerge,...
%    meanOligoWithFollowerBefMerge,infoMergeAftBef,infoSplitBefAft] = ...
%    analyzeMasterMergingSplittingRelFollower(tracksMF,tracksF,timeTolFromMerge)
%        
%
%INPUT  tracksMF, tracksF    : Output of getFollowerTracksFromMaster.
%                              MUST BE IN ALTERNATIVE FORMAT (see
%                              getFollowerTracksFromMaster for format details).
%       timeTolFromMerge     : Allowed time difference (in frames) between
%                              last presence of follower with a merging
%                              segment and its merge time.
%                              Optional. Default: 10.
%
%OUTPUT fracFollowerBefAmongMergesWithFollowerAft: Among the merge events
%                        with follower presence after the merge, the
%                        fraction of merges with follower presence also
%                        before the merge.
%       fracMonoNoFollowerBefMerge: Among the merge events with follower
%                        presence after the merge and no follower presence
%                        before the merge, the fraction of segments
%                        classified as monomers. The two entries refer to
%                        the first and second segments. While both are
%                        considered to have no follower before the merge
%                        because of the time tolerance, the first segment
%                        might have overall higher follower presence than
%                        the second segment (but not close enough to the
%                        merge event).
%       fracMonoWithFollowerBefMerge: Among the merge events with follower
%                        presence after the merge and with follower presence
%                        before the merge, the fraction of segments
%                        classified as monomers. The two entries refer to
%                        the first and second segments. The first segment
%                        is the one with follower presence close to the
%                        merge event. The second segment generally has no
%                        follower presence. In the rate cases where both
%                        segments have follower close to the merge event,
%                        the first segment has more reliable follower
%                        presence.
%       meanOligoNoFollowerBefMerge: Mean oligomeric state for same
%                        segments as fracMonoNoFollowerBefMerge.
%       meanOligoWithFollowerBefMerge: Mean oligomeric state for same
%                        segments as fracMonoWithFollowerBefMerge.
%       infoMergeAftBef: 3x4x(number of merges) array storing follower
%                        information. For each merge (indicated by 3rd
%                        dimension), rows include information for segment
%                        after merging (row 1), the segment before merging
%                        with more consistent follower presence (row 2),
%                        and the segment before merging with less
%                        consistent follower presence (row 3). Columns
%                        store the merging time, first appearance of
%                        follower relative to merging time, last appearance
%                        of follower relative to merging time, and fraction
%                        of time between first and last appearance where
%                        follower is in fact present. A higher fraction
%                        indicates more consistent follower presence.
%       infoSplitBefAft: Same as inforMergeAftBef, but for splits. The main
%                        difference is that row 1 is for segment before
%                        splitting, while rows 2 and 3 are for segments
%                        after splitting.
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

%% Input

if nargin < 3 || isempty(timeTolFromMerge)
    timeTolFromMerge = 10;
end

%determine tracks with > 1 segment
numSegments = getNumSegmentsPerTrack(tracksMF);
indx2plus = find(numSegments > 1);

%% Survey of merging and splitting events

%allocate memory and initialize indexing variables
[infoMergeAftBef,infoSplitBefAft] = deal(NaN(3,5,length(indx2plus)));
iGlobalMerge = 0;
iGlobalSplit = 0;

%go over tracks with M/S and analyze follower presence relative to M/S time
for iTrack = indx2plus'
    
    %collect follower information (condensed form of what is done in
    %separateMasterTracksBasedOnFollower function)
    [~,followerIndxMat] = convStruct2MatIgnoreMS(tracksF(iTrack));
    followerIndxMat = (followerIndxMat ~= 0);
    sumFollowerMat = sum(followerIndxMat,2);
    segLft = getTrackSEL(tracksMF(iTrack).tracksCoordAmpCG);
    segLft = segLft(:,3);
    indx0 = sumFollowerMat <= 2 & sumFollowerMat./segLft < 0.5; %unclear follower presence
    sumFollowerMat(indx0) = 0;
    numSegs = length(sumFollowerMat);
    [firstInstance,lastInstance] = deal(NaN(numSegs,1));
    for iSeg = find(sumFollowerMat'~=0)
        firstInstance(iSeg) = find(followerIndxMat(iSeg,:),1,'first');
        lastInstance(iSeg) = find(followerIndxMat(iSeg,:),1,'last');
    end
    first2lastInstance = lastInstance - firstInstance + 1;
    fracFollower = sumFollowerMat ./ first2lastInstance;
    fracFollower(isnan(fracFollower)) = 0;
    fracFollower(fracFollower<0.1) = 0; %unclear follower presence
    
    %get master oligomeric state
    oligoMaster = max(tracksMF(iTrack).aggregState,[],2);
    
    %find M/S events
    seqOfEvents = tracksMF(iTrack).seqOfEvents;
    tmp = seqOfEvents(~isnan(seqOfEvents(:,4)),:);
    eventInfo = [tmp(1:2:end,1:3) tmp(2:2:end,3) tmp(1:2:end,4)];
    
    %collect follower information before and after events
    %also collect master oligomeric state information
    for iEvent = 1 : size(eventInfo,1)
        
        %if event is a merge
        if eventInfo(iEvent,2) == 2
            
            %get segments involved and merge time
            iSegBef1 = eventInfo(iEvent,3);
            iSegBef2 = eventInfo(iEvent,4);
            iSegAft  = eventInfo(iEvent,5);
            mergeTime = eventInfo(iEvent,1);
            
            %in case seg2 has larger fraction of follower, switch seg1 and seg2
            if fracFollower(iSegBef2) > fracFollower(iSegBef1)
                tmp = iSegBef1;
                iSegBef1 = iSegBef2;
                iSegBef2 = tmp;
            end
            
            %collect information
            iGlobalMerge = iGlobalMerge + 1;
            infoMergeAftBef(1,:,iGlobalMerge) = [mergeTime firstInstance(iSegAft)-mergeTime lastInstance(iSegAft)-mergeTime fracFollower(iSegAft) oligoMaster(iSegAft)];
            infoMergeAftBef(2,:,iGlobalMerge) = [mergeTime firstInstance(iSegBef1)-mergeTime lastInstance(iSegBef1)-mergeTime fracFollower(iSegBef1) oligoMaster(iSegBef1)];
            infoMergeAftBef(3,:,iGlobalMerge) = [mergeTime firstInstance(iSegBef2)-mergeTime lastInstance(iSegBef2)-mergeTime fracFollower(iSegBef2) oligoMaster(iSegBef2)];
            
        else %if event is a split
            
            %get segments involved and split time
            iSegAft1 = eventInfo(iEvent,3);
            iSegAft2 = eventInfo(iEvent,4);
            iSegBef  = eventInfo(iEvent,5);
            splitTime = eventInfo(iEvent,1);
            
            %in case seg2 has larger fraction of follower, switch seg1 and seg2
            if fracFollower(iSegAft2) > fracFollower(iSegAft1)
                tmp = iSegAft1;
                iSegAft1 = iSegAft2;
                iSegAft2 = tmp;
            end
            
            %collect information
            iGlobalSplit = iGlobalSplit + 1;
            infoSplitBefAft(1,:,iGlobalSplit) = [splitTime firstInstance(iSegBef)-splitTime lastInstance(iSegBef)-splitTime fracFollower(iSegBef) oligoMaster(iSegBef)];
            infoSplitBefAft(2,:,iGlobalSplit) = [splitTime firstInstance(iSegAft1)-splitTime lastInstance(iSegAft1)-splitTime fracFollower(iSegAft1) oligoMaster(iSegAft1)];
            infoSplitBefAft(3,:,iGlobalSplit) = [splitTime firstInstance(iSegAft2)-splitTime lastInstance(iSegAft2)-splitTime fracFollower(iSegAft2) oligoMaster(iSegAft2)];
            
        end
        
        
    end
    
end

%remove unfilled entries
infoMergeAftBef = infoMergeAftBef(:,:,~isnan(squeeze(infoMergeAftBef(1,1,:))));
infoSplitBefAft = infoSplitBefAft(:,:,~isnan(squeeze(infoSplitBefAft(1,1,:))));

%% Analysis

% % remove entries with no follower presence after cleanup
% infoMergeAftBef = infoMergeAftBef(:,:,squeeze(sum(infoMergeAftBef(:,4,:),1))~=0);
% infoSplitBefAft = infoSplitBefAft(:,:,squeeze(sum(infoSplitBefAft(:,4,:),1))~=0);

%find tracks that have follower after merging
indxFollowerAftMerge = squeeze(infoMergeAftBef(1,4,:)) > 0;
infoMergeAftBefFollower = infoMergeAftBef(:,:,indxFollowerAftMerge);
numFollowerAftMergeTot = sum(indxFollowerAftMerge);

%among those, determine when they have follower before merging (NaN means
%no follower)
timeFollowerBefMerge = [squeeze(infoMergeAftBefFollower(2,3,:)) squeeze(infoMergeAftBefFollower(3,3,:))];

%also get the oligomeric states before merging
oligoStateBefMerge = [squeeze(infoMergeAftBefFollower(2,5,:)) squeeze(infoMergeAftBefFollower(3,5,:))];

%determine whether there is follower before merging
timeFollowerBefMergeMax = max(timeFollowerBefMerge,[],2);
indxWithFollowerBefMerge = find( timeFollowerBefMergeMax >= -timeTolFromMerge );
indxNoFollowerBefMerge = setdiff((1:numFollowerAftMergeTot)',indxWithFollowerBefMerge);
numWithFollowerBefMerge = length(indxWithFollowerBefMerge);

%calculate fraction of tracks, among those with follower after merging,
%that have follower before merging
fracFollowerBefAmongMergesWithFollowerAft = numWithFollowerBefMerge / numFollowerAftMergeTot;

%get oligomeric states before merging, for merges with or without follower
%before merging

%without follower
oligoStateNoFollower = oligoStateBefMerge(indxNoFollowerBefMerge,:);

%with follower
%first column is for segment with follower closer to the merge time
oligoStateWithFollower = ones(numWithFollowerBefMerge,2);
for iTmp = 1 : numWithFollowerBefMerge
    
    iMerge = indxWithFollowerBefMerge(iTmp);
    iCol = find( timeFollowerBefMerge(iMerge,:) == timeFollowerBefMergeMax(iMerge) );
    oligoStateWithFollower(iTmp,:) = [oligoStateBefMerge(iMerge,iCol(1)) oligoStateBefMerge(iMerge,3-iCol(1))]; %Sometimes, both segments have follower until the same point. In this case, put in first column the segment with more reliable follower presence. Thus the iCol(1).
    
end

%get fraction of monomers in each case
fracMonoNoFollowerBefMerge = sum(oligoStateNoFollower==1,1) ./ length(indxNoFollowerBefMerge); %Indicate dimension of sum too?
fracMonoWithFollowerBefMerge = sum(oligoStateWithFollower==1,1) ./ numWithFollowerBefMerge;

%get mean oligomeric state in each case
meanOligoNoFollowerBefMerge = mean(oligoStateNoFollower,1); %Indicate dimension of mean?
meanOligoWithFollowerBefMerge = mean(oligoStateWithFollower,1);


