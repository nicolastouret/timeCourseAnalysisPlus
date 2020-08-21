function tracksReform = reformCompTracksRemoveArtifacts( tracksIn, timeThreshold)
%reformCompTracksRemoveArtifacts will call the function
%removeSplitMergeArtifactsChronologicalLuciana and reform the compTracks based in
%the output of this function
%
% SYNOPSIS [ tracksReform ] = reformCompTracksRemoveArtifacts( tracksIn )
%
% INPUT
%               tracks     : Output of trackCloseGapsKalman.
%
%       
%       timeThreshold:       Minimum time between a split and a re-merge or time
%                            after a segment that is borm to merge. To be use if
%                            removePotArtifacts is called.
%                            Default: 1.
%
%
% OUTPUT
%          tracksReform     : tracks after removing all artifacts for
%          merging and splitting
%
% Luciana de Oliveira, August 2017
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

% load compTracks
compTracks=tracksIn;
%initiate the compTrack
tracksReform=tracksIn;
% go over all compTracks
iTrackNew=length(compTracks)+1;


for iTrack = 1: length(compTracks) 
    % load seqOfEvents
   seqOfEventsIn = compTracks(iTrack).seqOfEvents;

%calculate length seqOfEvents

lengthSeqOfEvents= size(seqOfEventsIn,1);

% if the seqOfEvents have more than one segment

if lengthSeqOfEvents > 2
  %calculate number os segments
  tracksFeatIndxCG=compTracks(iTrack).tracksFeatIndxCG;
  numberOfSegments=size(tracksFeatIndxCG,1);
  
  
    % initiate indexRemove
  indexRemove=1;
  segRemove=0;
    indexRemoveRow=1;
  seqOfEventsArtifacts  = removeSplitMergeArtifactsChronologicalLuciana(seqOfEventsIn,timeThreshold);

    % Each row in the 
for segIndex=1:numberOfSegments
    
   
    % deterine if the segment is part of the compTrack
    condSegment=find(~isnan((seqOfEventsArtifacts(:,4))) & seqOfEventsArtifacts(:,3)==segIndex|seqOfEventsArtifacts(:,4)==segIndex, 1);
    
    % if the segment is not conected with any other it will be removed and
    % allocate as a new compTrack
    
    if isempty(condSegment)
         % save segmentNumber to remove from compTrack
    segRemove(indexRemove)=segIndex;
            % find where the segment starts and end
        startIndex= find( seqOfEventsArtifacts(:,2)==1 &  seqOfEventsArtifacts(:,3)==segIndex);
        endIndex= find( seqOfEventsArtifacts(:,2)==2 &  seqOfEventsArtifacts(:,3)==segIndex);
      
    % save the rows that need to be removed from seqOfEvents
    rowsRemove(indexRemoveRow)=startIndex;
    rowsRemove(indexRemoveRow+1)=endIndex;
        % add it to tracksReform as a individual track 
        tracksReform(iTrackNew).seqOfEvents=[seqOfEventsArtifacts(startIndex,:);seqOfEventsArtifacts(endIndex,:) ];
                
        %reallocate tracksFeatIndxCG and tracksCoordAmpCG 
    tracksReform(iTrackNew).tracksFeatIndxCG = compTracks(iTrack).tracksFeatIndxCG(segIndex,:);
    tracksReform(iTrackNew).tracksCoordAmpCG = compTracks(iTrack).tracksCoordAmpCG(segIndex,:);
   
       
    % if this track don't start in the first frame, remove the empty frames
    % in both tracksFeatIndxCG and tracksCoordAmpCG 
    %frame where the segment was borm
    timeSegStartSegment1=compTracks(iTrack).seqOfEvents(1,1);
    timeSegStartSegment2=tracksReform(iTrackNew).seqOfEvents(1,1);
    
    
    
    
    if timeSegStartSegment1~=timeSegStartSegment2
        
        % calculate the time difference
        timeDiff=timeSegStartSegment2-timeSegStartSegment1;
        tracksReform(iTrackNew).tracksFeatIndxCG = tracksReform(iTrackNew).tracksFeatIndxCG(:,timeDiff+1:end);
        tracksReform(iTrackNew).tracksCoordAmpCG = tracksReform(iTrackNew).tracksCoordAmpCG(:,8*(timeDiff)+1:end);
    end
    %increase index remove
    
   indexRemove=indexRemove+1;
   indexRemoveRow=indexRemoveRow+2;
    %increase iTrackNew
    iTrackNew=iTrackNew+1;
    end
end

if segRemove
 % remove this segment from compTrack
    tracksReform(iTrack).tracksFeatIndxCG(segRemove,:) =[];
    tracksReform(iTrack).tracksCoordAmpCG(segRemove,:) =[];
    tracksReform(iTrack).seqOfEvents(rowsRemove,:)=[];
     
   
%     
%     % remove this segments from seqOfEvents
%     tracksReform(iTrack).seqOfEvents(segRow,:)
    clear segRemove rowsRemove
end

end %if length seqOfEvents is greater than 2

end %for loop tracks
% some compTracks can be empty after the reform, this is step is to remove a compTrack if it is empty
tracksReform=tracksReform(~cellfun(@isempty,{tracksReform.seqOfEvents}));
% 
%    for i=1: length(tracksReform)
% seqOfEvents=tracksReform(i).seqOfEvents;
% if size (seqOfEvents,1)==2
% timeIni=seqOfEvents(1,1);
% timeEnd=seqOfEvents(2,1);
% timeSegment=timeEnd-timeIni;
% if timeSegment<=1
%   tracksReform(i)=[];  
% end
% end
% end
end

