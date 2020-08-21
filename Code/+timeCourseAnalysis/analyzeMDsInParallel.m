function [ redoMLanalysis ] = analyzeMDsInParallel( CMLs , doNewAnalysis, channels)
%ANALYZEMDSINPARALLEL Anlyze all the MovieData objects in a CombinedMovieList array
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

%Original: Probably Mark Kittisopikul
%Modified: K. Jaqaman, Feb. 2019

redoMLanalysis = false;

% Analyze all MovieData in parallel first
MLs = [CMLs.movieLists_];
MDs = [MLs.movies_];
numMDs = length(MDs);
saveFileArray = cell(1,numMDs);
channelArray{1} = channels;
channelArray = repmat(channelArray,1,numMDs);

% Get TimeCourseAnalysis processes or add them
procs = timeCourseAnalysis.getMovieDataTimeCourseAnalysisProcess(MDs);

if(doNewAnalysis)
    
    % Do per movie part of time course analysis (Alternatively, proc.run)
    MDsummaries = parcellfun_progress(@resultsIndTimeCoursePerMovie,MDs,saveFileArray,channelArray,'UniformOutput',false,'Heading','MovieData analysis');
    %     MDsummaries = cellfun(@resultsIndTimeCoursePerMovie,MDs,saveFileArray,channelArray,'UniformOutput',false); %keep for future debugging purposes to avoid parallelization problem
    
    % Assign results into processes
    [procs.summary_] = MDsummaries{:};
    cellfun(@save,MDs);
    redoMLanalysis = true;
    
else
    
    resSummary = arrayfun(@(proc) proc.summary_,procs,'UniformOutput',false);
    todo = cellfun('isempty',resSummary);
    
    if(any(todo))
        
        MDsummaries(todo) = parcellfun_progress(@resultsIndTimeCoursePerMovie,MDs(todo),saveFileArray(todo),channelArray(todo),'UniformOutput',false,'Heading','MovieData analysis');
        
        % Assign results into processes
        [procs(todo).summary_] = MDsummaries{todo};
        cellfun(@save,MDs);
        redoMLanalysis = true;
        
    end
    
end

if(redoMLanalysis)

    % Reset all MovieList level analysis if that needs to be redone
    % Potential efficiency gain: Figure out which ML needs to be reanalyzed
    % by mapping todo back to the ML
    procs = timeCourseAnalysis.getMovieObjectTimeCourseAnalysisProcess(MLs);
    [procs.summary_] = deal([]);

end

end

