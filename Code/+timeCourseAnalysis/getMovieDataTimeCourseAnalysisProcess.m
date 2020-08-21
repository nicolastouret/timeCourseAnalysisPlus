function [ timeCourseAnalysisProcess ] = getMovieDataTimeCourseAnalysisProcess( MD , outputArray)
%getMovieDataTimeCourseAnalysisProcess Obtains a TimeCourseAnalysisProcess
%associated with a MovieData
% MD : MovieData , or an array of MovieData, or a cell array
% outputArray : true if an array of processes should be output (default)
%               or false if a cell array match the size of MD should be
%               output
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
    if(nargin < 2)
        outputArray = true;
    end

    % Deal with cell or arrays of MovieData
    if(iscell(MD))
        timeCourseAnalysisProcess = cellfun( ...
            @timeCourseAnalysis.getMovieDataTimeCourseAnalysisProcess, ...
            MD, ...
            'UniformOutput',false);
        if(outputArray)
            timeCourseAnalysisProcess = [timeCourseAnalysisProcess{:}];
        end
        return;
    elseif(~isscalar(MD))
        timeCourseAnalysisProcess = arrayfun( ...
            @timeCourseAnalysis.getMovieDataTimeCourseAnalysisProcess, ...
            MD, ...
            'UniformOutput',false);
        if(outputArray)
            timeCourseAnalysisProcess = [timeCourseAnalysisProcess{:}];
        end
        return;
    end

    % Check if TimeCourseAnalysisProcess already exists
    idx = MD.getProcessIndex('TimeCourseAnalysisProcess');
    if(isempty(idx))
        % If not, create a new one
        timeCourseAnalysisProcess = TimeCourseAnalysisProcess(MD);
        MD.addProcess(timeCourseAnalysisProcess);
    else
        timeCourseAnalysisProcess = MD.getProcess(idx);
    end
    timeCourseAnalysisProcess.sanityCheck();
end

