%nested function for above: checking MD has necessary processes
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
function [] = MLCheck(ML, parameter)
    nMD = numel(ML.movies_);
    for iMD = 1:nMD
        if isempty(ML.movies_{iMD}.getProcessIndex('MotionAnalysisProcess'))
            error('timeCourseAnalysis:MotionAnalysisProcessMissing', ...
                ['MovieData ' ML.movies_{iMD}.movieDataFileName_ ...
                '\n at ' ML.movies_{iMD}.movieDataPath_ ...
                '\n does not contain MotionAnalysisProcess']);
        end
        if parameter.doPartition ...
                && isempty(ML.movies_{iMD}.getProcessIndex('PartitionAnalysisProcess'))
            error('timeCourseAnalysis:PartitionAnalysisProcessMissing', ... 
                ['MovieData ' ML.movies_{iMD}.movieDataFileName_ ...
                '\n at ' ML.movies_{iMD}.movieDataPath_ ...
                '\n does not contain PartitionAnalysisProcess']);
        end
    end
end

