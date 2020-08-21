%% Time Course Analysis (MD-level)
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
function [MDExtra] = MDAnalyze(MD) %#ok<INUSD>
    %Need blank if not used
    MDExtra.blank = [];
    %{
    %loads 'partitionResult'
    if analysisPara.doPartition
        load(MD.processes_{channel, MD.getProcessIndex('PartitionAnalysisProcess')}.outFilePaths_{1});
    end
    %loads 'tracks' and 'diffAnalysisRes'
    if analysisPara.doPartition
        load(MD.processes_{channel, MD.getProcessIndex('MotionAnalysisProcess')}.outFilePaths_{1});
    end
    %}
end
