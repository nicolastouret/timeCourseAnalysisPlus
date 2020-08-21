function [ ] = confirmDelete( clusterProfileName, ID )
%timeCourseAnalysis.link.confirmDelete Confirm job deletion before doing so
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

    if(strcmp(questdlg('Are you sure you want to delete the job?'),'Yes'))
        delete(findJob(parcluster(clusterProfileName),'ID',num2str(ID)))
        disp('Job Deleted');
    else
        disp('Job NOT deleted');
    end
end

