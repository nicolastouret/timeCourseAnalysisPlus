function [ infoLink ] = info( clusterProfileName, ID )
%timeCourseAnalysis.link.info Generate link to find the job at a certain
%cluster and job ID
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

    infoLink = ['<a href="matlab: ' ...
        'findJob(parcluster(''' clusterProfileName '''),' ...
        '''ID'',' num2str(ID) ...
        '),' ...
        'disp(''To delete job, run: delete(ans)'');' ...
        'disp(timeCourseAnalysis.link.info(''' clusterProfileName ''',' num2str(ID) ' ));' ...
        'disp(timeCourseAnalysis.link.diary(''' clusterProfileName ''',' num2str(ID) ' ));' ...
        '">Click to Display Job Information</a>'];

end

