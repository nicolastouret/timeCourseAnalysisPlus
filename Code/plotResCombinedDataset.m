function plotResCombinedDataset(varPerCond,timePerCond,namePerCond,...
    colorPerCond,figNameList,dir2save,yaxisUnits,shiftNegTime)

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

numCond = length(namePerCond);
numFig = length(figNameList);

%check whether input includes standard deviations
if size(varPerCond{1},3) == 1
    stdFlag = 0;
else
    stdFlag = 1;
end

%% Plot

%plot each variable separately
for iVar = 1 : numFig
    
    h = figure('Name',figNameList{iVar}); hold on
    
    time0Info = NaN(numCond,2);
    
    for iCond = 1 : numCond
        
        %get this condition's time points
        time4plot = timePerCond{iCond};
        
        %shift time to start at 0, store index of original 0
        if shiftNegTime
            indx0 = find(time4plot==0);
            if indx0 > 1
                time4plot = time4plot - time4plot(1);
                time0Info(iCond,:) = [time4plot(indx0) varPerCond{iCond}(indx0,iVar,1)];
            end            
        end
        
        %make plot
        plot(time4plot,varPerCond{iCond}(:,iVar,1),'Color',colorPerCond{iCond});
        if stdFlag
            myErrorbar(time4plot,varPerCond{iCond}(:,iVar,1),varPerCond{iCond}(:,iVar,2)./sqrt(varPerCond{iCond}(:,iVar,3)));
        end
        
    end
    
    legend(namePerCond)
    
    %indicate original 0 if time was shifted
    if shiftNegTime
        plot(time0Info(:,1),time0Info(:,2),'ko','MarkerSize',10)
    end
    
    xlabel('Time (min)')
    ylabel([figNameList{iVar} ' ' yaxisUnits])
    
    savefig(h,fullfile(dir2save,figNameList{iVar}));
    
end