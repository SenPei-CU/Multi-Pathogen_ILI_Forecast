function modelfit(region,season,ftime)
%fit flu using MIF, output the best fit inital condition and parameters
%ftime: the forecast week should be between week 4 to week 35
%if flu activity is low, there is no fitting
seasons=[1997:2007,2010:2013];
regions={'National','Region 1','Region 2','Region 3','Region 4','Region 5',...
    'Region 6','Region 7','Region 8','Region 9'};
pathogens={'AH1','AH3','B','RSV','PIV12','PIV3'};
%15 seasons, 10 regions, 6 pathogens
load scale
load signals
startweek=4;
endweek=35;
ICs_flu=zeros(7,3);
for pid=1:3%fit flu
    [ic,~]=ICfit_flu(season,region,scale(region,pid),pid,ftime);
    ICs_flu(:,pid)=ic;
end
save('modelfit.mat','ICs_flu');
