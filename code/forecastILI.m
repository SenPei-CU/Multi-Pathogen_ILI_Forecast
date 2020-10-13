function forecastILI(regionid,seasonid,ftime)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code fore paper: Sen Pei & Jeffrey Shaman, Aggregating forecasts of
% multiple respiratory pathogens supports more accurate forecasting of
% influenza-like illness, PLoS Computational Biology. 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%regionid, e.g., regionid=1 means National
%seasonid, e.g., seasonid=1 means 1997
%ftime, forecast week, from 4 to 35
%see input details below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seasons=[1997:2007,2010:2013];
% regions={'National','Region 1','Region 2','Region 3','Region 4','Region 5',...
%     'Region 6','Region 7','Region 8','Region 9'};
% pathogens={'AH1','AH3','B','RSV','PIV12','PIV3'};
% ftime: forecast week - between week 4 to week 35
% week 1 = EPIWEEK 40, week 2 = EPIWEEK 41, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%output
%Aggregation: ILI prediction
%forecastens: forecasts for six pathogens
%wposts: weighting factor for each pathogen
%modelfit: initial condition for flu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit circulating flu strain
% output best fitting initial condition and parameters
modelfit(regionid,seasonid,ftime);
%generate ensemble predictions for each pathogen
prepareensemble(regionid,seasonid,ftime);
%estimate multiplicative factors using MCMC
getweight(regionid,seasonid,ftime);
%aggregate and post-process
aggregation(regionid,seasonid,ftime);

