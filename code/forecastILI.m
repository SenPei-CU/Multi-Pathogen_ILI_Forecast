function forecastILI(region,season,ftime)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code fore paper: Sen Pei & Jeffrey Shaman, Aggregating forecasts of
% multiple respiratory pathogens supports more accurate forecasting of
% influenza-like illness, PLoS Computational Biology. 2020.
% Programmed by Sen Pei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seasons=[1997:2007,2010:2013];
% regions={'National','Region 1','Region 2','Region 3','Region 4','Region 5',...
%     'Region 6','Region 7','Region 8','Region 9'};
% pathogens={'AH1','AH3','B','RSV','PIV12','PIV3'};
% ftime: forecast week - between week 4 to week 35
% week 1 = EPIWEEK 40, week 2 = EPIWEEK 41, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit circulating flu strain
% output best fitting initial condition and parameters
modelfit(region,season,ftime);
%generate ensemble predictions for each pathogen
prepareensemble(region,season,ftime);
%estimate multiplicative factors using MCMC
getweight(region,season,ftime);
%aggregate and post-process
aggregation(region,season,ftime);

