function [ic,pred]=ICfit_flu(season,region,scale,pid,ftime)
seasons=[1997:2007,2010:2013];
regions={'National','Region 1','Region 2','Region 3','Region 4','Region 5',...
    'Region 6','Region 7','Region 8','Region 9'};
pathogens={'AH1','AH3','B','RSV','PIV12','PIV3'};
%%%%%%%%%%%%%%%%%%%
load signals%national+10 regions
load AH
load xmax_flu
load xmin_flu
load oevpara
xmin=xmin_flu(:,region,pid);
xmax=xmax_flu(:,region,pid);
%%%%%%%%%%%%%%%%check if there is activity
threshold=0.01;
activity=0;
obs=signals(1:ftime,pid+2,season,region);
obs=obs-min(obs);
minidx=find(obs==min(obs));
minidx=minidx(end);
onsetweek=1;
for i=max(3,minidx):ftime
    if (obs(i-2)>=threshold)&&(obs(i-1)>=threshold)&&(obs(i)>=threshold)
        activity=1;
        onsetweek=i-2;
        break;
    end
end
timethreshold=[20,15,30];
if (activity==0)&&(ftime<=timethreshold(pid))%no activity
    ic=zeros(7,1);
    obs=signals(1:ftime,pid+2,season,region);
    pred=analogues(obs,signals,pid,season,region,scale);
else%there is activity
    backstep=4;
    AH=AH(:,region);
    obs=signals(max(1,onsetweek-backstep):ftime,pid+2,season,region)*scale;%starting from 4 weeks before onset
    ts=285+7*(max(1,onsetweek-backstep)-1);
    num_times=length(obs);
    oevbase=oevpara(1,pid,region)^2;
    oevfactor=oevpara(2,pid,region)^2;
    obs_var=zeros(1,num_times);
    for tt=1:num_times
        ave=obs(tt);
        obs_var(tt) = oevbase+ave^2/oevfactor;
    end
    %%%%%%%%%%%%%%%%%%%%
    num_ic=1;
    ic=zeros(7,num_ic);
    for id=1:num_ic
        num_ens=300;
        Iter=10;
        num_fit=1;
        paras=zeros(4,num_fit);
        for i=1:num_fit
            [para,~]=MIF_flu(obs,obs_var,ts,num_ens,Iter,num_times,AH,xmin,xmax);
            paras(:,i)=para;
        end
        theta=mean(paras,2);
        %%%%%%%%%%%%%%%%%%
        [ic(:,id),pred]=fit_flu(theta,obs,ts,num_times,AH);
        %pred starts from four weeks before onset
        temp=pred;
        pred=zeros(52,1);
        pred(max(1,onsetweek-backstep):52,1)=temp(1:52-max(1,onsetweek-backstep)+1);
        temp=signals(1:ftime,pid+2,season,region)*scale;
        pred(1:max(1,onsetweek-backstep)-1,1)=temp(1:max(1,onsetweek-backstep)-1,1);
    end
end


function pred=analogues(obs,signals,pid,season,region,scale)
if pid==1
    v=14; alpha=0.3; base=1e-3;
elseif pid==2
    v=14; alpha=0.4; base=1e-3;
elseif pid==3
    v=14; alpha=1; base=1e-3;
end
n=length(obs);
candidates=squeeze(signals(1:52,pid+2,:,region));
candidates(:,season)=[];%remove cnt season
dist=zeros(size(candidates,2),1);
for i=1:size(candidates,2)
    temp=candidates(1:n,i);
    diff=obs-temp;
    diff(diff==0)=base;
    dist(i)=sum(diff.^2);
end
dist(:,2)=1:size(candidates,2);
dist=sortrows(dist,1);
idx=dist(1:v,2);
pred=zeros(52,1);
for i=1:v
    w=dist(i,1)^(-alpha)/sum(dist(1:v,1).^(-alpha));
    pred=pred+w*candidates(:,idx(i));
end
pred(1:n)=obs;
pred=pred*scale;