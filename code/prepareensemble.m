function prepareensemble(region,season,ftime)
%generate ensemble predictions for each pathogen
seasons=[1997:2007,2010:2013];
regions={'National','Region 1','Region 2','Region 3','Region 4','Region 5',...
    'Region 6','Region 7','Region 8','Region 9'};
pathogens={'AH1','AH3','B','RSV','PIV12','PIV3'};
%15 seasons, 10 regions, 6 pathogens
load scale
load signals
load modelfit
load oevpara
load optimalptb%for MIF
load bestptb%for analogues
load AH
load xmin_flu
load xmax_flu
xmin=xmin_flu;
xmax=xmax_flu;

num_times=40;
num_ens=1000;
forecastens=cell(6,1);
for pid=1:6
    forecastens{pid}=zeros(num_times-ftime+1,num_ens);
end

%for flu
for pid=1:3
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
    if (activity==1)%there is onset
        %MIF
        IC=ICs_flu(:,pid);
        z=optimalptb(:,pid,region);
        sigma0=z(1);
        kp=z(2);
        sigmay=oevpara(1,pid,region);
        ky=oevpara(2,pid,region);
        backstep=4;
        Tstart=max(1,onsetweek-backstep);
        pred=perturb_flu(Tstart,IC,ftime,num_times-ftime,AH(:,region),scale(region,pid),num_ens,sigma0,kp,xmin(:,region,pid),xmax(:,region,pid));
        temp=pred(ftime:num_times,:);
        temp=temp+randn(size(temp,1),num_ens).*(sqrt(sigmay^2+temp.^2/ky^2));
        temp=max(temp,0);
        forecastens{pid}=temp;
    else
        %analogues
        obs=signals(1:ftime,pid+2,season,region);
        ptb=bestptb(pid,region);
        pred=analogues_flu(obs,signals,pid,season,region,ptb,num_ens);
        temp=pred(ftime:num_times,:);
        forecastens{pid}=temp;
    end
end

% for nonflu
for pid=4:6
    obs=signals(1:ftime,pid+2,season,region);
    ptb=bestptb(pid,region);
    pred=analogues_nonflu(obs,signals,pid,season,region,ptb,num_ens);
    temp=pred(ftime:num_times,:);
    forecastens{pid}=temp;
end

save('forecastens.mat','forecastens');


function pred=perturb_flu(Tstart,IC,ftime,num_times,AH,scale,num_ens,sigma0,kp,xmin,xmax)
ts=285+7*(Tstart-1);
N=1e5;dt=1;tmstep=7;discrete=0;
pred=zeros(52,num_ens);
%intergrate from Tstart to ftime
for t=1:ftime-Tstart
    tcnt=ts+(t-1)*tmstep;
    IC=SIRS_AH(IC,tcnt,dt,tmstep,N,AH,discrete);
    pred(Tstart+t,:)=IC(3,:)/scale;
end

%generate perturbed IC
%find eigenvector
b=log(IC(4)-IC(5)); a=-180;
AHtemp=[AH;AH];
BT1=exp(a*AHtemp+b)+IC(5);
BETA=BT1/IC(7);
L=IC(6);D=IC(7);
S=IC(1);I=IC(2);
ts=tcnt+tmstep;
sigmaS=S;
sigmaI=I;

M=[-2*(1/L+BETA(ts)*I/N),BETA(ts)*I/N*sigmaS/sigmaI-(1/L+BETA(ts)*S/N)*sigmaI/sigmaS;BETA(ts)*I/N*sigmaS/sigmaI-(1/L+BETA(ts)*S/N)*sigmaI/sigmaS,2*(BETA(ts)*S/N-1/D)];
[E,Lambda]=eig(M);
if Lambda(1,1)>Lambda(2,2)
    e1=E(:,1);
else
    e1=E(:,2);
end
if e1(1)<0
    e1=-e1;
end
lambda1=max(Lambda(1,1),Lambda(2,2));


ptb=sigma0*exp(-kp*lambda1);


ptl=0.025:0.95/(num_ens-1):0.975;
rnd = norminv(ptl,0,ptb);

x=IC*ones(1,num_ens);

x(1,:)=x(1,:)+sigmaS*e1(1)*rnd;
x(2,:)=x(2,:)+sigmaI*e1(2)*rnd;

x=checkbound(x,xmin,xmax);
ts=tcnt+tmstep;
for t=1:num_times
    tcnt=ts+(t-1)*tmstep;
    x=SIRS_AH(x,tcnt,dt,tmstep,N,AH,discrete);
    pred(ftime+t,:)=x(3,:)/scale;
end

function x = checkbound(x,xmin,xmax)
N=1e5;
%S
% x(1,x(1,:)<xmin(1))=xmin(1);
% x(1,x(1,:)>xmax(1))=xmax(1);
x(1,x(1,:)<0)=0;
x(1,x(1,:)>N)=N;
%I
x(2,x(2,:)<0)=0;
x(2,x(2,:)>N)=N;
%obs
x(3,x(3,:)<0)=mean(x(3,:));
x(3,x(3,:)>N)=median(x(3,:));
%if R0max < R0min
x(4,x(4,:)<xmin(4))=xmin(4);
x(4,x(4,:)>xmax(4))=xmax(4);
x(5,x(5,:)<xmin(5))=xmin(5);
x(5,x(5,:)>xmax(5))=xmax(5);
x(4,x(4,:)<x(5,:))=x(5,x(4,:)<x(5,:))+0.01;
%L
x(6,x(6,:)<xmin(6))=xmin(6);
x(6,x(6,:)>xmax(6))=xmax(6);
%D
x(7,x(7,:)<xmin(7))=xmin(7);
x(7,x(7,:)>xmax(7))=xmax(7);

function pred=analogues_flu(obs,signals,pid,season,region,ptb,num_ens)
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
w=zeros(v,1);
predave=zeros(52,1);
for i=1:v
    w(i)=dist(i,1)^(-alpha)/sum(dist(1:v,1).^(-alpha));
    predave=predave+w(i)*candidates(:,idx(i));
end
y = randsample(idx,num_ens,true,w);

sample=candidates(:,y);

sample=sample.*(1+randn(52,num_ens)*0.05);

total=sum(sample)';total(:,2)=(1:num_ens)';
total=sortrows(total,1);
index=total(:,2)';

sample=sample(:,index);

pred=predave*ones(1,num_ens)+ptb*(sample-predave*ones(1,num_ens));
pred=max(pred,0);

function pred=analogues_nonflu(obs,signals,pid,season,region,ptb,num_ens)
if pid==4%RSV
    v=14; alpha=1; base=1e-3;
elseif pid==5
    v=14; alpha=0.6; base=1e-3;
elseif pid==6
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
w=zeros(v,1);
predave=zeros(52,1);
for i=1:v
    w(i)=dist(i,1)^(-alpha)/sum(dist(1:v,1).^(-alpha));
    predave=predave+w(i)*candidates(:,idx(i));
end
y = randsample(idx,num_ens,true,w);

sample=candidates(:,y);

sample=sample.*(1+randn(52,num_ens)*0.05);

total=sum(sample)';total(:,2)=(1:num_ens)';
total=sortrows(total,1);
index=total(:,2)';

sample=sample(:,index);

pred=predave*ones(1,num_ens)+ptb*(sample-predave*ones(1,num_ens));
pred=max(pred,0);