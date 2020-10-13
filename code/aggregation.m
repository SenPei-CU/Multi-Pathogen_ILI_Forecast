function aggregation(region,season,ftime)
load signals
load mupost
load baseline
load bestpost
load('forecastens')
load('wposts')
z=bestpost(:,region);
%z=[ks;kc;k0;sigmaT;gamma;c];
pred=Aggregate(z,signals,forecastens,wposts,mupost(:,region),season,region,ftime);
%pred–row: num_times; column: num_ens
save('Aggregation.mat','pred')

%visualize forecast
ILI=signals(1:52,2,season,region);
hold on
plot(ftime:40,mean(pred(ftime:40,:),2),'LineWidth',2)
plot(1:40,ILI(1:40),'o','LineWidth',2)



function pred=Aggregate(z,signals,forecastens,wposts,mupost,season,region,ftime)
num_ens=size(wposts,2);
num_times=40;
ks=z(1);
kc=z(2);
k0=z(3);
sigmaT=z(4);
gamma=z(5);
c=z(6);
traj=zeros(num_times,num_ens,6);
%%%%%%%%%%%%%%%%%%%%%%latest week
for pid=1:6
    idx=randperm(num_ens);
    temp=signals(:,pid+2,season,region);
    traj(1:ftime-1,:,pid)=temp(1:ftime-1)*ones(1,num_ens);
    temp=forecastens{pid};
    traj(ftime:num_times,:,pid)=temp(:,idx);
end
w=wposts(:,randperm(num_ens));
for pid=1:6
    traj(:,:,pid)=(ones(num_times,1)*w(pid,:)).*traj(:,:,pid);
end

%calculate current discrepancy
%%%%%%%%%%%%%%%%%%%%%%%%%
obs=signals(ftime,2,season,region);
cntpred=sum(traj(ftime,:,:),3)+mupost(ftime)*ones(1,num_ens)+randn(1,num_ens).*abs(mupost(ftime)*ones(1,num_ens))/ks;

cntpred=mean(cntpred);
delta0=obs-cntpred;
%calculate season specific discrepancy
deltatf=delta0+randn(1,num_ens)*abs(delta0)/k0;
deltaT=randn(1,num_ens)*sigmaT;
deltat=zeros(num_times,num_ens);
for t=ftime+1:num_times
    deltat(t,:)=deltatf+(deltaT-deltatf)*((t-ftime)/(num_times-ftime))^gamma;
end

%raw aggregation
pred=sum(traj,3);
%add systematic bias
pred=pred+mupost(1:num_times)*ones(1,num_ens)+randn(num_times,num_ens).*abs(mupost(1:num_times)*ones(1,num_ens))/ks;
%add current bias
pred=pred+deltat+randn(num_times,num_ens).*abs(deltat)/kc;
%redistribute
pred=mean(pred,2)*ones(1,num_ens)+c*(pred-mean(pred,2)*ones(1,num_ens));

ILI=signals(1:ftime,2,season,region);
pred(1:ftime,:)=ILI*ones(1,num_ens);
pred=max(pred,0);