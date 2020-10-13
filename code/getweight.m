function getweight(region,season,ftime)
load signals
load oevpara
load offset
oevpara(1,:,:)=offset;
seasons=[1997:2007,2010:2013];
regions={'National','Region 1','Region 2','Region 3','Region 4','Region 5',...
    'Region 6','Region 7','Region 8','Region 9'};
pathogens={'AH1','AH3','B','RSV','PIV12','PIV3'};

num_sample=1000;

ILI=signals(1:ftime,2,season,region);
v=signals(1:ftime,3:8,season,region);
wposts=MCMC(oevpara(:,:,region),v,ILI,num_sample);
        
save('wposts.mat','wposts');

function wpost=MCMC(oevpara,v,ILI,num_sample)
%smooth v
ftime=size(v,1);
vs=zeros(ftime,6);
vs(1,:)=v(1,:);
vs(ftime,:)=v(ftime,:);
for i=2:ftime-1
    vs(i,:)=(v(i-1,:)+v(i,:)+v(i+1,:))/3;
end
%calculate OEV
oevs=zeros(ftime,6);
for i=1:ftime
    for pid=1:6
        oevs(i,pid)=oevpara(1,pid)^2+vs(i,pid)^2/oevpara(2,pid)^2;
    end
end
wcnt=1*rand(6,1);
Iter=2e4;
wrec=zeros(6,Iter);
cnt=0;
loglikelihoodcnt=calculatell(wcnt,ILI,vs,oevs);
llrec=zeros(Iter,1);
arec=zeros(Iter,1);
while cnt<Iter
    cnt=cnt+1;
    %update each w
    for i=1:6
        %update w(i)
        wnew=wcnt;
        wnew(i)=wnew(i)+0.01*randn();
        wnew(i)=max(0.01,wnew(i));
        wnew(i)=min(1,wnew(i));
        loglikelihoodnew=calculatell(wnew,ILI,vs,oevs);
        a=exp(loglikelihoodnew-loglikelihoodcnt);
        if (a>1)||(rand()<a)
            wcnt=wnew;
            loglikelihoodcnt=loglikelihoodnew;
        end
        
    end
    wrec(:,cnt)=wcnt;
    llrec(cnt)=loglikelihoodcnt;
    arec(cnt)=a;
end
burnin=0.5;
wrec=wrec(:,round(size(wrec,2)*burnin):size(wrec,2));
skip=10;
wrec=wrec(:,1:skip:size(wrec,2));
wpost=wrec(:,(size(wrec,2)-num_sample+1:size(wrec,2)));



function loglikelihood=calculatell(wcnt,ILI,vs,oevs)
%calculate aggregated distribution
mu=vs*wcnt;
oev=oevs*(wcnt.^2);
sigma=sqrt(oev);
loglikelihood=0;
for t=1:size(ILI,1)
    loglikelihood=loglikelihood+log(normpdf(ILI(t),mu(t),sigma(t))*0.01);
end

