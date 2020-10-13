function [para,theta]=MIF_flu(obstruth,obs_var,ts,num_ens,Iter,num_times,AH,xmin,xmax)
dt=1;
tmstep=7;
discrete=0;
%%%initialization
N=1e5;
%%%%%%%MIF setting
num_var=size(xmin,1);%number of variables
theta=zeros(4,Iter+1);%parameters at each iteration

sig=zeros(1,Iter);
alp=0.9;%alpha
SIG=(xmax-xmin).^2/4/num_times;%variance of theta

%%%start iteration for Iter round
for n=1:Iter
    xprior=zeros(num_var,num_ens,num_times);
    xpost=zeros(num_var,num_ens,num_times);
    sig(n)=alp^(n-1);
    So=zeros(num_var,num_ens);
    %%%generate new ensemble members using multivariate normal distribution
    if (n==1)
        Sigma=diag(sig(n)^2*SIG);
    else 
        Sigma=diag(sig(n)^2*SIG);
    end
    if (n==1)
        %first guess of state space
        x0=lhsu(xmin,xmax,num_ens);
        So=x0;
        theta(:,1)=mean(x0(4:7,:),2);
    else
        x0=lhsu(xmin,xmax,num_ens);
        So(1:3,:)=x0(1:3,:);%initialize variables
        So(4:7,:)=mvnrnd(theta(:,n)',Sigma(4:7,4:7),num_ens)';%parameters
        So(3,:)=0;%observations
    end
    %%%correct lower/upper bounds of the parameters
    So=checkbound(So,xmin,xmax);
    %start looping
    for t=1:num_times
        tcnt=ts+(t-1)*tmstep;
        %integrate the model
        So=SIRS_AH(So,tcnt,dt,tmstep,N,AH,discrete);
        xprior(:,:,t)=So;
        obs=So(3,:);
        %EAKF
        prior_var = var(obs);
        post_var = prior_var*obs_var(t)/(prior_var+obs_var(t));
        if prior_var==0
            post_var=0;
            prior_var=1e-3;
        end
        prior_mean = mean(obs);
        post_mean = post_var*(prior_mean/prior_var + obstruth(t)/obs_var(t));
        alpha = (obs_var(t)/(obs_var(t)+prior_var)).^0.5;
        dy = post_mean + alpha*((obs)-prior_mean)-obs;
        rr=zeros(1,size(So,1));
        for j=1:size(So,1)
            A=cov(So(j,:),obs);
            rr(j)=A(2,1)/prior_var;
        end
        dx=rr'*dy;
        So=So+dx;
        So=checkbound(So,xmin,xmax);
        xpost(:,:,t) = So;
    end
    temp=squeeze(mean(xpost(4:7,:,:),2));
    theta(:,n+1)=mean(temp,2);
end

para=theta(:,Iter+1);



function s=lhsu(xmin,xmax,nsample)
% LHS from uniform distribution
nvar=length(xmin);
ran=rand(nsample,nvar);
s=zeros(nsample,nvar);
for j=1: nvar
   idx=randperm(nsample);
   P =(idx'-ran(:,j))/nsample;
   s(:,j) = xmin(j) + P.* (xmax(j)-xmin(j));
end
s=s';

function x = checkbound(x,xmin,xmax)
N=1e5;
%S
x(1,x(1,:)<xmin(1))=xmin(1);
x(1,x(1,:)>xmax(1))=xmax(1);
%I
x(2,x(2,:)<0)=mean(x(2,:));
x(2,x(2,:)>N)=median(x(2,:));
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

