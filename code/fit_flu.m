function [bestx,pred]=fit_flu(theta,obs,ts,num_times,AH)
%find initial condition (S, I) for flu, grid search
N=1e5;
Sl=0.3*N; Su=1*N;
Sstep=500;
Il=0; Iu=100;
Istep=(Iu-Il)/5;
num_S=(Su-Sl)/Sstep;
num_I=(Iu-Il)/Istep;

minMAE=1e5;
bestx=zeros(7,1);
for i=1:num_S
    for j=1:num_I
        x=[Sl+i*Sstep;Il+j*Istep+1;0;theta];
        MAE=getMAE(x,obs,ts,AH,num_times);
        if MAE<minMAE
            minMAE=MAE;
            bestx=x;
        end
    end
end
pred=plotfit(bestx,ts,AH);


function MAE=getMAE(x,obs,ts,AH,num_times)
N=1e5;dt=1;tmstep=7;discrete=0;
pred=zeros(num_times,1);
for t=1:num_times
    tcnt=ts+(t-1)*tmstep;
    x=SIRS_AH(x,tcnt,dt,tmstep,N,AH,discrete);
    pred(t)=x(3);
end
%MAE
MAE=sum(abs(obs(max(1,num_times-8):num_times)-pred(max(1,num_times-8):num_times)));

function pred=plotfit(x,ts,AH)
num_times=52;
N=1e5;dt=1;tmstep=7;discrete=0;
pred=zeros(num_times,1);
for t=1:num_times
    tcnt=ts+(t-1)*tmstep;
    x=SIRS_AH(x,tcnt,dt,tmstep,N,AH,discrete);
    pred(t)=x(3);
end
