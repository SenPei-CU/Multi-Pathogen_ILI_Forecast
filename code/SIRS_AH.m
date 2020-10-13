function x=SIRS_AH(x,ts,dt,tmstep,N,AH,discrete)
%x=[S,I,obs,R0max,R0min,L,D]
if discrete==0
AH(length(AH)+1:2*length(AH))=AH;
num_ens=size(x,2);
BT1=zeros(length(AH),num_ens);
BETA=zeros(length(AH),num_ens);
for i=1:num_ens
  b=log(x(4,i)-x(5,i)); a=-180;
  BT1(:,i)=exp(a*AH+b)+x(5,i);
  BETA(:,i)=BT1(:,i)/x(7,i);
end
BETA=[BETA;BETA];
Sr=zeros(3,num_ens,tmstep+1);
Incidence=zeros(tmstep+1,num_ens);
Sr(:,:,1)=x(1:3,:);
L=x(6,:);
D=x(7,:);
%start integration
tcnt=0;
for t=ts+dt:dt:ts+tmstep*dt
    tcnt=tcnt+1;
    betat=floor(t);
    Eimmloss=dt*((N*ones(1,num_ens)-Sr(1,:,tcnt)-Sr(2,:,tcnt))./L);
    Einf=min(dt*(BETA(betat,:).*Sr(1,:,tcnt).*Sr(2,:,tcnt)/N),Sr(1,:,tcnt));
    Erecov=min(dt*(Sr(2,:,tcnt)./D),Sr(2,:,tcnt));
    Eimmloss=max(Eimmloss,0);Einf=max(Einf,0);Erecov=max(Erecov,0);
    sk1=Eimmloss-Einf;
    ik1=Einf-Erecov;
    ik1i=Einf;
    Ts1=Sr(1,:,tcnt)+(sk1/2);
    Ti1=Sr(2,:,tcnt)+(ik1/2);
    Eimmloss=dt*((N*ones(1,num_ens)-Ts1-Ti1)./L);
    Einf=min(dt*(BETA(betat,:).*Ts1.*Ti1/N),Ts1);
    Erecov=min(dt*(Ti1./D),Ti1);
    Eimmloss=max(Eimmloss,0);Einf=max(Einf,0);Erecov=max(Erecov,0);
    sk2=Eimmloss-Einf;
    ik2=Einf-Erecov;
    ik2i=Einf;
    Ts2=Sr(1,:,tcnt)+(sk2/2);
    Ti2=Sr(2,:,tcnt)+(ik2/2);
    Eimmloss=dt*((N*ones(1,num_ens)-Ts2-Ti2)./L);
    Einf=min(dt*(BETA(betat,:).*Ts2.*Ti2/N),Ts2);
    Erecov=min(dt*(Ti2./D),Ti2);
    Eimmloss=max(Eimmloss,0);Einf=max(Einf,0);Erecov=max(Erecov,0);
    sk3=Eimmloss-Einf;
    ik3=Einf-Erecov;
    ik3i=Einf;
    Ts3=Sr(1,:,tcnt)+(sk3);
    Ti3=Sr(2,:,tcnt)+(ik3);
    Eimmloss=dt*((N*ones(1,num_ens)-Ts3-Ti3)./L);
    Einf=min(dt*(BETA(betat,:).*Ts3.*Ti3/N),Ts3);
    Erecov=min(dt*(Ti3./D),Ti3);
    Eimmloss=max(Eimmloss,0);Einf=max(Einf,0);Erecov=max(Erecov,0);
    sk4=Eimmloss-Einf;
    ik4=Einf-Erecov;
    ik4i=Einf;
    Sr(1,:,tcnt+1)=Sr(1,:,tcnt)+(sk1/6+sk2/3+sk3/3+sk4/6);
    Sr(2,:,tcnt+1)=Sr(2,:,tcnt)+(ik1/6+ik2/3+ik3/3+ik4/6);
    Incidence(tcnt+1,:)=(ik1i/6+ik2i/3+ik3i/3+ik4i/6);
end
x(1,:)=Sr(1,:,tcnt+1);
x(2,:)=Sr(2,:,tcnt+1);
x(3,:)=sum(Incidence)/N;
end

if discrete==1

AH(length(AH)+1:2*length(AH))=AH;
num_ens=size(x,2);
BT1=zeros(length(AH),num_ens);
BETA=zeros(length(AH),num_ens);
for i=1:num_ens
  b=log(x(4,i)-x(5,i)); a=-180;
  BT1(:,i)=exp(a*AH+b)+x(5,i);
  BETA(:,i)=BT1(:,i)/x(7,i);
end
BETA=[BETA;BETA];
Sr=zeros(3,num_ens,tmstep+1);
Incidence=zeros(tmstep+1,num_ens);
Sr(:,:,1)=x(1:3,:);
L=x(6,:);
D=x(7,:);
%start integration
tcnt=0;
for t=ts+dt:dt:ts+tmstep*dt
    tcnt=tcnt+1;
    betat=floor(t);
    Eimmloss=dt*((N*ones(1,num_ens)-Sr(1,:,tcnt)-Sr(2,:,tcnt))./L);
    Einf=min(dt*(BETA(betat,:).*Sr(1,:,tcnt).*Sr(2,:,tcnt)/N),Sr(1,:,tcnt));
    Erecov=min(dt*(Sr(2,:,tcnt)./D),Sr(2,:,tcnt));
    Eimmloss=max(Eimmloss,0);Einf=max(Einf,0);Erecov=max(Erecov,0);
    l=poissrnd([Eimmloss;Einf;Erecov]);
    sk1=l(1,:)-l(2,:);
    ik1=l(2,:)-l(3,:);
    ik1i=l(2,:);
    Ts1=Sr(1,:,tcnt)+(sk1/2);
    Ti1=Sr(2,:,tcnt)+(ik1/2);
    Eimmloss=dt*((N*ones(1,num_ens)-Ts1-Ti1)./L);
    Einf=min(dt*(BETA(betat,:).*Ts1.*Ti1/N),Ts1);
    Erecov=min(dt*(Ti1./D),Ti1);
    Eimmloss=max(Eimmloss,0);Einf=max(Einf,0);Erecov=max(Erecov,0);
    l=poissrnd([Eimmloss;Einf;Erecov]);
    sk2=l(1,:)-l(2,:);
    ik2=l(2,:)-l(3,:);
    ik2i=l(2,:);
    Ts2=Sr(1,:,tcnt)+(sk2/2);
    Ti2=Sr(2,:,tcnt)+(ik2/2);
    Eimmloss=dt*((N*ones(1,num_ens)-Ts2-Ti2)./L);
    Einf=min(dt*(BETA(betat,:).*Ts2.*Ti2/N),Ts2);
    Erecov=min(dt*(Ti2./D),Ti2);
    Eimmloss=max(Eimmloss,0);Einf=max(Einf,0);Erecov=max(Erecov,0);
    l=poissrnd([Eimmloss;Einf;Erecov]);
    sk3=l(1,:)-l(2,:);
    ik3=l(2,:)-l(3,:);
    ik3i=l(2,:);
    Ts3=Sr(1,:,tcnt)+(sk3);
    Ti3=Sr(2,:,tcnt)+(ik3);
    Eimmloss=dt*((N*ones(1,num_ens)-Ts3-Ti3)./L);
    Einf=min(dt*(BETA(betat,:).*Ts3.*Ti3/N),Ts3);
    Erecov=min(dt*(Ti3./D),Ti3);
    Eimmloss=max(Eimmloss,0);Einf=max(Einf,0);Erecov=max(Erecov,0);
    l=poissrnd([Eimmloss;Einf;Erecov]);
    sk4=l(1,:)-l(2,:);
    ik4=l(2,:)-l(3,:);
    ik4i=l(2,:);
    Sr(1,:,tcnt+1)=Sr(1,:,tcnt)+(sk1/6+sk2/3+sk3/3+sk4/6);
    Sr(2,:,tcnt+1)=Sr(2,:,tcnt)+(ik1/6+ik2/3+ik3/3+ik4/6);
    Incidence(tcnt+1,:)=(ik1i/6+ik2i/3+ik3i/3+ik4i/6);
end
x(1,:)=Sr(1,:,tcnt+1);
x(2,:)=Sr(2,:,tcnt+1);
x(3,:)=sum(Incidence)/N;
end    
    
end
