function [xahomo,xs]=self_loop(it,N,g,bc,xs0,sw0,sdw,soc)
% fESEs,fOTEs,fOTuEs,fOTdEs
% xOPLm1,xOPLm2,xOPLm3,xOPLm4
% swAvg,swOP
% it,N = # of iterations,sites
% 
% g = interaction strength
w=.02; % Matsubara freqency
ky=.001;
%xs0=linspace(.04,.04,N).*(-.5*cos(ky));
%sw0=zeros(1,N);
%bc=1; % boundary (0=open, 1=periodic)

sw=zeros(it,N);
xs=zeros(it,N);
%swAvg=zeros(1,it+1);
xahomo=zeros(1,it+1);


xahomo(1)=sum(xs0)/length(xs0);
%swAvg(1)=sum(sw0)/length(sw0);
%[sw(1,:),xs(1,:)]=selfcons_OP_sw(xs0,sw0,N,g,bc,ky);
[sw(1,:),xs(1,:)]=self_OP(xs0,sw0,N,g,bc,ky,sdw,soc);
sw(1,:)=zeros(1,N);
%swAvg(2)=sum(sw(1,:))/length(sw(1,:));
xahomo(2)=sum(xs(1,:))/length(xs(1,:));
%swAvg(2)=sum(sw(1,10:N-10))/length(sw(1,10:N-10));

for i=2:it
%[sw(i,:),xs(i,:)]=selfcons_OP_sw(xs(i-1,:),sw(i-1,:),N,g,bc,ky);
[sw(i,:),xs(i,:)]=self_OP(xs(i-1,:),sw(i-1,:),N,g,bc,ky,sdw,soc);
sw(i,:)=zeros(1,N);
%swi=sw(i,1:N);
xsi=xs(i,1:N);
%swAvg(i+1)=sum(swi)/length(swi);
xahomo(i+1)=sum(xsi)/length(xsi);
%swAvg(i+1)=sum(sci(10:N-10))/length(sci(10:N-10));
end

[fESExs,fOTEs,fOTExs,fOTuEs,fOTdEs]=selfcons_OP_oddw(xs(it,:),sw(it,:),N,g,w,bc,ky,sdw,soc);

end
