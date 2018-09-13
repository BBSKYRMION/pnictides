function [fESEs,fOTEs,fOTuEs,fOTdEs]=selfcons_OP_oddw(xs,sc,N,g,w,bc,ky)
%sdw=.14;
%soc=.06;
sdw=.001;
soc=.0001;
mu=1.1;
%ky=.1;
NE=N; % # of energies to sum over
%xs=zeros(1,N);
H=pnict_5orb_Hamiltonian(N,ky,mu,xs,sc,sdw,soc,bc,g);
[Ud,Es]=eig(H);
Es=diag(sort(diag(Es),'ascend')); 
[c, ind]=sort(diag(Es),'ascend'); 
Ud=Ud(:,ind);

psi=zeros(N+bc,4,10,NE);

Fuidj_w_lt_a=zeros(N,N,10,NE);
Fdiuj_w_lt_a=zeros(N,N,10,NE);
Fuiuj_w_lt_a=zeros(N,N,10,NE);
Fdidj_w_lt_a=zeros(N,N,10,NE);
Fuidj_w_a=zeros(N,N,NE);
Fdiuj_w_a=zeros(N,N,NE);
Fuiuj_w_a=zeros(N,N,NE);
Fdidj_w_a=zeros(N,N,NE);
Fuidj_w=zeros(N,N);
Fdiuj_w=zeros(N,N);
Fuiuj_w=zeros(N,N);
Fdidj_w=zeros(N,N);
Fuidj_mw_lt_a=zeros(N,N,10,NE);
Fdiuj_mw_lt_a=zeros(N,N,10,NE);
Fuiuj_mw_lt_a=zeros(N,N,10,NE);
Fdidj_mw_lt_a=zeros(N,N,10,NE);
Fuidj_mw_a=zeros(N,N,NE);
Fdiuj_mw_a=zeros(N,N,NE);
Fuiuj_mw_a=zeros(N,N,NE);
Fdidj_mw_a=zeros(N,N,NE);
Fuidj_mw=zeros(N,N);
Fdiuj_mw=zeros(N,N);
Fuiuj_mw=zeros(N,N);
Fdidj_mw=zeros(N,N);
fESExs=zeros(1,N-1+bc);
fESEs=zeros(1,N);
fOTEs=zeros(1,N);

 
for lt=1:10 % lambda and tau dof
   for a=1:NE % energies
      for ms=1:4 % mu and sigma dof
         for n=1:N % sites
psi(n,ms,lt,a)=Ud((n-1)*40+(ms-1)*10+lt,a);
         end
             if bc==1
             psi(N+1,ms,lt,a)=psi(1,ms,lt,a);
             else
             end
       end
    end
end
              for i=1:N+bc
                  for j=1:N+bc
                      for a=1:NE    
                            for lt=1:10
    Fuidj_w_lt_a(i,j,lt,a) = psi(i,1,lt,a)*conj(psi(j,3,lt,a))/(1i*w-Es(a));   
	Fdiuj_w_lt_a(i,j,lt,a) = -psi(i,2,lt,a)*conj(psi(j,4,lt,a))/(1i*w-Es(a));
	Fuiuj_w_lt_a(i,j,lt,a) = -psi(i,1,lt,a)*conj(psi(j,4,lt,a))/(1i*w-Es(a));
    Fdidj_w_lt_a(i,j,lt,a) = psi(i,2,lt,a)*conj(psi(j,3,lt,a))/(1i*w-Es(a));  
    Fuidj_mw_lt_a(i,j,lt,a) = psi(i,1,lt,a)*conj(psi(j,3,lt,a))/(-1i*w-Es(a));   
	Fdiuj_mw_lt_a(i,j,lt,a) = -psi(i,2,lt,a)*conj(psi(j,4,lt,a))/(-1i*w-Es(a));
	Fuiuj_mw_lt_a(i,j,lt,a) = -psi(i,1,lt,a)*conj(psi(j,4,lt,a))/(-1i*w-Es(a));
    Fdidj_mw_lt_a(i,j,lt,a) = psi(i,2,lt,a)*conj(psi(j,3,lt,a))/(-1i*w-Es(a));  
                            end    
    Fuidj_w_a(i,j,a) =sum(Fuidj_w_lt_a(i,j,:,a));  
    Fdiuj_w_a(i,j,a) =sum(Fdiuj_w_lt_a(i,j,:,a)); 
    Fuiuj_w_a(i,j,a) =sum(Fuiuj_w_lt_a(i,j,:,a)); 
    Fdidj_w_a(i,j,a) =sum(Fdidj_w_lt_a(i,j,:,a)); 
    Fuidj_mw_a(i,j,a) =sum(Fuidj_mw_lt_a(i,j,:,a));  
    Fdiuj_mw_a(i,j,a) =sum(Fdiuj_mw_lt_a(i,j,:,a)); 
    Fuiuj_mw_a(i,j,a) =sum(Fuiuj_mw_lt_a(i,j,:,a)); 
    Fdidj_mw_a(i,j,a) =sum(Fdidj_mw_lt_a(i,j,:,a));  
                       end         
    Fuidj_w(i,j) =sum(Fuidj_w_a(i,j,:));  
    Fdiuj_w(i,j) =sum(Fdiuj_w_a(i,j,:)); 
    Fuiuj_w(i,j) =sum(Fuiuj_w_a(i,j,:)); 
    Fdidj_w(i,j) =sum(Fdidj_w_a(i,j,:)); 
    Fuidj_mw(i,j) =sum(Fuidj_mw_a(i,j,:));  
    Fdiuj_mw(i,j) =sum(Fdiuj_mw_a(i,j,:)); 
    Fuiuj_mw(i,j) =sum(Fuiuj_mw_a(i,j,:)); 
    Fdidj_mw(i,j) =sum(Fdidj_mw_a(i,j,:)); 
                  end
              end
for i=1:N+bc-1
fESExs(i)=(Fuidj_w(i,i+1)+Fuidj_w(i+1,i)-Fdiuj_w(i,i+1)-Fdiuj_w(i+1,i))/4;
fESEs(i)=(Fuidj_w(i,i)-Fdiuj_w(i,i)+Fuidj_mw(i,i)-Fdiuj_mw(i,i))/4;
fOTEs(i)=(Fuidj_w(i,i)+Fdiuj_w(i,i)-Fuidj_mw(i,i)-Fdiuj_mw(i,i))/4;
end
fESEs(N)=(Fuidj_w(N,N)-Fdiuj_w(N,N)+Fuidj_mw(N,N)-Fdiuj_mw(N,N))/4;
fOTEs(N)=(Fuidj_w(N,N)+Fdiuj_w(N,N)-Fuidj_mw(N,N)-Fdiuj_mw(N,N))/4;
fOTuEs(N)=(Fuiuj_w(N,N)-Fuiuj_mw(N,N))/2;
fOTdEs(N)=(Fdidj_w(N,N)-Fdidj_mw(N,N))/2;
end