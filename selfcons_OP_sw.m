function [fESEs,fESExs]=selfcons_OP_sw(xs,sw,N,g,bc,ky)
%sdw=.14;
%soc=.06;
sdw=.001;
soc=.0001;
mu=1.1;
%ky=.1;
NE=N/2; % # of energies to sum over

H=pnict_5orb_Hamiltonian(N,ky,mu,xs,sw,sdw,soc,bc,g);
[Ud,Es]=eig(H);
Es=diag(sort(diag(Es),'ascend')); 
[c, ind]=sort(diag(Es),'ascend'); 
Ud=Ud(:,ind);

psi=zeros(N+bc,4,10,NE);

Fuidj_lt_a=zeros(N,N,10,NE);
Fdiuj_lt_a=zeros(N,N,10,NE);
Fuiuj_lt_a=zeros(N,N,10,NE);
Fdidj_lt_a=zeros(N,N,10,NE);
Fuidj_a=zeros(N,N,NE);
Fdiuj_a=zeros(N,N,NE);
Fuiuj_a=zeros(N,N,NE);
Fdidj_a=zeros(N,N,NE);
Fuidj=zeros(N,N);
Fdiuj=zeros(N,N);
Fuiuj=zeros(N,N);
Fdidj=zeros(N,N);
fESExs=zeros(1,N-1+bc);
fESEs=zeros(1,N);

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
    Fuidj_lt_a(i,j,lt,a) = psi(i,1,lt,a)*conj(psi(j,3,lt,a));   
	Fdiuj_lt_a(i,j,lt,a) = -psi(i,2,lt,a)*conj(psi(j,4,lt,a));
	Fuiuj_lt_a(i,j,lt,a) = -psi(i,1,lt,a)*conj(psi(j,4,lt,a));
    Fdidj_lt_a(i,j,lt,a) = psi(i,2,lt,a)*conj(psi(j,3,lt,a));   
                            end    
    Fuidj_a(i,j,a) =sum(Fuidj_lt_a(i,j,:,a));  
    Fdiuj_a(i,j,a) =sum(Fdiuj_lt_a(i,j,:,a)); 
    Fuiuj_a(i,j,a) =sum(Fuiuj_lt_a(i,j,:,a)); 
    Fdidj_a(i,j,a) =sum(Fdidj_lt_a(i,j,:,a)); 
                       end         
    Fuidj(i,j) =sum(Fuidj_a(i,j,:));  
    Fdiuj(i,j) =sum(Fdiuj_a(i,j,:)); 
    Fuiuj(i,j) =sum(Fuiuj_a(i,j,:)); 
    Fdidj(i,j) =sum(Fdidj_a(i,j,:)); 
                  end
              end
for i=1:N+bc-1
fESExs(i)=(Fuidj(i,i+1)+Fuidj(i+1,i)-Fdiuj(i,i+1)-Fdiuj(i+1,i))/4;
fESEs(i)=(Fuidj(i,i)-Fdiuj(i,i))/2;
end
fESEs(N)=(Fuidj(N,N)-Fdiuj(N,N))/2;
fESExs(N)=(Fuidj(N,N+1)+Fuidj(N+1,N)-Fdiuj(N,N+1)-Fdiuj(N+1,N))/4;
end