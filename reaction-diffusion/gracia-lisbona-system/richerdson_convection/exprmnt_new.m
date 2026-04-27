clc
clear all
format long e 
for ep_counter=1:8

    e2=2^(-2*(ep_counter));
    e1=e2*(2^(-2*(ep_counter-1)));
for N_counter =1:4
    
    
    Nx=64*2^(N_counter-1);
    Nt=2*4^(N_counter-1);
    %k=0.1*2^(-(N_counter-1));
    %Nt=1/k;
    k=1/Nt;
sig2=min(1/4,sqrt(e2)*log(Nx));
sig1=min(sig2/2,sqrt(e1)*log(Nx));
sig2;
sig1;
 
 x=zeros(1,Nx+1);
x_old=zeros(1,Nx+1);
U1_full=zeros(Nt+1,Nx+1);
U2_full=zeros(Nt+1,Nx+1);
U1=zeros(Nt+1,Nx+1);
err1=zeros(Nt+1,Nx+1);
err1_sup=zeros(1,Nx+1);
t=zeros(1,Nt+1);
u1_d=zeros(1,Nx+1);
u1_intl=zeros(1,Nx+1);
u1_sol=zeros(1,Nx+1);

U1_pre=zeros(Nt+1,Nx+1);

err2=zeros(Nt+1,Nx+1);
err2_sup=zeros(1,Nx+1);

u2_d=zeros(1,Nx+1);
u2_intl=zeros(1,Nx+1);
u2_sol=zeros(1,Nx+1);

U2_pre=zeros(Nt+1,Nx+1);


% construction of mesh points along  spatial variable and time variable
[x,h]=mesh_generation(sig1,sig2,Nx); 

for j=1:Nt+1
    t(j)=((j-1))/Nt;
end

%%%%%%%%%%%%%%%%%%%
% numerical approximation of initial and boundary value conditions and
% nonfomogeneous term
 for l=1:3
     if (l==2)
         Nx=2*Nx;
         Nt=2*Nt;
         k=1/Nt;
         sig2=min(1/4,sqrt(e2)*log(Nx/2));
         sig1=min(sig2/2,sqrt(e1)*log(Nx/2));
         [x,h]=mesh_generation(sig1,sig2,Nx);
  
  
         U1_full=zeros(Nt+1,Nx+1);
         U2_full=zeros(Nt+1,Nx+1);
         U1_sec=zeros(Nt+1,Nx+1);
         U2_sec=zeros(Nt+1,Nx+1);
   
         err1=zeros(Nt/2+1,Nx/2+1);  
         err2=zeros(Nt/2+1,Nx/2+1);%to match the dimension
            err_sup1=zeros(1,Nx/2+1);
            err_sup2=zeros(1,Nx/2+1);
            u1_d=zeros(1,Nx+1);
            u2_d=zeros(1,Nx+1);
            u1_intl=zeros(1,Nx+1);
             u2_intl=zeros(1,Nx+1);
            u1_sol=zeros(1,Nx+1);
            u2_sol=zeros(1,Nx+1);
     end

  if(l==3)
         Nx=2*Nx;
         Nt=2*Nt;
         k=1/Nt;
         sig2=min(1/4,sqrt(e2)*log(Nx/4));
         sig1=min(sig2/2,sqrt(e1)*log(Nx/4));
         [x,h]=mesh_generation(sig1,sig2,Nx);
  
  
         U1_full=zeros(Nt+1,Nx+1);
         U2_full=zeros(Nt+1,Nx+1);
%    U1_sec=zeros((T/del_t)+1,N+1);
%    U2_sec=zeros((T/del_t)+1,N+1);
   
         err1=zeros(Nt/2+1,Nx/2+1);  
         err2=zeros(Nt/2+1,Nx/2+1);%to match the dimension
            err_sup1=zeros(1,Nx/2+1);
            err_sup2=zeros(1,Nx/2+1);
            u1_d=zeros(1,Nx+1);
            u2_d=zeros(1,Nx+1);
             u1_intl=zeros(1,Nx+1);
            u2_intl=zeros(1,Nx+1);
            u1_sol=zeros(1,Nx+1);
            u2_sol=zeros(1,Nx+1);
  end
     
  time_counter=1;
r_count=1;
x;

while(time_counter<=Nt)
        if(time_counter==1)
            for j=1:Nx+1
                u1_d(1,j)=0;
             end
              U1_full(r_count,:)=u1_d;
                
              for j=1:Nx+1
                  u2_d(1,j)=0;
              end
               U2_full(r_count,:)=u2_d;
         else
                        clear u1_d;
                        clear u2_d;
                u1_d(1,:)=u1_intl(1,:);
                u2_d(1,:)=u2_intl(1,:);
        end
 
        [u1_sol,u2_sol]= solution(e1,e2,x,Nx,k,u1_d,u2_d,time_counter); 

                time_counter=time_counter+1;
                r_count=time_counter;
                time_counter;
                U1_full(r_count,:)=u1_sol(1,:);
                U2_full(r_count,:)=u2_sol(1,:);
                u1_intl(1,:)=u1_sol(1,:);
                u2_intl(1,:)=u2_sol(1,:);
    

end
if(l==1)
    U1_pre=U1_full;
    U2_pre=U2_full;
end
if(l==2)
   U1_extrp=(2*U1_full(1:2:Nt+1,1:2:Nx+1)-U1_pre);
   U1=U1_full(1:2:Nt+1,1:2:Nx+1);
   U1_sec=U1_full;
   
   U2_extrp=(2*U2_full(1:2:Nt+1,1:2:Nx+1)-U2_pre);
   U2=U2_full(1:2:Nt+1,1:2:Nx+1);
   U2_sec=U2_full;
 
end
if (l==3)
   U1_error=(2*U1_full(1:2:Nt+1,1:2:Nx+1)-U1_sec);
   err1_extrp=abs(U1_error(1:2:Nt/2+1,1:2:Nx/2+1)-U1_extrp);
   
   U2_error=(2*U2_full(1:2:Nt+1,1:2:Nx+1)-U2_sec);
   err2_extrp=abs(U2_error(1:2:Nt/2+1,1:2:Nx/2+1)-U2_extrp);
end

 end
Nx=Nx/4;
Nt=Nt/4;
k=4*k;

err1_normal=abs(U1-U1_pre);
err2_normal=abs(U2-U2_pre);

err1_sup_extrp=max(err1_extrp);
max_err1_extrp=max(err1_sup_extrp);
error1_extrp(ep_counter,N_counter)=max_err1_extrp;

err2_sup_extrp=max(err2_extrp);
max_err2_extrp=max(err2_sup_extrp);
error2_extrp(ep_counter,N_counter)=max_err2_extrp;


err1_sup_normal=max(err1_normal);
max_err1_normal=max(err1_sup_normal);
error1_normal(ep_counter,N_counter)=max_err1_normal;

err2_sup_normal=max(err2_normal);
max_err2_normal=max(err2_sup_normal);
error2_normal(ep_counter,N_counter)=max_err2_normal;
end
end

for j=1:ep_counter
    for i=1:N_counter-1
         conv1_extrp(j,i)=log2(error1_extrp(j,i)/error1_extrp(j,i+1));
         conv2_extrp(j,i)=log2(error2_extrp(j,i)/error2_extrp(j,i+1));
    end
end
for j=1:ep_counter
    for i=1:N_counter-1
         conv1_normal(j,i)=log2(error1_normal(j,i)/error1_normal(j,i+1));
         conv2_normal(j,i)=log2(error2_normal(j,i)/error2_normal(j,i+1));
    end
end
conv1_normal
error2_normal
conv2_normal

error1_extrp
conv1_extrp
error2_extrp
conv2_extrp
