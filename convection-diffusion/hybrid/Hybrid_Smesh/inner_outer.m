clc
clear all
format short e
%e=input('Give the perturbation coefficient  ');
C=1.001;
m=2;
T=input('Give a value of T');
delay_time=input('Give a value of Delaied time');
% inner=zeros(1,5);
% outer=zeros(1,5);
for ep_counter=4:4:8
%     e=exp(-ep_counter);
    e=10^(-ep_counter);
 for N_t_counter=4:5
   % N_t_counter=6;
    N=2^(3+N_t_counter);
    %N=input('Give a value of N');
disp(['Row ' num2str(ep_counter) ' Column ' num2str(N_t_counter)]);

% del_t=0.1*2^(-(N_t_counter-1));

%del_t=0.8/N;
M=N^2;
del_t=T/M;

p=delay_time/del_t;
gr_fact=T/delay_time;
% sig=input('Give the value of the fitting factor');
sig_o=4.2;
sig=min(.5,sig_o*e*log(N));
x=zeros(1,N+1);
U_full=zeros((T/del_t)+1,N+1);
U=zeros((T/del_t)+1,N+1);
err=zeros((T/del_t)+1,N+1);
err_sup=zeros(1,N+1);
t=zeros(1,(T/del_t)+1);
u_d=zeros(p+1,N+1);
u_intl=zeros(p+1,N+1);
u_sol=zeros(p+1,N+1);
w_sin=zeros(p+1,N+1);

r_count=0;
for i=1:(T/del_t)+1
    t(i)=(i-1)*del_t;
end
for i=1:N/2
   h(i)=2*(1-sig)/N;
end

for i=N/2+1:N
   h(i)=2*sig/N;
end

x(1)=0;
for i=1:N
   x(i+1)=x(i)+h(i);
end
 
flag=1;
counter=0;


time_counter=1;

        


while(time_counter<=gr_fact)
    
         if(time_counter==1)
            for n=1:p+1
                for j=1:N+1
                     k=n-1;
                     u_d(n,j)=exp(-(-k*del_t))*(exp(-1/e)+(1-exp(-1/e))*x(j)-exp(-(1-x(j))/e));
                end
            end

        else
            for n=1:p+1
                u_d(n,:)=u_intl(n,:);
            end
            u_d;
        end
        
    u_d;
    u_sol=solution(e,N,p,x,del_t,u_d,time_counter);
    
    
    

            for n=1:p+1
                U_full(r_count+n,:)=u_sol(n,:);
                u_intl(n,:)=u_sol(p+1-(n-1),:);
            end
 u_intl;
            r_count=time_counter*p;

    time_counter=time_counter+1;

end
x;
counter;
for j=1:gr_fact*p+1
    for i=1:N+1
%         U(gr_fact*p+1-(j-1),i)=exp(-((j-1)*del_t+x(i)/sqrt(e)));
          U(j,i)=exp(-((j-1)*del_t))*(exp(-1/e)+(1-exp(-1/e))*x(i)-exp(-(1-x(i))/e));
%           U(j,i)=exp(-((j-1)*del_t+x(i)/sqrt(e)));
    end
end
% for i=1:N+1
%     plot(P(i,:),0:counter-1);
% end

 err=abs(U-U_full);

% U_out=U(:,1:N/2+1);
% U_in=U(:,N/2+2:N);


%surf(x,t,U_full);

% surf(x,t,U)
if(ep_counter==4 ||ep_counter==8 )
    err1=err(:,1:N/2+1);
    err_sup1=max(err1);
    max_err1=max(err_sup1)
    err2=err(:,N/2+2:N+1);
    err_sup2=max(err2);
    max_err2=max(err_sup2)
%     inner(1,N_t_counter)=max_err1;
%     outer(1,N_t_counter)=max_err2;
end
 err_sup=max(err);
max_err=max(err_sup);
error(ep_counter,N_t_counter)=max_err;
 end
% if(ep_counter==4 ||ep_counter==8 )
% inner
% outer
% end

end
error
% for j=1:10
%     for i=1:4
%         conv(j,i)=log2(error(j,i)/error(j,i+1));
%     end
% end

conv
       
% hold off
counter;
