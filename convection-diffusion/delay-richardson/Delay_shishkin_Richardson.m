clc
clear all
%e=input('Give the perturbation coefficient  ');

T=input('Give a value of T');
delay_time=input('Give a value of Delaied time');

for ep_counter=1:8
%     e=2^(-2*(ep_counter-1));
% e=exp(-ep_counter);
 e=10^(-ep_counter);
for N_t_counter= 1:6
    N=2^(3+N_t_counter);
    %N=input('Give a value of N');
disp(['Row ' num2str(ep_counter) ' Column ' num2str(N_t_counter)]);

del_t=0.1*2^(-(N_t_counter-1));
%  del_t=1/N;
p=delay_time/del_t;
gr_fact=T/delay_time;
% sig=input('Give the value of the fitting factor');

sig_o=2.2;
% alpha=0.4672;
% sig=min(.5,2*e*log(N)/alpha);
sig=min(.5,sig_o*e*log(N));
% sig=2.2*e*log(N);
% sig=.5;
x=zeros(1,N+1);
x_old=zeros(1,N+1);
U_full=zeros((T/del_t)+1,N+1);
U=zeros((T/del_t)+1,N+1);
err=zeros((T/del_t)+1,N+1);
err_sup=zeros(1,N+1);
t=zeros(1,(T/del_t)+1);
u_d=zeros(p+1,N+1);
u_intl=zeros(p+1,N+1);
u_sol=zeros(p+1,N+1);
w_sin=zeros(p+1,N+1);
U_pre=zeros((T/del_t)+1,N+1);

r_count=0;
for i=1:(T/del_t)+1
    t(i)=(i-1)*del_t;
end
% for i=1:N+1
%     if (i<=N/2+1)    % adjusted for indexing
%         x(i)=(i-1)*2*(1-sig)/N;
%     else
%         x(i)=1-sig+((i-1)-N/2)*2*sig/N;
%     end
%    
% end
% x
% for i=1:N+1
%     if(i<=N/2+1)
%         h(i)=2*(1-sig)/N;
%     else
%         h(i)=2*sig/N;
%     end
% end

for i=1:N/2
   h(i)=2*(1-sig)/N;
end

H1=2*(1-sig)/N;


for i=N/2+1:N
   h(i)=2*sig/N;
end

H2=2*sig/N;

x(1)=0;
for i=1:N
   x(i+1)=x(i)+h(i);
end
x(N+1)=1; 

x_old=x; %storing for use it in exact solution
counter=0;


time_counter=1;

        

for err_counter=1:3
        if(err_counter==2)
            N=2*N;
            del_t=del_t/2;
%             del_t=1/N;
            p=delay_time/del_t;
            h=zeros(1,N);
            x=zeros(1,N+1);
            sig=min(.5,sig_o*e*log(N/2)); %
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
            x(N+1)=1;
            
            U_full=zeros((T/del_t)+1,N+1);
            U_sec=zeros((T/del_t)+1,N+1);
%             U=zeros((T/del_t)+1,N+1);
            err=zeros((T/(2*del_t))+1,N/2+1);    %to match the dimension
            err_sup=zeros(1,N/2+1);
            u_d=zeros(p+1,N+1);
            u_intl=zeros(p+1,N+1);
            u_sol=zeros(p+1,N+1);%to match the dimension 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%             h=zeros(1,2*N);
%             x=zeros(1,2*N+1);
%             for i=1:N
%                 h(i)=(H1)/2;
%             end
% 
%             for i=N+1:2*N
%                 h(i)=(H2)/2;
%             end
% 
%             x(1)=0;
%             for i=1:2*N
%                 x(i+1)=x(i)+h(i);
%             end
%             x(2*N+1)=1;
%             
%             U_full=zeros((2*T/del_t)+1,2*N+1);
% %             U=zeros((T/del_t)+1,N+1);
%             err=zeros((T/(del_t))+1,N+1);    %to match the dimension
%             err_sup=zeros(1,N+1);
%             p=2*p;
%             N=2*N;
%             del_t=del_t/2;
%             u_d=zeros(p+1,N+1);
%             u_intl=zeros(p+1,N+1);
%             u_sol=zeros(p+1,N+1);%to match 
        end
        if(err_counter==3)
           N=2*N;
            del_t=del_t/2;
            p=delay_time/del_t;
            h=zeros(1,N);
            x=zeros(1,N+1);
            sig=min(.5,sig_o*e*log(N/4)); %
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
            x(N+1)=1;
            
            U_full=zeros((T/del_t)+1,N+1);
            
%             U=zeros((T/del_t)+1,N+1);
            err=zeros((T/(2*del_t))+1,N/2+1);    %to match the dimension
            err_sup=zeros(1,N/2+1);
            u_d=zeros(p+1,N+1);
            u_intl=zeros(p+1,N+1);
            u_sol=zeros(p+1,N+1);%to match the dimension  
        end
time_counter=1;
r_count=0;
x;

while(time_counter<=gr_fact)
    
    x;
    
         if(time_counter==1)
            u_d=zeros(p+1,N+1);
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
if(err_counter==1)
    U_pre=U_full;
    x;
    
%     T
%     N
%     p
%     del_t
end
if(err_counter==2)
        
%     T
%     N
%     p
%     del_t
    x(1:2:N+1);
   U_extrp=(2*U_full(1:2:(T/del_t)+1,1:2:N+1)-U_pre);
   U=U_full(1:2:(T/del_t)+1,1:2:N+1);
   U_sec=U_full;
%       %%%%%%%%%%%%%%%for surface plot%%%%%%%%%%%%%%%%%%%%
%    for i=1:(T/del_t)+1
%     t(i)=(i-1)*del_t;
%    end
%   
%    if(ep_counter==8 && N_t_counter==3)
%    figure(1);   
%    surf(x(1:2:N+1),t(1:2:(T/del_t)+1),U_extrp);
%    figure(2);
%    surf(x(1:2:N+1),t(1:2:(T/del_t)+1),U_pre); 
%    end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
end
if(err_counter==3)
   U_error=(2*U_full(1:2:(T/del_t)+1,1:2:N+1)-U_sec);
   err_extrp=abs(U_error(1:2:T/(2*del_t)+1,1:2:N/2+1)-U_extrp);

end
end
p=p/4;
N=N/4;
del_t=4*del_t;
x;
counter;
% for j=1:gr_fact*p+1
%     for i=1:N+1
%         k=j-1;
%           U(j,i)=exp(-(k*del_t))*(exp(-1/e)+(1-exp(-1/e))*x_old(i)-exp(-(1-x_old(i))/e));
%     end
% end

U_full;
U;
% err_extrp=abs(U_error-U_extrp);
err_normal=abs(U-U_pre);


err_sup_extrp=max(err_extrp);
max_err_extrp=max(err_sup_extrp);
error_extrp(ep_counter,N_t_counter)=max_err_extrp;


err_sup_normal=max(err_normal);
max_err_normal=max(err_sup_normal);
error_normal(ep_counter,N_t_counter)=max_err_normal;

end

end
error_extrp
for j=1:5
    for i=1:4
        conv_extrp(j,i)=log2(error_extrp(j,i)/error_extrp(j,i+1));
%         conv_extrp(j,i)=log(error_extrp(j,i)/error_extrp(j,i+1))/log((2*log(2^(3+i)))/log(2*(2^(3+i))));
    end
end
conv_extrp

error_normal
for j=1:5
    for i=1:4
        conv_normal(j,i)=log2(error_normal(j,i)/error_normal(j,i+1));
%        conv_normal(j,i)=log(error_normal(j,i)/error_normal(j,i+1))/log((2*log(2^(3+i)))/log(2*(2^(3+i))));
    end
end

conv_normal
       
% hold off
counter;
