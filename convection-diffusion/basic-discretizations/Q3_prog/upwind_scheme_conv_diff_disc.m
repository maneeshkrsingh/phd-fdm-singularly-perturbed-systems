clc
clear all
format short
% %e=input('Give the perturbation coefficient  ');
% T=input('Give a value of T');
% delay_time=input('Give a value of Delaied time');
% hold on
for ep_counter=1:9
  %   e=2^(-5*(ep_counter-1));
     e=2^(-2*(ep_counter-1)); 
    %e=10^(-ep_counter);
for N_t_counter=1:7
    N=2^(3+N_t_counter);
    %N=input('Give a value of N');
disp(['Row ' num2str(ep_counter) ' Column ' num2str(N_t_counter)]);

% del_t=0.1*2^(-(N_t_counter-1));

% sig=input('Give the value of the fitting factor');
dis=0.5;

  sig_o1=2.2;
  sig_o2=2.2;
sig1=min(dis/2,sig_o1*e*log(N));
sig2=min((1-dis)/2,sig_o2*e*log(N));
 %sig1=dis/2;
 %sig2=(1-dis)/2;
U_full=zeros(1,N+1);
U=zeros(1,N+1);

u_sol=zeros(1,N+1);
%[x h]=bsmesh(sig_o1,sig_o2,sig1,sig2,dis,N,e);  %BS mesh
[x h]=mesh(sig1,sig2,dis,N,e);         %Shishkin mesh
flag=1;
counter=0;



for err_counter=1:2
        if(err_counter==2)
            N=2*N;
           sig1=min(dis/2,sig_o1*e*log(N/2));
            sig2=min((1-dis)/2,sig_o2*e*log(N/2)); %
             %sig1=dis/2;
             %sig2=(1-dis)/2;
           % [x h]=bsmesh(sig_o1,sig_o2,sig1,sig2,dis,N,e);         %BS mesh
            [x h]=mesh(sig1,sig2,dis,N,e);        %Shishkin mesh
            

            err=zeros(1,N/2+1);    %to match the dimension
            err_sup=zeros(1,N/2+1);

            u_sol=zeros(1,N+1);%to match the dimension 
        end
time_counter=1;
r_count=0;
x;

    

    u_sol=solution(e,N,x);
    if(err_counter==1 && ep_counter==3 && N_t_counter==4)
   % figure(1);
        plot(x,u_sol)
        hold on
    end
    
    if(err_counter==1 && ep_counter==7 && N_t_counter==4)
   % figure(1);
        plot(x,u_sol,'*')
        hold off
    end



if(err_counter==1)
    U_pre=u_sol;
else
    err=abs(u_sol(1,1:2:N+1)-U_pre);
 %   err_sup=max(err);
    max_err=max(err);
%     if(ep_counter==4 && N_t_counter==4)
%     figure(2);        
%         plot(x(1:2:N+1),err)
%     end
end

end
U_full;
x;
counter;

% if(N_t_counter==1)
%     U_pre=U_full;
%     x
% else
%     x
% %     U_post=U_full;
%     err=abs(U_full(1:2:(T/del_t)+1,1:2:N+1)-U_pre);
%      err_sup=max(err);
%     max_err=max(err_sup);
%     error(ep_counter,N_t_counter-1)=max_err;
%     U_pre=U_full;
% end
error(ep_counter,N_t_counter)=max_err;

end
% clear Matrix U_pre;
% clear Matrix err
% if ep_counter==2
% plot(x,U_full(321,:),'r')
% end
% if ep_counter==3
% plot(x,U_full(321,:),'g')
% end
% if ep_counter==4
% plot(x,U_full(321,:))
% end
end
error
% hold off
for j=1:9
    for i=1:6
        conv(j,i)=log2(error(j,i)/error(j,i+1));
    end
end

conv
%        
% hold off
counter;
