clc
clear all
format short

for ep_counter=1:11
    %e=2^(-5*(ep_counter-1));
     e=10^-(ep_counter-1);
for N_t_counter=1:7
    N=2^(3+N_t_counter);

disp(['Row ' num2str(ep_counter) ' Column ' num2str(N_t_counter)]);

dis=(1/3);

%  sig_o=2.1;
% sig1=min(dis/2,sig_o*e*log(N));
% sig2=min((1-dis)/2,sig_o*e*log(N));

sig1=dis/2;
sig2=(1-dis)/2;


u_sol=zeros(1,N+1);

%[x h]=bsmesh(sig_o,sig1,sig2,dis,N,e);
[x h]=mesh(sig1,sig2,dis,N,e);               %Shishkin mesh

for err_counter=1:2
        if(err_counter==2)
            N=2*N;
%            sig1=min(dis/2,sig_o*e*log(N/2));
%             sig2=min((1-dis)/2,sig_o*e*log(N/2)); %
            sig1=dis/2;
             sig2=(1-dis)/2;

%              sig2=(1-dis)/2;
           % [x h]=bsmesh(sig_o,sig1,sig2,dis,N,e);
            [x h]=mesh(sig1,sig2,dis,N,e);                 %Shishkin mesh

            err=zeros(1,N/2+1);    %to match the dimension
            err_sup=zeros(1,N/2+1);

            u_sol=zeros(1,N+1);%to match the dimension 
        end

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
    max_err=max(err);
%     if(ep_counter==4 && N_t_counter==4)
%     figure(2);        
%         plot(x(1:2:N+1),err)
%     end
end

end


error(ep_counter,N_t_counter)=max_err;

end

end
error
for j=1:9
    for i=1:6
        conv(j,i)=log2(error(j,i)/error(j,i+1));
    end
end

conv


