clc
clear all
%e=input('Give the perturbation coefficient  ');
C=1.1;
m=2;
T=1;

for ep_counter=8:8
    e2=2^(-(ep_counter-1));
    e1=e2*(2^(-2*(ep_counter-1)));
for N_t_counter= 1:2
    N=64*2^(N_t_counter-1);
    %N=input('Give a value of N');
%disp(['Row ' num2str(ep_counter) ' Column ' num2str(N_t_counter)]);
del_x=1/N;
del_t=0.1*2^(-(N_t_counter-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%
x=zeros(1,N+1);
x_old=zeros(1,N+1);
Y=zeros(1,N+1);
h=zeros(1,N);
U1_full=zeros((T/del_t)+1,N+1);
V1_full=zeros((T/del_t)+1,N+1);
U2_full=zeros((T/del_t)+1,N+1);
V2_full=zeros((T/del_t)+1,N+1);
U=zeros((T/del_t)+1,N+1);
err=zeros((T/del_t)+1,N+1);
err_sup=zeros(1,N+1);
t=zeros(1,(T/del_t)+1);
u1_d=zeros(1,N+1);
v1_d=zeros(1,N+1);
u1_intl=zeros(1,N+1);
v1_intl=zeros(1,N+1);
u1_sol=zeros(1,N+1);
v1_smo=zeros(1,N+1);

u2_d=zeros(1,N+1);
v2_d=zeros(1,N+1);
u2_intl=zeros(1,N+1);
v2_intl=zeros(1,N+1);
u2_sol=zeros(1,N+1);
v2_smo=zeros(1,N+1);
w1_sin=zeros((T/del_t)+1,N+1);
w2_sin=zeros((T/del_t)+1,N+1);

u1_double_d=zeros(1,2*N+1);
u2_double_d=zeros(1,2*N+1);
u1_double_intl=zeros(1,2*N+1);
u2_double_intl=zeros(1,2*N+1);
U1_double_full=zeros((T/del_t)+1,2*N+1);
U2_double_full=zeros((T/del_t)+1,2*N+1);
U1=zeros(1,2*N+1);
U2=zeros(1,2*N+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%
r_count=1;
for i=1:N+1
    x(i)=del_x*(i-1);
end
for i=2:N+1
    j=i-1;
    h(j)=x(i)-x(i-1);
end

%flag=1;
counter=0;

time_counter=1;

    
 while(time_counter<=(T/del_t)+1)
          for l=1:2
              if(l==1)
                flag=1;
                while(flag==1)    
                    if(time_counter==1)
                         for j=1:N+1
                             u1_d(1,j)=0;
                         end
                         U1_full(r_count,:)=u1_d;
                
                        for j=1:N+1
                         u2_d(1,j)=0;
                        end
                        U2_full(r_count,:)=u2_d;
                
                
                        for j=1:N+1
                            v1_d(1,j)=0;
                        end
                        V1_full(r_count,:)=v1_d;
                
                        for j=1:N+1
                          v2_d(1,j)=0;
                        end
                        V2_full(r_count,:)=v2_d;
                    else
                        clear u1_d;
                        clear u2_d;
                        u1_d(1,:)=u1_intl(1,:);
                        v1_d(1,:)=v1_intl(1,:);
                        u2_d(1,:)=u2_intl(1,:);
                        v2_d(1,:)=v2_intl(1,:);

                    end
                    [u1_sol,u2_sol]=solution(e1,e2,x,N,del_t,u1_d,u2_d,time_counter);
                    [v1_smo,v2_smo]=smooth(e1,e2,x,N,del_t,v1_d,v2_d,time_counter);
            

                    W1=u1_sol-v1_smo;
                    W2=u2_sol-v2_smo;
   
%                     w1_sin=u1_sol-v1_smo;
%                     w2_sin=u2_sol-v2_smo;

%                     W1=w1_sin(round((T/del_t)),:);
%                     W2=w2_sin(round((T/del_t)),:);
    
                    for i=2:N+1
                         j=i-1;
                         h(j)=x(i)-x(i-1);
                    end 
                    a_dis=0;
       
                     a_dis=h(1)*((abs(2*(((W1(3)-W1(2))/h(2))-((W1(2)-W1(1))/h(1)))/(h(2)+h(1))))^(1/2)+(abs(2*(((W2(3)-W2(2))/h(2))-((W2(2)-W2(1))/h(1)))/(h(2)+h(1))))^(1/2)); %h(i) has been changed for adjusting index
       
                    a_dis=a_dis+h(N)*((abs(2*(((W1(N+1)-W1(N))/h(N))-((W1(N)-W1(N-1))/h(N-1)))/(h(N)+h(N-1))))^(1/2)+(abs(2*(((W2(N+1)-W2(N))/h(N))-((W2(N)-W2(N-1))/h(N-1)))/(h(N)+h(N-1))))^(1/2)); %h(i) has been changed for adjusting index
       
                    for i=3:N
                        a_dis=a_dis+h(i-1)*(((abs(2*(((W1(i)-W1(i-1))/h(i-1))-((W1(i-1)-W1(i-2))/h(i-2)))/(h(i-1)+h(i-2))))^(1/2)+(abs(2*(((W1(i+1)-W1(i))/h(i))-((W1(i)-W1(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2))/2+((abs(2*(((W2(i)-W2(i-1))/h(i-1))-((W2(i-1)-W2(i-2))/h(i-2)))/(h(i-1)+h(i-2))))^(1/2)+(abs(2*(((W2(i+1)-W2(i))/h(i))-((W2(i)-W2(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2))/2); %h(i) has been changed for adjusting index
                    end
                    M=zeros(N+1,1);
                    M(1)=M(2);
                    M(N+1)=M(N);
                    a_dis;
                    for i=2:N
                      M(i)=a_dis+(abs(2*(((W1(i+1)-W1(i))/h(i))-((W1(i)-W1(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2)+(abs(2*(((W2(i+1)-W2(i))/h(i))-((W2(i)-W2(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2);
                    end
                    H=zeros(N,1);
                    for i=2:N+1
                     j=i-1;
                     H(j)=((M(i-1)+M(i))/2)*h(j);           
                    end
                     L=zeros(N+1,1);
                     L(1)=0;
                     for i=1:N
                       sum=0;
                          for j=1:i
                          sum=sum+H(j);
                          end
                          L(i+1)=sum;
                     end
                     C_tol=N*max(H)/L(N+1);
                     if(C_tol<=C)
                        flag=0;
                        x_old=x;
                     else
                        for i=1:N+1
                            Y(i)=(i-1)*L(N+1)/N;
                        end
                        new_x=interp1(L,x,Y);
                        x=new_x;
                     end
        %%%%%%%%%%%%%%%%%%%%%
%                      F1=U1;
%                      F2=U2;
                end
      %  time_counter=1;
                
                
              else
                 for i=1:N+1
                    x_double(2*i-1)=x(i);
                 end
                 for i=1:N
                     x_double(2*i)=x(i)+0.5*(x(i+1)-x(i));
                 end
                 P=2*N;
                  
                 for i=2:P+1
                    j=i-1;
                    h_double(j)=x_double(i)-x_double(i-1);
                 end 
                 if(time_counter==1)
                     for j=1:P+1
                        u1_double_d(1,j)=0;
                     end
                     U1_double_full(r_count,:)=u1_double_d;
                     for j=1:P+1
                        u2_double_d(1,j)=0;
                     end
                     U2_double_full(r_count,:)=u2_double_d;
                 else
                        clear u1_double_d;
                        clear u2_double_d;
                        u1_double_d(1,:)=u1_double_intl(1,:);
                        u2_double_d(1,:)=u2_double_intl(1,:);

                 end
                 [U1,U2]=solution(e1,e2,x_double,P,del_t/2,u1_double_d,u2_double_d,time_counter); 
                 G1=U1;
                 G2=U2;
              end
          end
                u1_intl(1,:)=u1_sol(1,:);
                v1_intl(1,:)=v1_smo(1,:);
                u2_intl(1,:)=u2_sol(1,:);
                v2_intl(1,:)=v2_smo(1,:);
                time_counter=time_counter+1;
                r_count=time_counter;
                time_counter;
                U1_full(r_count,:)=u1_sol(1,:);
                V1_full(r_count,:)=v1_smo(1,:);
                U2_full(r_count,:)=u2_sol(1,:);
                V2_full(r_count,:)=v2_smo(1,:);
                
                u1_double_intl(1,:)=U1;
                u2_double_intl(1,:)=U2;
                U1_double_full(r_count,:)=U1;
                U2_double_full(r_count,:)=U2;
 end

 
% counter;
% for j=1:gr_fact*p+1
%     for i=1:N+1
% %         U(gr_fact*p+1-(j-1),i)=exp(-((j-1)*del_t+x(i)/sqrt(e)));
%           U(j,i)=exp(-((j-1)*del_t))*(exp(-1/e)+(1-exp(-1/e))*x(i)-exp(-(1-x(i))/e));
%     end
% end
% 
% U_full;
% U;
%  err=abs(U-U_full);
% 
% 
% 
% %surf(x,t,U_full);
% 
% % surf(x,t,U)
% 
%  err_sup=max(err);
% max_err=max(err_sup);
% error(ep_counter,N_t_counter)=max_err;
end

end
error
for j=1:10
    for i=1:5
        conv(j,i)=log2(error(j,i)/error(j,i+1));
    end
end
conv