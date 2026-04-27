clc
clear all
format short
e=input('Give the perturbation coefficient  ');
C=1.001;
m=2;
N=input('Give a value of N');
del_x=1/N;
T=input('Give a value of T');
delay_time=input('Give a value of Delaied time');
del_t=0.1;
p=delay_time/del_t;
gr_fact=T/delay_time;
r_count=0;
for i=1:N+1
    x(i)=del_x*(i-1);
end
for i=2:N+1
    j=i-1;
    h(j)=x(i)-x(i-1);
end
w_sin=zeros(p+1,N+1);
flag=1;
counter=0;


time_counter=1;
while(time_counter<=gr_fact)
    
    while(flag==1)
%         P(:,counter+1)=x;
        counter=counter+1;
        counter;
        if(time_counter==1)
            for n=1:p+1
                for j=1:N+1
                     k=n-1;
                     u_d(n,j)=exp(-(-k*del_t))*(exp(-1/e)+(1-exp(-1/e))*x(j)-exp(-(1-x(j))/e));
                end
            end
            for n=1:p+1
                for j=1:N+1
                    k=n-1;
                    v_d(n,j)=exp(-(-k*del_t))*(exp(-1/e)+(1-exp(-1/e))*x(j)-exp(-(1-x(j))/e));
                end
            end
        else
            for n=1:p+1
                u_d(n,:)=interp1(x_old,u_intl(n,:),x);
                v_d(n,:)=interp1(x_old,v_intl(n,:),x);
            end
        end
        u_d;
        v_d;
        u_sol=solution(e,N,p,x,del_t,u_d,time_counter);
        v_smo=smooth(e,N,p,x,del_t,v_d,time_counter);
        u_sol;
        v_smo;
        w_sin=u_sol-v_smo;

        w=w_sin(round(p-1),:);
        
        for i=2:N+1
            j=i-1;
            h(j)=x(i)-x(i-1);
        end
        h;
        a_dis=0;
       
        a_dis=h(1)*(abs(2*(((w(3)-w(2))/h(2))-((w(2)-w(1))/h(1)))/(h(2)+h(1))))^(1/2); %h(i) has been changed for adjusting index
       
        a_dis=a_dis+h(N)*(abs(2*(((w(N+1)-w(N))/h(N))-((w(N)-w(N-1))/h(N-1)))/(h(N)+h(N-1))))^(1/2); %h(i) has been changed for adjusting index
       
        for i=3:N
            a_dis=a_dis+h(i-1)*((((abs(2*(((w(i)-w(i-1))/h(i-1))-((w(i-1)-w(i-2))/h(i-2)))/(h(i-1)+h(i-2))))^(1/2)+(abs(2*(((w(i+1)-w(i))/h(i))-((w(i)-w(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2)))/2); %h(i) has been changed for adjusting index
        end
        M=zeros(N+1,1);
        a_dis;
        for i=2:N
            M(i)=a_dis+(abs(2*(((w(i+1)-w(i))/h(i))-((w(i)-w(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/2);
        end

        M(1)=M(2);
        M(N+1)=M(N);
        M;

        for i=2:N+1
            j=i-1;
            H(j)=((M(i-1)+M(i))/2)*h(j);
        end

        H;
        h;
        L(1)=0;
        for i=1:N
            sum=0;
            for j=1:i
                sum=sum+H(j);
            end
            L(i+1)=sum;
        end

        C_tol=N*max(H)/L(N+1);
        C_tol;
        C;
        x;
            
        if(C_tol<=C)
            flag=0;
            for n=1:p+1
                U_full(r_count+n,:)=u_sol(n,:);
                u_intl(n,:)=u_sol(p+1-(n-1),:);
                v_intl(n,:)=v_smo(p+1-(n-1),:);
            end
            u_sol;
            u_intl;
            x_old=x;
            r_count=time_counter*p;
        else
            for i=1:N+1
                Y(i)=(i-1)*L(N+1)/N;
            end
            L;
            x;
            Y;
            new_x=interp1(L,x,Y);
            x=new_x;
        end
        u_sol;

        x(1)=0;
       % x(N+1)=1;
    

%U;
    end
    time_counter=time_counter+1;
    flag=1;
    
end
x;
% hold off
% u
% U
% x
% 
% err=abs(U-u)
% err_sup=max(err)
% P
%  hold on
% P;
counter;
for j=1:gr_fact*p+1
    for i=1:N+1
%         U(gr_fact*p+1-(j-1),i)=exp(-((j-1)*del_t+x(i)/sqrt(e)));
          U(j,i)=exp(-((j-1)*del_t))*(exp(-1/e)+(1-exp(-1/e))*x(i)-exp(-(1-x(i))/e));
    end
end
% for i=1:N+1
%     plot(P(i,:),0:counter-1);
% end
x
U_full
U
 err=abs(U-U_full)
err_sup=max(err);
max_err=max(err_sup)
% hold off
counter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc
% clear all
% format long
% e=input('Give the perturbation coefficient  ');
% C=1.001;
% m=2;
% N=input('Give a value of N');
% del_x=1/N;
% T=input('Give a value of T');
% delay_time=input('Give a value of Delaied time');
% del_t=0.1;
% p=delay_time/del_t;
% gr_fact=T/delay_time;
% r_count=0;
% 
% for i=1:(T/del_t)+1
%     t(i)=(i-1)*del_t;
% end
% 
% for i=1:N+1
%     x(i)=del_x*(i-1);
% end
% for i=2:N+1
%     j=i-1;
%     h(j)=x(i)-x(i-1);
% end
% w_sin=zeros(p+1,N+1);
% flag=1;
% counter=0;
% 
% 
% time_counter=1;
% while(time_counter<=gr_fact)
%     
%     while(flag==1)
% %         P(:,counter+1)=x;
%         counter=counter+1;
%         counter;
%         if(time_counter==1)
%             for n=1:p+1
%                 for j=1:N+1
%                      k=n-1;
%                      u_d(n,j)=exp(-(-k*del_t))*(exp(-1/e)+(1-exp(-1/e))*x(j)-exp(-(1-x(j))/e));
%                     % u_d(n,j)=exp(-(-k*del_t+x(j)/sqrt(e)));
%                 end
%             end
%             for n=1:p+1
%                 for j=1:N+1
%                     k=n-1;
%                     v_d(n,j)=exp(-(-k*del_t))*(exp(-1/e)+(1-exp(-1/e))*x(j)-exp(-(1-x(j))/e));
%                 end
%             end
%         else
%             for n=1:p+1
%                 u_d(n,:)=interp1(x_old,u_intl(n,:),x);
%                 v_d(n,:)=interp1(x_old,v_intl(n,:),x);
%             end
%         end
%         u_d;
%         v_d
%         u_sol=solution(e,N,p,x,del_t,u_d,time_counter);
%         v_smo=smooth(e,N,p,x,del_t,v_d,time_counter)
%         u_sol;
%         v_smo;
%         w_sin=u_sol-v_smo;
% 
%         w=w_sin(round(p-1),:);
%         
%         for i=2:N+1
%             j=i-1;
%             h(j)=x(i)-x(i-1);
%         end
%         h;
%         a_dis=0;
%        
%         a_dis=h(1)*(abs(2*(((w(3)-w(2))/h(2))-((w(2)-w(1))/h(1)))/(h(2)+h(1))))^(1/m); %h(i) has been changed for adjusting index
%        
%         a_dis=a_dis+h(N)*(abs(2*(((w(N+1)-w(N))/h(N))-((w(N)-w(N-1))/h(N-1)))/(h(N)+h(N-1))))^(1/m); %h(i) has been changed for adjusting index
%        
%         for i=3:N
%             a_dis=a_dis+h(i-1)*((((abs(2*(((w(i)-w(i-1))/h(i-1))-((w(i-1)-w(i-2))/h(i-2)))/(h(i-1)+h(i-2))))^(1/m)+(abs(2*(((w(i+1)-w(i))/h(i))-((w(i)-w(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/m)))/2); %h(i) has been changed for adjusting index
%         end
%         M=zeros(N+1,1);
%         a_dis;
%         for i=2:N
%             M(i)=a_dis+(abs(2*(((w(i+1)-w(i))/h(i))-((w(i)-w(i-1))/h(i-1)))/(h(i)+h(i-1))))^(1/m);
%         end
% 
%         M(1)=M(2);
%         M(N+1)=M(N);
%         M;
% 
%         for i=2:N+1
%             j=i-1;
%             H(j)=((M(i-1)+M(i))/2)*h(j);
%         end
% 
%         H;
%         h;
%         L(1)=0;
%         for i=1:N
%             sum=0;
%             for j=1:i
%                 sum=sum+H(j);
%             end
%             L(i+1)=sum;
%         end
% 
%         C_tol=N*max(H)/L(N+1);
%         C_tol;
%         C;
%         x;
%             
%         if(C_tol<=C)
%             flag=0;
%             for n=1:p+1
%                 U_full(r_count+n,:)=u_sol(n,:);
%                 u_intl(n,:)=u_sol(p+1-(n-1),:);
%                 v_intl(n,:)=v_smo(p+1-(n-1),:);
%             end
%             u_sol;
%             u_intl;
%             x_old=x;
%             r_count=time_counter*p;
%         else
%             for i=1:N+1
%                 Y(i)=(i-1)*L(N+1)/N;
%             end
%             L;
%             x
%             Y;
%             new_x=interp1(L,x,Y);
%             x=new_x;
%         end
%         u_sol;
% 
%         x(1)=0;
%         x(N+1)=1;
% 
% %U;
%     end
%     time_counter=time_counter+1;
%     flag=1;
%     
% end
% x;
% % hold off
% % u
% % U
% % x
% % 
% % err=abs(U-u)
% % err_sup=max(err)
% % P
% %  hold on
% % P;
% counter;
% for j=1:gr_fact*p+1
%     for i=1:N+1
% %         U(gr_fact*p+1-(j-1),i)=exp(-((j-1)*del_t+x(i)/sqrt(e)));
%           U(j,i)=exp(-((j-1)*del_t))*(exp(-1/e)+(1-exp(-1/e))*x(i)-exp(-(1-x(i))/e));
%     end
% end
% % for i=1:N+1
% %     plot(P(i,:),0:counter-1);
% % end
% U_full;
% U;
% 
% 
% % surf(x,t,U_full)
% 
% surf(x,t,U)
%  err=abs(U-U_full);
% err_sup=max(err);
% max_err=max(err_sup)
% % hold off
% counter;


