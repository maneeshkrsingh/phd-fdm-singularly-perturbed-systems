clc
clear all
format short e 
%  e=2^(-30);
tao=1;
% Nx=input('Nx=');
% Nt=input('Nt=');
for d1=1:14
    
%     e(d1)=10^(-d1);
    e(d1)=2^(-2*(d1-1));
    e(d1);
    m1=exp(-1/e(d1));
    m2=1-exp(-1/e(d1));


for p=1:6
    Nx1=16*2^(p-1);
%     Nt1=16*2^(p-1);
    Nt1=10*2^(p-1);



% %%%%%%%
s1=min(1/2,tao*e(d1)*log(Nx1));
% s1

% 
% F=zeros(Nx1+1,Nt1+1);
for l=1:1

    Nx=Nx1*2^(l-1);
%     Nt=Nt1*2^(l-1);
    Nt=Nt1*2^(l-1);

A=zeros(Nx-1);
b=zeros(Nx-1,1);
d=zeros(Nx-1,1);
u=zeros(Nx+1,Nt+1);
u_exct=zeros(Nx+1,Nt+1);
f=zeros(Nx+1,Nt+1);
a=zeros(Nx+1,1);
x=zeros(Nx+1,1);
h=zeros(Nx,1);
t=zeros(Nt+1,1);







for i=1:(Nx/2)
      h(i)=(2*(1-s1))/(Nx);
end

for i=(Nx/2)+1:Nx
    h(i)=(2*s1)/(Nx);
 end


x(1)=0;
 for i=1:Nx
    x(i+1)=x(i)+h(i);
end

for j=1:Nt+1
    t(j)=(1*(j-1))/Nt;
end
t;
k=1/Nt;

for i=1:Nx+1
    u(i,1)=m1+(m2*x(i))-exp((x(i)-1)/e(d1));
end
for j=1:Nt+1
    u(1,j)=0;
    u(Nx+1,j)=0;
end
%%%%%%%%%%%%%%%%%
% for i=1:Nx+1
%     a(i)=1+(x(i)*(1-x(i)));
% end

for j=1:Nt+1
    for i=1:Nx+1
        f(i,j)=exp(-t(j))*(-m1+exp((x(i)-1)/e(d1))+m2-((x(i))*m2));
    end
end

for j=1:Nt+1
    for i=1:Nx+1
        u_exct(i,j)=exp(-t(j))*(m1+(m2*x(i))-exp((x(i)-1)/e(d1)));
    end
end
%%%%%%%%%%%%%%%%%
for i=2:Nx-2


%        if(i<=(Nx/2)) 
%            disp('yes')
        A(i,i-1)=-((2*e(d1)*k)/((h(i)+h(i+1))*h(i)))-((k)/(h(i)));
        A(i,i)=1+((2*e(d1)*k)/((h(i)+h(i+1))*h(i+1)))+((2*e(d1)*k)/((h(i)+h(i+1))*h(i)))+(k/(h(i)));
        A(i,i+1)=-((2*e(d1)*k)/((h(i)+h(i+1))*h(i+1)));   
%        end 
   
%        else
% %        if((Nx/2)+1<=i<=Nx-2)
%              disp('No')
%         A(i,i-1)=-((2*e(d1)*k)/((h(i)+h(i+1))*h(i)))-((k*a(i+1))/(h(i)+h(i+1)));
% %         A(i,i)=1+((2*e(d1)*k)/((h(i)+h(i+1))*h(i+1)))+((2*e(d1)*k)/((h(i)+h(i+1))*h(i)));
%          A(i,i)=1+((2*e(d1)*k)/((h(i)*h(i+1))));
%         A(i,i+1)=-((2*e(d1)*k)/((h(i)+h(i+1))*h(i+1)))+((k*a(i+1))/(h(i)+h(i+1)));
%        end 
%       

end
i=1;
% disp('yes1')
% %         A(i,i)=(1/2)+((2*e(d1)*k)/((h(i)+h(i+1))*h(i+1)))+((2*e(d1)*k)/((h(i)+h(i+1))*h(i)))+((k*(a(i)+a(i+1)))/(2*h(i)));
%           A(i,i)=(1/2)+((2*e(d1)*k)/((h(i)*h(i+1))))+((k*(a(i)+a(i+1)))/(2*h(i)));
%         A(i,i+1)=-((2*e(d1)*k)/((h(i)+h(i+1))*h(i+1)));
          A(i,i)=1+((2*e(d1)*k)/((h(i)+h(i+1))*h(i+1)))+((2*e(d1)*k)/((h(i)+h(i+1))*h(i)))+(k/(h(i)));
        A(i,i+1)=-((2*e(d1)*k)/((h(i)+h(i+1))*h(i+1)));   

i=Nx-1; 
%     disp('No1')
%         A(i,i-1)=-((2*e(d1)*k)/((h(i)+h(i+1))*h(i)))-((k*a(i+1))/(h(i)+h(i+1)));
% %         A(i,i)=1+((2*e(d1)*k)/((h(i)+h(i+1))*h(i+1)))+((2*e(d1)*k)/((h(i)+h(i+1))*h(i)));
%         A(i,i)=1+((2*e(d1)*k)/((h(i)*h(i+1))));
        A(i,i-1)=-((2*e(d1)*k)/((h(i)+h(i+1))*h(i)))-((k)/(h(i)));
        A(i,i)=1+((2*e(d1)*k)/((h(i)+h(i+1))*h(i+1)))+((2*e(d1)*k)/((h(i)+h(i+1))*h(i)))+(k/(h(i)));

%       for j=2:Nt+1 %check it
for j=2:Nt+1   %check itNt+1
     for i=2:Nx-2
             b(i)=(k*f(i+1,j))+u(i+1,j-1);
%          if(i<=(Nx/2)) 
%                disp('Yes2')
%              b(i)=((k/2)*(f(i+1,j)+f(i,j)))+((1/2)*(u(i+1,j-1)+u(i,j-1)));
% %          end 
%    
%          else
% %         if((Nx/2)+1<=i<=Nx-2) 
%     disp('No2')
%               b(i)=(k*f(i+1,j))+u(i+1,j-1);  
% 
%         end 
            
     end
         

%          b(1)=((k/2)*(f(2,j)+f(1,j)))+((1/2)*(u(2,j-1)+u(1,j-1)));
%          b(Nx-1)=(k*f(Nx,j))+u(Nx,j-1);
         b(1)=(k*f(2,j))+u(2,j-1);
         b(Nx-1)=(k*f(Nx,j))+u(Nx,j-1);
         
         
         
          d=A\b;


         for i=2:Nx
              u(i,j)=d(i-1); %storing the value of d(each time level) in column of u.
         end

end
% if(l==1)
%        F=u;
% end
% if(l==2)
%        G=u;
% end

   


       

end

   max_err(d1,p)=max(max(abs(u_exct-u)));
   if(p>1)
   cnv_rt(d1,p-1)=log2( max_err(d1,p-1)/ max_err(d1,p));
   end

end
max_err
% 
   cnv_rt
%     figure(d1)
%      surf(x,t,u')

end
%       surf(x,t,u');
%     plot(x,u(:,Nx+1))

     