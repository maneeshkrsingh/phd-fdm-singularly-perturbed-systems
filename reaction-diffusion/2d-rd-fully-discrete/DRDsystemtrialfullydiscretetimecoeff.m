  clc
  clear all
  
  format short e 
  
  
  tic

   e1=2^(-30);
   e2=2^(-18);

   
   
  for N_counter =1:1
    N1=128*2^(N_counter-1);
    
    Nt1=512*4^(N_counter-1);
 
  sigx1=sqrt(e1)*log(N1);
  sigx2=sqrt(e2)*log(N1);
  sigy1=sqrt(e1)*log(N1);
  sigy2=sqrt(e2)*log(N1);

for l=1:1
     N=N1*2^(l-1);
     Nt=Nt1*2^(l-1); 
   

   % N=16;
    
   % Nt=4;
    
%     
%   sigx1=min(1/8,1*sqrt(e1)*log(N));
%   sigx2=min(1/4,1*sqrt(e2)*log(N));
%   sigy1=min(1/8,1*sqrt(e1)*log(N));
%   sigy2=min(1/4,1*sqrt(e2)*log(N));


 x=zeros(N+1,1);
 y=zeros(N+1,1);
 t=zeros(Nt+1,1);
 hy=zeros(N,1);
 hx=zeros(N,1);
%  u_vect=zeros(2*(N-1),2*(N-1),1);
 u1=zeros(N+1,N+1,Nt+1);
 u2=zeros(N+1,N+1,Nt+1);
A=sparse((N-1)*(N-1));
C1=sparse((N-1)*(N-1));
C2=sparse((N-1)*(N-1));
D=sparse((N-1)*(N-1));
%M=sparse((2*N-2)*(2*N-2));
 B1=sparse((N-1)*(N-1),1);
 B2=sparse((N-1)*(N-1),1);
 P=zeros(2*N-2,1);

   % u(:,:,1)=u_d_mat(:,:,1);

%for n=1:p
%     u_vect=zeros((N-1)*(N-1),1);

for n=1:Nt+1
    t(n)=((n-1))/Nt;
end
 del_t=1/Nt;

%%%%%%%%%%%%%%%%%%%%%% x direction %%%%%%%%
for i=1:(N/8)
    hx(i)=8*sigx1/N;
end

for i=N/8+1:(N/4)
    hx(i)=8*(sigx2-sigx1)/N;
end


for i=(N/4)+1:3*N/4
    hx(i)=2*(1-2*sigx2)/N;
end

for i=3*N/4+1:(7*N/8)
    hx(i)=8*(sigx2-sigx1)/N;
end

for i=7*N/8+1:(N)
    hx(i)=8*sigx1/N;
end

x(1)=0;
for i=1:N-1
     x(i+1)=x(i)+hx(i);
     
 end
x(N+1)=1;
%%%%%%%%%%%%%%%%%%%%%% y direction %%%%%%%%

for j=1:(N/8)
    hy(j)=8*sigy1/N;
end

for j=N/8+1:(N/4)
    hy(j)=8*(sigy2-sigy1)/N;
end


for j=(N/4)+1:3*N/4
    hy(j)=2*(1-2*sigy2)/N;
end

for j=3*N/4+1:(7*N/8)
    hy(j)=8*(sigy2-sigy1)/N;
end

for j=7*N/8+1:(N)
    hy(j)=8*sigy1/N;
end

y(1)=0;
for j=1:N-1
     y(j+1)=y(j)+hy(j);
end
 y(N+1)=1;

 
 %%%%%%%%%%%%%%%%%%%%%intial conditions and boundary conditon for both steps
%%%%%%%%%%%%%%%%%%%%%

for i=1:N+1
    for j=1:N+1
        
        %u1(i,j,1)=0;              %%%%%%%%%%%%%%%% %u1(x,y,0)=0
        %u2(i,j,1)=0;              %%%%%%%%%%%%%%%% %u2(x,y,0)=0
        u1(i,j,1)=x(i)+y(j);              
        u2(i,j,1)=(x(i)^2)*(y(j)^2);  
    
    end
end

%%%%%%%%%%%%%%% Homogeneous boundary condition%%%%%%%%%%%%%

for j=1:N+1
    for k=1:Nt+1

           u1(1,j,k)=0;          %%%%%%%%%%%%%%%% %u1(0,y,t)=0
           u2(1,j,k)=0;           %%%%%%%%%%%%%%%% %u2(0,y,t)=0
           u1(N+1,j,k)=0;        %%%%%%%%%%%%%%%% %u2(1,y,t)=0
           u2(N+1,j,k)=0;        %%%%%%%%%%%%%%%% %u2(1,y,t)=0
    end
end

% %%%%%%%%%%%%%%% Non Homogeneous boundary condition for x variable %%%%%%%%%%%%%
% for j=1:N+1
%     for k=1:Nt+1
% 
%            u1(1,j,k)=(t(1)^2)*(1-t(k))*(y(j));       %%%%%%%%%%%%%%%% %u1(0,y,t)=u1(0,y,t)=t^{2}(1-t)(x+y) take x=0
%            u2(1,j,k)=(t(1)^3)*(1-t(k)^2)*(+y(j));     %%%%%%%%%%%%%%%% %u2(0,y,t)=u2(0,y,t)=t^{3}(1-t^{2})(x+y) take x=0
%            u1(N+1,j,k)=(t(k)^2)*(1-t(k))*(1+y(j));      %%%%%%%%%%%%%%%% %u2(1,y,t)=0
%            u2(N+1,j,k)=(t(k)^3)*(1-t(k)^2)*(1+y(j));    %%%%%%%%%%%%%%%% %u2(1,y,t)=0
%     end
% end
% 
% 
% %%%%%%%%%%%%%%% Non Homogeneous boundary condition for y variable %%%%%%%%%%%%%
% for i=1:N+1
%     for k=1:Nt+1
% 
%            u1(1,j,k)=(t(1)^2)*(1-t(k))*(x(i));        %%%%%%%%%%%%%%%%% u1(x,0,t)=t^{2}(1-t)(x+y) take y=0
%            u2(1,j,k)=(t(1)^3)*(1-t(k)^2)*(x(i));      %%%%%%%%%%%%%%%%% u2(x,0,t)=t^{3}(1-t^{2})(x+y) take y=0
%            u1(N+1,j,k)=(t(k)^2)*(1-t(k))*(x(i)+1);     %%%%%%%%%%%%%%%%% u2(x,1,t)=0
%            u2(N+1,j,k)=(t(k)^3)*(1-t(k)^2)*(x(i)+1);   %%%%%%%%%%%%%%%%% u2(x,1,t)=0
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 for n=2:Nt+1
 
 

%%%%%%%%%%%%%%%%%%%% %%%%%%%   matrix A     %%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%


  for i=2:N
        fx=(x(i+1)-x(i)); %h(i+1)
        bx=(x(i)-x(i-1));  %h(i)
        cx=(x(i+1)-x(i-1));  %h(i+1)+h(i)
    for j=2:N
        k=j-1;
        fy=(y(j+1)-y(j));   %h(j+1)
        by=(y(j)-y(j-1));    %h(j)
        cy=(y(j+1)-y(j-1));   %h(j+1)+h(j)
        if(i==2)
            if(k==1)
                A(k,k)=del_t*((2*e1/((bx+fx)*bx))+(2*e1/((bx+fx)*fx))+(2*e1/((by+fy)*by))+(2*e1/((by+fy)*fy))+2+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
                A(k,k+1)=-del_t*2*e1/(fy*(by+fy));
                A(k,k+(N-1))=-del_t*2*e1/(fx*(bx+fx));          %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            if(k~=1 && k~=N-1)
                A(k,k-1)=-del_t*2*e1/(by*(by+fy));
                A(k,k)=del_t*((2*e1/((bx+fx)*bx))+(2*e1/((bx+fx)*fx))+(2*e1/((by+fy)*by))+(2*e1/((by+fy)*fy))+2+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
                A(k,k+1)=-del_t*2*e1/(fy*(by+fy));
                A(k,k+(N-1))=-del_t*2*e1/(fx*(bx+fx));          %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            if(k==N-1)
                A(k,k-1)=-del_t*2*e1/(fy*(by+fy));
                A(k,k)=del_t*((2*e1/((bx+fx)*bx))+(2*e1/((bx+fx)*fx))+(2*e1/((by+fy)*by))+(2*e1/((by+fy)*fy))+2+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
                A(k,k+(N-1))=-del_t*2*e1/(fx*(bx+fx));           %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            
        end
        
        if(i~=2 &&  i~=N)
            if(k==1)             
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k)=del_t*((2*e1/((bx+fx)*bx))+(2*e1/((bx+fx)*fx))+(2*e1/((by+fy)*by))+(2*e1/((by+fy)*fy))+2+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=-del_t*2*e1/(fy*(by+fy));
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=-del_t*2*e1/(bx*(bx+fx));                  %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k+(N-1))=-del_t*2*e1/(fx*(bx+fx));                      %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))         
            end
            if(k~=1 && k~=N-1)
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=-del_t*2*e1/(by*(by+fy));             
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k)=del_t*((2*e1/((bx+fx)*bx))+(2*e1/((bx+fx)*fx))+(2*e1/((by+fy)*by))+(2*e1/((by+fy)*fy))+2+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=-del_t*2*e1/(fy*(by+fy));
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=-del_t*2*e1/(bx*(bx+fx));                         %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k+(N-1))=-del_t*2*e1/(fx*(bx+fx));                         %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            if(k==N-1)
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=-del_t*2*e1/(by*(by+fy));           
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k)=del_t*((2*e1/((bx+fx)*bx))+(2*e1/((bx+fx)*fx))+(2*e1/((by+fy)*by))+(2*e1/((by+fy)*fy))+2+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=-del_t*2*e1/(bx*(bx+fx));                           %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k+(N-1))=-del_t*2*e1/(fx*(bx+fx));                           %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            
        end 
        
        if(i==N)
            if(k==1)             
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k)=del_t*((2*e1/((bx+fx)*bx))+(2*e1/((bx+fx)*fx))+(2*e1/((by+fy)*by))+(2*e1/((by+fy)*fy))+2+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=-del_t*2*e1/(fy*(by+fy));
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=-del_t*2*e1/(bx*(bx+fx));          %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
            end
            if(k~=1 && k~=N-1)
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=-del_t*2*e1/(by*(by+fy));                                 
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k)=del_t*((2*e1/((bx+fx)*bx))+(2*e1/((bx+fx)*fx))+(2*e1/((by+fy)*by))+(2*e1/((by+fy)*fy))+2+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=-del_t*2*e1/(fy*(by+fy));
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=-del_t*2*e1/(bx*(bx+fx));          %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
            end
            if(k==N-1)
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=-del_t*2*e1/(by*(by+fy));                                              
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k)=del_t*((2*e1/((bx+fx)*bx))+(2*e1/((bx+fx)*fx))+(2*e1/((by+fy)*by))+(2*e1/((by+fy)*fy))+2+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
             A((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=-del_t*2*e1/(bx*(bx+fx));          %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
            end
        end
    end
  end
  
  A;
  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%   matrix D     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  for i=2:N
        fx=(x(i+1)-x(i)); %h(i+1)
        bx=(x(i)-x(i-1));  %h(i)
        cx=(x(i+1)-x(i-1));  %h(i+1)+h(i)
    for j=2:N
        k=j-1;
        fy=(y(j+1)-y(j));   %h(j+1)
        by=(y(j)-y(j-1));    %h(j)
        cy=(y(j+1)-y(j-1));   %h(j+1)+h(j)
        if(i==2)
            if(k==1)
                D(k,k)=del_t*((2*e2/((bx+fx)*bx))+(2*e2/((bx+fx)*fx))+(2*e1/((by+fy)*by))+(2*e2/((by+fy)*fy))+4+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
                D(k,k+1)=-del_t*2*e2/(fy*(by+fy));
                D(k,k+(N-1))=-del_t*2*e2/(fx*(bx+fx));          %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            if(k~=1 && k~=N-1)
                D(k,k-1)=-del_t*2*e2/(by*(by+fy));
                D(k,k)=del_t*((2*e2/((bx+fx)*bx))+(2*e2/((bx+fx)*fx))+(2*e1/((by+fy)*by))+(2*e2/((by+fy)*fy))+4+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
                D(k,k+1)=-del_t*2*e2/(fy*(by+fy));
                D(k,k+(N-1))=-del_t*2*e2/(fx*(bx+fx));          %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            if(k==N-1)
                D(k,k-1)=-del_t*2*e2/(fy*(by+fy));
                D(k,k)=del_t*((2*e2/((bx+fx)*bx))+(2*e2/((bx+fx)*fx))+(2*e1/((by+fy)*by))+(2*e2/((by+fy)*fy))+4+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
                D(k,k+(N-1))=-del_t*2*e2/(fx*(bx+fx));           %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            
        end
        
        if(i~=2 &&  i~=N)
            if(k==1)             
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k)=del_t*((2*e2/((bx+fx)*bx))+(2*e2/((bx+fx)*fx))+(2*e2/((by+fy)*by))+(2*e2/((by+fy)*fy))+4+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=-del_t*2*e2/(fy*(by+fy));
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=-del_t*2*e2/(bx*(bx+fx));                  %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k+(N-1))=-del_t*2*e2/(fx*(bx+fx));                      %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))         
            end
            if(k~=1 && k~=N-1)
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=-del_t*2*e2/(by*(by+fy));             
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k)=del_t*((2*e2/((bx+fx)*bx))+(2*e2/((bx+fx)*fx))+(2*e2/((by+fy)*by))+(2*e2/((by+fy)*fy))+4+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=-del_t*2*e2/(fy*(by+fy));
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=-del_t*2*e2/(bx*(bx+fx));                         %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k+(N-1))=-del_t*2*e2/(fx*(bx+fx));                         %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            if(k==N-1)
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=-del_t*2*e2/(by*(by+fy));           
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k)=del_t*((2*e2/((bx+fx)*bx))+(2*e2/((bx+fx)*fx))+(2*e2/((by+fy)*by))+(2*e2/((by+fy)*fy))+4+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=-del_t*2*e2/(bx*(bx+fx));                           %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k+(N-1))=-del_t*2*e2/(fx*(bx+fx));                           %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            
        end 
        
        if(i==N)
            if(k==1)             
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k)=del_t*((2*e2/((bx+fx)*bx))+(2*e2/((bx+fx)*fx))+(2*e2/((by+fy)*by))+(2*e2/((by+fy)*fy))+4+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=-del_t*2*e2/(fy*(by+fy));
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=-del_t*2*e2/(bx*(bx+fx));          %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
            end
            if(k~=1 && k~=N-1)
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=-del_t*2*e2/(by*(by+fy));                                 
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k)=del_t*((2*e2/((bx+fx)*bx))+(2*e2/((bx+fx)*fx))+(2*e2/((by+fy)*by))+(2*e2/((by+fy)*fy))+4+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=-del_t*2*e2/(fy*(by+fy));
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=-del_t*2*e2/(bx*(bx+fx));          %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
            end
            if(k==N-1)
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=-del_t*2*e2/(by*(by+fy));                                               
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k)=del_t*((2*e2/((bx+fx)*bx))+(2*e2/((bx+fx)*fx))+(2*e2/((by+fy)*by))+(2*e2/((by+fy)*fy))+4+(t(n)^2)*(x(i)^2)*(y(j)^2))+1;
             D((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=-del_t*2*e2/(bx*(bx+fx));          %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
            end
        end
    end
  end  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% C1 matrix  %%%%%%%%%%%%%%
  
  for i=2:N
        fx=(x(i+1)-x(i)); %h(i+1)
        bx=(x(i)-x(i-1));  %h(i)
        cx=(x(i+1)-x(i-1));  %h(i+1)+h(i)
    for j=2:N
        k=j-1;
        fy=(y(j+1)-y(j));   %h(j+1)
        by=(y(j)-y(j-1));    %h(j)
        cy=(y(j+1)-y(j-1));   %h(j+1)+h(j)
        if(i==2)
            if(k==1)
                C1(k,k)=-del_t*(1+t(n)*(x(i)^2)*(y(j)^2));
                C1(k,k+1)=0;
                C1(k,k+(N-1))=0;         
            end
            if(k~=1 && k~=N-1)
                C1(k,k-1)=0;
                C1(k,k)=-del_t*(1+t(n)*(x(i)^2)*(y(j)^2));
                C1(k,k+1)=0;
                C1(k,k+(N-1))=0;          %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            if(k==N-1)
                C1(k,k-1)=0;
                C1(k,k)=-del_t*(1+t(n)*(x(i)^2)*(y(j)^2));
                C1(k,k+(N-1))=0;           %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            
        end
        
        if(i~=2 &&  i~=N)
            if(k==1)             
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k)=-del_t*(1+t(n)*(x(i)^2)*(y(j)^2));
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=0;
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=0;                  %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k+(N-1))=0;                      %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))         
            end
            if(k~=1 && k~=N-1)
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=0;             
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k)=-del_t*(1+t(n)*(x(i)^2)*(y(j)^2));
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=0;
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=0;                         %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k+(N-1))=0;                         %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            if(k==N-1)
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=0;           
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k)=-del_t*(1+t(n)*(x(i)^2)*(y(j)^2));
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=0;                           %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k+(N-1))=0;                           %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            
        end 
        
        if(i==N)
            if(k==1)             
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k)=-del_t*(1+t(n)*(x(i)^2)*(y(j)^2));
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=0;
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=0;          %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
            end
            if(k~=1 && k~=N-1)
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=0;                                 
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k)=-del_t*(1+t(n)*(x(i)^2)*(y(j)^2));
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=0;
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=0;          %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
            end
            if(k==N-1)
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=0;                                                
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k)=-del_t*(1+t(n)*(x(i)^2)*(y(j)^2));
             C1((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=0;          %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
            end
        end
    end
  end  
  
  
  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% C2 matrix  %%%%%%%%%%%%%%
  
  for i=2:N
        fx=(x(i+1)-x(i)); %h(i+1)
        bx=(x(i)-x(i-1));  %h(i)
        cx=(x(i+1)-x(i-1));  %h(i+1)+h(i)
    for j=2:N
        k=j-1;
        fy=(y(j+1)-y(j));   %h(j+1)
        by=(y(j)-y(j-1));    %h(j)
        cy=(y(j+1)-y(j-1));   %h(j+1)+h(j)
        if(i==2)
            if(k==1)
                C2(k,k)=-del_t*t(n)*((x(i)^2)*(y(j)^2));
                C2(k,k+1)=0;
                C2(k,k+(N-1))=0;         
            end
            if(k~=1 && k~=N-1)
                C2(k,k-1)=0;
                C2(k,k)=-del_t*t(n)*((x(i)^2)*(y(j)^2));
                C2(k,k+1)=0;
                C2(k,k+(N-1))=0;          %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            if(k==N-1)
                C2(k,k-1)=0;
                C2(k,k)=-del_t*t(n)*((x(i)^2)*(y(j)^2));
                C2(k,k+(N-1))=0;           %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            
        end
        
        if(i~=2 &&  i~=N)
            if(k==1)             
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k)=-del_t*t(n)*((x(i)^2)*(y(j)^2));
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=0;
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=0;                  %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k+(N-1))=0;                      %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))         
            end
            if(k~=1 && k~=N-1)
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=0;             
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k)=-del_t*t(n)*((x(i)^2)*(y(j)^2));
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=0;
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=0;                         %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k+(N-1))=0;                         %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            if(k==N-1)
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=0;           
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k)=-del_t*t(n)*((x(i)^2)*(y(j)^2));
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=0;                           %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k+(N-1))=0;                           %%%%%%%%%%%%%%%%%%%%%% + N-1 means u(x(i+1),y(j))
            end
            
        end 
        
        if(i==N)
            if(k==1)             
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k)=-del_t*t(n)*((x(i)^2)*(y(j)^2));
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=0;
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=0;          %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
            end
            if(k~=1 && k~=N-1)
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=0;                                 
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k)=-del_t*t(n)*((x(i)^2)*(y(j)^2));
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k+1)=0;
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=0;          %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
            end
            if(k==N-1)
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k-1)=0;                                                
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k)=-del_t*t(n)*((x(i)^2)*(y(j)^2));
             C2((i-2)*(N-1)+k,(i-2)*(N-1)+k-(N-1))=0;          %%%%%%%%%%%%%%%%%%%%%% - (N-1) means u(x(i-1),y(j))
            end
        end
    end
  end  
  
  
  
 M=[A,C2;C1,D]; 
  
  
%   m=n-1;
%    % l=0;
  for i=2:N
        fx=(x(i+1)-x(i));
        bx=(x(i)-x(i-1));
        cx=(x(i+1)-x(i-1));
    for j=2:N
        k=j-1;
        fy=(y(j+1)-y(j));
        by=(y(j)-y(j-1));
        cy=(y(j+1)-y(j-1));
    %    l=l+1;
        if(i==2)
            if(k==1)
                B1(k)=del_t*exp(t(n))*(x(i)*y(j)*(1-x(i)*(1-y(j))))+u1(i,j,n);         %BCs are zero
            end
            if(k~=1 && k~=N-1)
                B1(k)=del_t*exp(t(n))*(x(i)*y(j)*(1-x(i)*(1-y(j))))+u1(i,j,n);        %BCs are zero
            end
            if(k==N-1)
                B1(k)=del_t*exp(t(n))*(x(i)*y(j)*(1-x(i)*(1-y(j))))+u1(i,j,n);        %BCs are zero
            end  
        end
        
        if(i~=2 &&  i~=N)
            if(k==1)
                B1((i-2)*(N-1)+k)=del_t*exp(t(n))*(x(i)*y(j)*(1-x(i)*(1-y(j))))+u1(i,j,n);           %BCs are zero
            end
            if(k~=1 && k~=N-1)
                 B1((i-2)*(N-1)+k)=del_t*exp(t(n))*(x(i)*y(j)*(1-x(i)*(1-y(j))))+u1(i,j,n);  
            end
            if(k==N-1)
                B1((i-2)*(N-1)+k)=del_t*exp(t(n))*(x(i)*y(j)*(1-x(i)*(1-y(j))))+u1(i,j,n);          %BCs are zero
            end  
        end
        
        
        if(i==N)
            if(k==1)
                B1((i-2)*(N-1)+k)=del_t*exp(t(n))*(x(i)*y(j)*(1-x(i)*(1-y(j))))+u1(i,j,n);       %BCs are zero
            end
            if(k~=1 && k~=N-1)
                B1((i-2)*(N-1)+k)=del_t*exp(t(n))*(x(i)*y(j)*(1-x(i)*(1-y(j))))+u1(i,j,n);         %BCs are zero
            end
            if(k==N-1)
                 B1((i-2)*(N-1)+k)=del_t*exp(t(n))*(x(i)*y(j)*(1-x(i)*(1-y(j))))+u1(i,j,n);        %BCs are zero
            end 
        end 
    end
  end
  
  
  
   for i=2:N
        fx=(x(i+1)-x(i));
        bx=(x(i)-x(i-1));
        cx=(x(i+1)-x(i-1));
    for j=2:N
        k=j-1;
        fy=(y(j+1)-y(j));
        by=(y(j)-y(j-1));
        cy=(y(j+1)-y(j-1));
    %    l=l+1;
        if(i==2)
            if(k==1)
                B2(k)=del_t*exp(t(n))*((x(i)^2)*(y(j)^2)*(1-x(i)*(1-y(j))))+u2(i,j,n);         %BCs are zero
            end
            if(k~=1 && k~=N-1)
                B2(k)=del_t*exp(t(n))*((x(i)^2)*(y(j)^2)*(1-x(i)*(1-y(j))))+u2(i,j,n);          %BCs are zero
            end
            if(k==N-1)
                B2(k)=del_t*exp(t(n))*((x(i)^2)*(y(j)^2)*(1-x(i)*(1-y(j))))+u2(i,j,n);         %BCs are zero
            end  
        end
        
        if(i~=2 &&  i~=N)
            if(k==1)
                B2((i-2)*(N-1)+k)=del_t*exp(t(n))*((x(i)^2)*(y(j)^2)*(1-x(i)*(1-y(j))))+u2(i,j,n);            %BCs are zero
            end
            if(k~=1 && k~=N-1)
                 B2((i-2)*(N-1)+k)=del_t*exp(t(n))*((x(i)^2)*(y(j)^2)*(1-x(i)*(1-y(j))))+u2(i,j,n);   
            end
            if(k==N-1)
                B2((i-2)*(N-1)+k)=del_t*exp(t(n))*((x(i)^2)*(y(j)^2)*(1-x(i)*(1-y(j))))+u2(i,j,n);            %BCs are zero
            end  
        end
        
        
        if(i==N)
            if(k==1)
                B2((i-2)*(N-1)+k)=del_t*exp(t(n))*((x(i)^2)*(y(j)^2)*(1-x(i)*(1-y(j))))+u2(i,j,n);         %BCs are zero
            end
            if(k~=1 && k~=N-1)
                B2((i-2)*(N-1)+k)=del_t*exp(t(n))*((x(i)^2)*(y(j)^2)*(1-x(i)*(1-y(j))))+u2(i,j,n);         %BCs are zero
            end
            if(k==N-1)
                 B2((i-2)*(N-1)+k)=del_t*exp(t(n))*((x(i)^2)*(y(j)^2)*(1-x(i)*(1-y(j))))+u2(i,j,n);          %BCs are zero
            end 
        end 
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
   B=[B1;B2];
  
    P=M\B;
   
   for i=2:N
        for j=2:N
           u1(i,j,n)=P((i-2)*(N-1)+(j-1));
        end
   end 
    
     for i=N+1:2*N-2
        for j=N+1:2*N-1
           u2(i+1-N,j+1-N,n)=P((i-2)*(N-1)+(j-1));
        end
    end 
   
   
  end
%   
   if(l==1)
          F1=u1;
          F2=u2;
         
    end
%        if(l==2)
%           G1=u1;
%           G2=u2;
%           
%        end
       
end
  %%%%%%%%%%%%%%%%%%%%%
 
    %max_err(N_counter)=max(max(max(abs(F-G(1:2:N+1,1:2:N+1,1:2:Nt+1)))));
%     err1=max(max(abs(F1-G1(1:2:N+1,1:2:N+1,1:2:Nt+1))));
%      err_sup1=max(err1);
%          max_err1=max(err_sup1);
%       %error1(ep_counter,N_counter)=max_err1;
%        error1(N_counter)=max_err1;
%     
%     
%     err2=max(max(abs(F2-G2(1:2:N+1,1:2:N+1,1:2:Nt+1))));
%     err_sup2=max(err2);
%          max_err2=max(err_sup2);
%        %error2(ep_counter,N_counter)=max_err2;
%        error2(N_counter)=max_err2;
   
  end
  
%  for i=1:N_counter-1
%          %cnv_rt1(j,i)=log2(error1(j,i)/error1(j,i+1));
%          cnv_rt1(i)=log2(error1(i)/error1(i+1));
% end
% %end
% %for j=1:ep_counter
%  for i=1:N_counter-1
%          %cnv_rt2(j,i)=log2(error2(j,i)/error2(j,i+1));
%         cnv_rt2(i)=log2(error2(i)/error2(i+1));
%  end 
  
%    error1
%   cnv_rt1
% % 
%    error2
%   cnv_rt2
% % B;
% 
% surf(x,y,u1(:,:,Nt+1)')

toc