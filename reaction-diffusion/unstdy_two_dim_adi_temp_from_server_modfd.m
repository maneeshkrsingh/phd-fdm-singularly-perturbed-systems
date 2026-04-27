  clc
  clear all
  format short e 

  e=2^(-10);
  bt1=1;
  bt2=1;
  
  c1=1;
  c2=exp(-1/sqrt(e))-1;

for l1=1:2
    Nx1=16*2^(l1-1);
    Ny1=16*2^(l1-1);
    Nt1=4*2^(l1-1);
  sx=min(1/2,2*(e^(1/(2))*log(Nx1)));
  sy=min(1/2,2*(e^(1/(2))*log(Ny1)));

  for l2=1:2
      Nx=Nx1*2^(l2-1);
      Ny=Ny1*2^(l2-1);
      Nt=Nt1*2^(l2-1);
  
 x=zeros(Nx+1,1);
 y=zeros(Ny+1,1);
 t=zeros(Nt+1,1);
 t1=zeros((2*Nt)+1,1);
 hy=zeros(Ny,1);
 hx=zeros(Nx,1);
 u=zeros(Nx+1,Ny+1,Nt+1);
 u_temp=zeros(Nx+1,Ny+1,(2*Nt)+1);
 b1=zeros(Nx-1,1);
 b2=zeros(Nx-1,1);
 d1=zeros(Nx-1,1);
 d2=zeros(Nx-1,1);
 A1=zeros(Nx-1,Nx-1);
 A2=zeros(Ny-1,Ny-1);

  
for k=1:Nt+1
    t(k)=((k-1))/Nt;
end
  ht=1/(Nt);
  ht1=ht/2;
  for k=1:(2*Nt)+1
      t1(k)=(k-1)/(2*Nt);
  end
  
for i=1:(Nx/2)+1
    x(i)=2*sx*(i-1)/Nx;
   
end

for i=(Nx/2)+2:Nx+1
    x(i)=sx+(2*(i-1-(Nx/2))*(1-sx)/Nx);
    
end
for i=1:Nx
    hx(i)=x(i+1)-x(i);
end



for j=1:(Ny/2)+1
    y(j)=2*sy*(j-1)/Ny;
   
end

for j=(Ny/2)+2:Ny+1
    y(j)=sy+(2*(j-1-(Ny/2))*(1-sy)/Ny);
    
end

for j=1:Ny
    hy(j)=y(j+1)-y(j);
end





for i=1:Nx+1
    for j=1:Ny+1
        
        u(i,j,1)=(c1+(c2*y(j))-exp(-y(j)/sqrt(e)))*(c1+(c2*x(i))-exp(-x(i)/sqrt(e)));
        u_temp(i,j,1)=(c1+(c2*y(j))-exp(-y(j)/sqrt(e)))*(c1+(c2*x(i))-exp(-x(i)/sqrt(e)));

    end
end

for j=1:Ny+1
    for k=1:Nt+1

           u(1,j,k)=0;
           u(Nx+1,j,k)=0;
    end
end

for i=1:Nx+1
    for k=1:Nt+1

         u(i,1,k)=0;
         u(i,Ny+1,k)=0;
    end
end

%%%%%%%%%%%%%%
for j=1:Ny+1
    for k=1:(2*Nt)+1
           u_temp(1,j,k)=0;
           u_temp(Nx+1,j,k)=0;

    end
    
end

for i=1:Nx+1
    for k=1:(2*Nt)+1

           u_temp(i,1,k)=0;
           u_temp(i,Ny+1,k)=0;
    end
end
%%%%%%%%%%%%%%%%%%


it=0;
% 
for k=2:(Nt)+1
% for k=2:6
it;

    for j=1:Ny-1
      for i=2:Nx-2

        A1(i,i-1)=-(e*ht)/((hx(i)+hx(i+1))*hx(i));
        A1(i,i)=1+(e*ht)/((hx(i)+hx(i+1))*hx(i+1))+(e*ht)/((hx(i)+hx(i+1))*hx(i))+((ht*(x(i+1))^(bt1))/(2*hx(i+1)))+(ht/4);
        A1(i,i+1)=-(e*ht)/((hx(i)+hx(i+1))*hx(i+1))-((ht*(x(i+1))^(bt1))/(2*hx(i+1)));
       
     end
     i=1;
        A1(i,i)=1+(e*ht)/((hx(i)+hx(i+1))*hx(i+1))+(e*ht)/((hx(i)+hx(i+1))*hx(i))+((ht*(x(i+1))^(bt1))/(2*hx(i+1)))+(ht/4);
        A1(i,i+1)=-(e*ht)/((hx(i)+hx(i+1))*hx(i+1))-((ht*(x(i+1))^(bt1))/(2*hx(i+1)));

     i=Nx-1; 
        A1(i,i-1)=-(e*ht)/((hx(i)+hx(i+1))*hx(i));
        A1(i,i)=1+(e*ht)/((hx(i)+hx(i+1))*hx(i+1))+(e*ht)/((hx(i)+hx(i+1))*hx(i))+((ht*(x(i+1))^(bt1))/(2*hx(i+1)))+(ht/4);

         
    for i=2:Nx-2
           
         b1(i)=((e*ht)/((hy(j)+hy(j+1))*hy(j)))*u(i+1,j,k-1)+(1-(e*ht)/((hy(j)+hy(j+1))*hy(j+1))-(e*ht)/((hy(j)+hy(j+1))*hy(j))-((ht*(y(j+1))^(bt2))/(2*hy(j+1)))-(ht/4))*u(i+1,j+1,k-1)+((e*ht)/((hy(j)+hy(j+1))*hy(j+1))+((ht*(y(j+1))^(bt2))/(2*hy(j+1))))*u(i+1,j+2,k-1)-(ht/4)*(exp(-t1(k-1+it)))*(-(exp(-x(i+1)/sqrt(e))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))-(exp(-y(j+1)/sqrt(e))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e))))+((x(i+1)^bt1)*(c2+(1/sqrt(e))*exp(-x(i+1)/sqrt(e)))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))+((y(j+1)^ bt2)*(c2+(1/sqrt(e))*exp(-y(j+1)/sqrt(e)))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e)))))-(ht/4)*(exp(-t1(k-1+it+1)))*(-(exp(-x(i+1)/sqrt(e))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))-(exp(-y(j+1)/sqrt(e))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e))))+((x(i+1)^bt1)*(c2+(1/sqrt(e))*exp(-x(i+1)/sqrt(e)))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))+((y(j+1)^ bt2)*(c2+(1/sqrt(e))*exp(-y(j+1)/sqrt(e)))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e)))));
    end

         b1(1)=((e*ht)/((hy(j)+hy(j+1))*hy(j)))*u(2,j,k-1)+(1-(e*ht)/((hy(j)+hy(j+1))*hy(j+1))-(e*ht)/((hy(j)+hy(j+1))*hy(j))-((ht*(y(j+1))^(bt2))/(2*hy(j+1)))-(ht/4))*u(2,j+1,k-1)+((e*ht)/((hy(j)+hy(j+1))*hy(j+1))+((ht*(y(j+1))^(bt2))/(2*hy(j+1))))*u(2,j+2,k-1)-(ht/4)*(exp(-t1(k-1+it)))*(-(exp(-x(2)/sqrt(e))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))-(exp(-y(j+1)/sqrt(e))*(c1+(c2*x(2))-exp(-x(2)/sqrt(e))))+((x(2)^bt1)*(c2+(1/sqrt(e))*exp(-x(2)/sqrt(e)))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))+((y(j+1)^bt2)*(c2+(1/sqrt(e))*exp(-y(j+1)/sqrt(e)))*(c1+(c2*x(2))-exp(-x(2)/sqrt(e)))))-(ht/4)*(exp(-t1(k-1+it+1)))*(-(exp(-x(2)/sqrt(e))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))-(exp(-y(j+1)/sqrt(e))*(c1+(c2*x(2))-exp(-x(2)/sqrt(e))))+((x(2)^ bt1)*(c2+(1/sqrt(e))*exp(-x(2)/sqrt(e)))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))+((y(j+1)^ bt2)*(c2+(1/sqrt(e))*exp(-y(j+1)/sqrt(e)))*(c1+(c2*x(2))-exp(-x(2)/sqrt(e)))))-(-(e*ht)/((hx(1)+hx(2))*hx(1)))*u(1,j+1,k);
         b1(Nx-1)=((e*ht)/((hy(j)+hy(j+1))*hy(j)))*u(Nx,j,k-1)+(1-(e*ht)/((hy(j)+hy(j+1))*hy(j+1))-(e*ht)/((hy(j)+hy(j+1))*hy(j))-((ht*(y(j+1))^(bt2))/(2*hy(j+1)))-(ht/4))*u(Nx,j+1,k-1)+((e*ht)/((hy(j)+hy(j+1))*hy(j+1))+((ht*(y(j+1))^(bt2))/(2*hy(j+1))))*u(Nx,j+2,k-1)-(ht/4)*(exp(-t1(k-1+it)))*(-(exp(-x(Nx)/sqrt(e))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))-(exp(-y(j+1)/sqrt(e))*(c1+(c2*x(Nx))-exp(-x(Nx)/sqrt(e))))+((x(Nx)^bt1)*(c2+(1/sqrt(e))*exp(-x(Nx)/sqrt(e)))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))+((y(j+1)^bt2)*(c2+(1/sqrt(e))*exp(-y(j+1)/sqrt(e)))*(c1+(c2*x(Nx))-exp(-x(Nx)/sqrt(e)))))-(ht/4)*(exp(-t1(k-1+it+1)))*(-(exp(-x(Nx)/sqrt(e))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))-(exp(-y(j+1)/sqrt(e))*(c1+(c2*x(Nx))-exp(-x(Nx)/sqrt(e))))+((x(Nx)^ bt1)*(c2+(1/sqrt(e))*exp(-x(Nx)/sqrt(e)))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))+((y(j+1)^ bt2)*(c2+(1/sqrt(e))*exp(-y(j+1)/sqrt(e)))*(c1+(c2*x(Nx))-exp(-x(Nx)/sqrt(e)))))-(-(e*ht)/((hx(Nx-1)+hx(Nx))*hx(Nx))-((ht*(x(Nx))^(bt1))/(2*hx(Nx))))*u(Nx+1,j+1,k);

       
          d1=A1\b1;  
          
         for i=2:Nx
              u_temp(i,j+1,2*(k-1))=d1(i-1); %storing the value of d1(each time level) in column in u.
         end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    for i=1:Nx-1
      for j=2:Ny-2

        A2(j,j-1)=-(e*ht)/((hy(j)+hy(j+1))*hy(j));
        A2(j,j)=1+(e*ht)/((hy(j)+hy(j+1))*hy(j+1))+(e*ht)/((hy(j)+hy(j+1))*hy(j))+((ht*(y(j+1))^(bt2))/(2*hy(j+1)))+(ht/4);
        A2(j,j+1)=-(e*ht)/((hy(j)+hy(j+1))*hy(j+1))-((ht*(y(j+1))^(bt2))/(2*hy(j+1)));
       
     end
     j=1;
        A2(j,j)=1+(e*ht)/((hy(j)+hy(j+1))*hy(j+1))+(e*ht)/((hy(j)+hy(j+1))*hy(j))+((ht*(y(j+1))^(bt2))/(2*hy(j+1)))+(ht/4);
        A2(j,j+1)=-(e*ht)/((hy(j)+hy(j+1))*hy(j+1))-((ht*(y(j+1))^(bt2))/(2*hy(j+1)));

     j=Ny-1; 
        A2(j,j-1)=-(e*ht)/((hy(j)+hy(j+1))*hy(j));
        A2(j,j)=1+(e*ht)/((hy(j)+hy(j+1))*hy(j+1))+(e*ht)/((hy(j)+hy(j+1))*hy(j))+((ht*(y(j+1))^(bt2))/(2*hy(j+1)))+(ht/4);


          
         
    for j=2:Ny-2
           
         b2(j)=((e*ht)/((hx(i)+hx(i+1))*hx(i)))*u_temp(i,j+1,2*(k-1))+(1-(e*ht)/((hx(i)+hx(i+1))*hx(i+1))-(e*ht)/((hx(i)+hx(i+1))*hx(i))-((ht*(x(i+1))^(bt1))/(2*hx(i+1)))-(ht/4))*u_temp(i+1,j+1,2*(k-1))+((e*ht)/((hx(i)+hx(i+1))*hx(i+1))+((ht*(x(i+1))^(bt1))/(2*hx(i+1))))*u_temp(i+2,j+1,2*(k-1))-(ht/4)*(exp(-t1(k-1+it+1)))*(-(exp(-x(i+1)/sqrt(e))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))-(exp(-y(j+1)/sqrt(e))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e))))+((x(i+1)^bt1)*(c2+(1/sqrt(e))*exp(-x(i+1)/sqrt(e)))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))+((y(j+1)^bt2)*(c2+(1/sqrt(e))*exp(-y(j+1)/sqrt(e)))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e)))))-(ht/4)*(exp(-t1(k-1+it+2)))*(-(exp(-x(i+1)/sqrt(e))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))-(exp(-y(j+1)/sqrt(e))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e))))+((x(i+1)^bt1)*(c2+(1/sqrt(e))*exp(-x(i+1)/sqrt(e)))*(c1+(c2*y(j+1))-exp(-y(j+1)/sqrt(e))))+((y(j+1)^bt2)*(c2+(1/sqrt(e))*exp(-y(j+1)/sqrt(e)))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e)))));
    end

         b2(1)=((e*ht)/((hx(i)+hx(i+1))*hx(i)))*u_temp(i,2,2*(k-1))+(1-(e*ht)/((hx(i)+hx(i+1))*hx(i+1))-(e*ht)/((hx(i)+hx(i+1))*hx(i))-((ht*(x(i+1))^(bt1))/(2*hx(i+1)))-(ht/4))*u_temp(i+1,2,2*(k-1))+((e*ht)/((hx(i)+hx(i+1))*hx(i+1))+((ht*(x(i+1))^(bt1))/(2*hx(i+1))))*u_temp(i+2,2,2*(k-1))-(ht/4)*(exp(-t1(k-1+it+1)))*(-(exp(-x(i+1)/sqrt(e))*(c1+(c2*y(2))-exp(-y(2)/sqrt(e))))-(exp(-y(2)/sqrt(e))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e))))+((x(i+1)^bt1)*(c2+(1/sqrt(e))*exp(-x(i+1)/sqrt(e)))*(c1+(c2*y(2))-exp(-y(2)/sqrt(e))))+((y(2)^bt2)*(c2+(1/sqrt(e))*exp(-y(2)/sqrt(e)))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e)))))-(ht/4)*(exp(-t1(k-1+it+2)))*(-(exp(-x(i+1)/sqrt(e))*(c1+(c2*y(2))-exp(-y(2)/sqrt(e))))-(exp(-y(2)/sqrt(e))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e))))+((x(i+1)^bt1)*(c2+(1/sqrt(e))*exp(-x(i+1)/sqrt(e)))*(c1+(c2*y(2))-exp(-y(2)/sqrt(e))))+((y(2)^bt2)*(c2+(1/sqrt(e))*exp(-y(2)/sqrt(e)))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e)))))-(-(e*ht)/((hy(1)+hy(2))*hy(1)))*u(i+1,1,k);
         b2(Ny-1)=((e*ht)/((hx(i)+hx(i+1))*hx(i)))*u_temp(i,Ny,2*(k-1))+(1-(e*ht)/((hx(i)+hx(i+1))*hx(i+1))-(e*ht)/((hx(i)+hx(i+1))*hx(i))-((ht*(x(i+1))^(bt1))/(2*hx(i+1)))-(ht/4))*u_temp(i+1,Ny,2*(k-1))+((e*ht)/((hx(i)+hx(i+1))*hx(i+1))+((ht*(x(i+1))^(bt1))/(2*hx(i+1))))*u_temp(i+2,Ny,2*(k-1))-(ht/4)*(exp(-t1(k-1+it+1)))*(-(exp(-x(i+1)/sqrt(e))*(c1+(c2*y(Ny))-exp(-y(Ny)/sqrt(e))))-(exp(-y(Ny)/sqrt(e))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e))))+((x(i+1)^bt1)*(c2+(1/sqrt(e))*exp(-x(i+1)/sqrt(e)))*(c1+(c2*y(Ny))-exp(-y(Ny)/sqrt(e))))+((y(Ny)^bt2)*(c2+(1/sqrt(e))*exp(-y(Ny)/sqrt(e)))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e)))))-(ht/4)*(exp(-t1(k-1+it+2)))*(-(exp(-x(i+1)/sqrt(e))*(c1+(c2*y(Ny))-exp(-y(Ny)/sqrt(e))))-(exp(-y(Ny)/sqrt(e))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e))))+((x(i+1)^bt1)*(c2+(1/sqrt(e))*exp(-x(i+1)/sqrt(e)))*(c1+(c2*y(Ny))-exp(-y(Ny)/sqrt(e))))+((y(Ny)^bt2)*(c2+(1/sqrt(e))*exp(-y(Ny)/sqrt(e)))*(c1+(c2*x(i+1))-exp(-x(i+1)/sqrt(e)))))-(-(e*ht)/((hy(Ny-1)+hy(Ny))*hy(Ny))-((ht*(y(Ny))^(bt2))/(2*hy(Ny))))*u(i+1,Ny+1,k);
 


          d2=A2\b2;  
          
         for j=2:Ny
              u(i+1,j,k)=d2(j-1); %storing the value of d1(each time level) in column in u.
         end
    end
     
  it=it+1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end
     if(l2==1)
         F=u;
     end
     if(l2==2)
         G=u;
     end



 
  end
   max_err(l1)=max(max(max(abs(F-G(1:2:(2*Nx1)+1,1:2:(2*Ny1)+1,1:2:(2*Nt1)+1)))))

   if(l1>1)
       cnv_rt(l1-1)=log2(max_err(l1-1)/max_err(l1))
   end
end
 surf(x,y,u(:,:,Nt+1)')
