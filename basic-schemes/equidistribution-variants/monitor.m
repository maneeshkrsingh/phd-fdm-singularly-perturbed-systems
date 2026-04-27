function [x,flag]=monitor(x,h,U1,U2,Nx,e1,e2)
C=1.1;
        M=zeros(Nx+1,1);
        M(1)=M(2);
        M(Nx+1)=M(Nx);
        for j=2:Nx-2
            i=j+1;
            A1=((2/(h(i)+h(i-1)))*sqrt(((U1(i+1)-U1(i))/h(i))-((U1(i)-U1(i-1))/h(i-1))));
            B1=e1^(1/4)*(2/(h(i)+h(i-1)))*sqrt(((U1(i+1))*(1/h(i))^2)+((U1(i))*(-(1/(h(i)^2))- (1/(h(i-1)^2))-(1/(h(i)*h(i-1)))))+((U1(i-1))*((1/(h(i)*h(i-1)))+(1/(h(i-1)^2))+(1/(h(i-2)*h(i-1)))))+U1(i-2)*(-1/(h(i-1)*h(i-2))));
            C1=((U1(i)-U1(i-1))/h(i-1));
            A2=((2/(h(i)+h(i-1)))*sqrt(((U2(i+1)-U2(i))/h(i))-((U2(i)-U2(i-1))/h(i-1))));
            B2=e2^(1/4)*(2/(h(i)+h(i-1)))*sqrt(((U2(i+1))*(1/h(i))^2)+((U2(i))*(-(1/(h(i)^2))- (1/(h(i-1)^2))-(1/(h(i)*h(i-1)))))+((U2(i-1))*((1/(h(i)*h(i-1)))+(1/(h(i-1)^2))+(1/(h(i-2)*h(i-1)))))+U2(i-2)*(-1/(h(i-1)*h(i-2))));
            C2=((U2(i)-U2(i-1))/h(i-1));
            M(i)=(1+A1+B1+C1+A2+B2+C2)^2;
        end
        H=zeros(Nx,1);
        for i=2:Nx+1
            j=i-1;
            H(j)=((M(i-1)+M(i))/2)*h(j);           
        end
        L=zeros(Nx+1,1);
        L(1)=0;
        for i=1:Nx
            sum=0;
            for j=1:i
                sum=sum+H(j);
            end
            L(i+1)=sum;
        end
         C_tol=Nx*max(H)/L(Nx+1);
        if(C_tol<=C)
            flag=0;
            x_old=x;
        else
            for i=1:Nx+1
                Y(i)=(i-1)*L(Nx+1)/Nx;
            end
            new_x=interp1(L,x,Y);
            x=new_x;
        end
