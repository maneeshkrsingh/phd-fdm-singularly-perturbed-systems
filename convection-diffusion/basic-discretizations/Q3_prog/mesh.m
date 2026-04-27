function [x,h]=mesh(sig1,sig2,dis,N,e)
% sig_o=1;
% sig1=min(dis/2,sig_o*e*log(N));
% sig2=min((1-dis)/2,sig_o*e*log(N));
x=zeros(1,N+1);
h=zeros(1,N);

% for i=1:(T/del_t)+1
%     t(i)=(i-1)*del_t;
% end
for i=1:N/4
   h(i)=4*(dis-sig1)/N;
end

for i=N/4+1:N/2
   h(i)=4*sig1/N;
end

for i=N/2+1:3*N/4
   h(i)=4*sig2/N;
end

for i=3*N/4+1:N
   h(i)=4*(1-dis-sig2)/N;
end

x(1)=0;
for i=1:N
   x(i+1)=x(i)+h(i);
end