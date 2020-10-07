clc, clear
mu=0.1; Beta=0.1; M=10; N=10;
A = round(rand(N));
A = triu(A) + triu(A,1)';
A = A - diag(diag(A));
for i=1:N
    for j=1:N       
            L(i)=1;
            w(i)=sum(A(i,:));
            r(i,j)=1-(1-(A(i,j)/w(i)))^L(i);
    end
end

p(1,:)=rand(N,1);
%%
q_ant=1;
for t=1:M 
    ro(t)=sum(p(t,:))/N;
    for i=1:N
        for j=1:N
            q(t,i)=(1-Beta*r(j,i)*p(t,j))*q_ant;
            q_ant=q(t,i); 
        end
        p(t+1,i)=(1-q(t,i))*(1-p(t,i))+(1-mu)*(p(t,i))+mu*(1-q(t,i))*(p(t,i));
    end
end



