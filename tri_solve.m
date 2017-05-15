function [RA,RB,n,X]=tri_solve(A,b)   %b为按行输入的向量，在函数中进行转置处理
    TAB = 'tri_solve:';
    B=[A b'];
    n=length(b);
    RA=rank(A);
	RB=rank(B);
    has_solution = (RB == RA);
    if ~has_solution
       sprintf('%s%s',TAB,'方程组无解');
    else
       if (RA~=n)
           sprintf('%s%s',TAB,'方程组有无穷多解');
       else
           sprintf('%s%s',TAB,'方程组有唯一解');
%            U=zeros(n,n);
           P=eye(n,n);
           X=zeros(n,1);
           y=zeros(n,1);
           b = b';
           
          % U(1,1:n) = A(1,1:n);
          % L(2:n,1) = A(2:n,1)/U(1,1);
           for r = 1:n
               for col_u = r:n
                    A(r,col_u)=A(r,col_u)-sum(A(r,1:r-1)*A(1:r-1,col_u));
               end
               if r<n
                   for row_l = r+1:n
                        A(row_l,r)=(A(row_l,r)-sum(A(row_l,1:r-1)*A(1:r-1,r)));%/U(r,r);
                   end
              
               
               row_max=find(abs(A(r+1:end, r))==max(abs(A(r+1:end, r))))+r-1;
              
               %列主元的交换
               swap_temp=A(row_max,:);
               A(row_max,:)=A(r,:);
               A(r,:)=swap_temp;
               
               swap_temp=P(row_max,:);
               P(row_max,:)=P(r,:);
               P(r,:)=swap_temp;

               for row_l = r+1:n
                    A(row_l,r)=A(row_l,r)/A(r,r);
               end
               end
           end
           U=triu(A);
           L=tril(A);
           L(logical(eye(n)))=1;
          % y = B(1:n,n+1);
          % A=B(1:n,1:n);
           P^-1
           pb = (P^-1)*b;
           y(1,1)=pb(1,1);
           for i=2:n
               y(i,1)=pb(i,1)-sum(L(i,1:i-1)*y(1:i-1,1));
           end
           
           X(n,1)=y(n,1)/U(n,n);
           for i=n-1:-1:1
               X(i,1)=(y(i,1)-sum(U(i,i+1:n)*X(i+1:n,1)))/U(i,i);
           end
       end
    end