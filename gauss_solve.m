function [RA,RB,n,X]=gauss_solve(A,b)   %b为按行输入的向量，在函数中进行转置处理
    TAB = 'gauss_solve:';
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
           X=zeros(n,1);
           
           %消去过程开始
           for step = 1:n-1
               %列主元的查找
               col_max=find(abs(B(step:end, step))==max(abs(B(step:end, step))))+step-1;
               %列主元的交换
               swap_temp=B(col_max,:);
               B(col_max,:)=B(step,:);
               B(step,:)=swap_temp;
               
               for row = step+1 : n
                   m = B(row,step)/B(step,step);
                   B(row,step:n+1)=B(row,step:n+1)-m.*B(step,step:n+1); %选用n+1的目的是同时求增广矩阵中的b矩阵
               end
               
           end
           %消去过程结束
           
           %回代过程开始
           %重写A矩阵和b矩阵
           A=B(1:n,1:n);
           b=B(1:n,n+1);
           X(n,1)=b(n,1)/A(n,n);
           for k=n-1:-1:1
               X(k,1)=(b(k,1)-sum(A(k,k+1:n)*X(k+1:n,1)))/A(k,k);
           end
           %回代过程结束
       end
    end
    disp('方程组的解为：');X