function [RA,RB,n,X]=tri_solve(A,b)   %b为按行输入的向量，在函数中进行转置处理
    TAB = 'tri_solve:';   %函数名称，列主元三角分解法
    B=[A b'];
    n=length(b);
    RA=rank(A);
	RB=rank(B);
    has_solution = (RB == RA);  %求矩阵的秩，判断是否相等
    if ~has_solution                   %矩阵的秩不相等，则方程组无解
       sprintf('%s%s',TAB,'方程组无解');
    else
       if (RA~=n)
           sprintf('%s%s',TAB,'方程组有无穷多解');
       else
           %当方程组有唯一解的时候开始运算
           sprintf('%s%s',TAB,'方程组有唯一解');
           %初始化解向量及中间变量
           P=eye(n,n);
           X=zeros(n,1);
           y=zeros(n,1);
           
           %三角分解过程开始，第r步运算开始
           for r = 1:n
               %列主元交换子程序，只剩下最后一行的时候不用交换
               if r<n
                   %计算L矩阵第r列的各行比较量Si
                   for row_l = r:n
                        B(row_l,r)=(B(row_l,r)-sum(B(row_l,1:r-1)*B(1:r-1,r)));
                   end
                   %获取各行比较量Si最大值
                   row_max=find(abs(B(r:end, r))==max(abs(B(r:end, r))))+r-1;

                   %进行列主元的交换，同时交换P矩阵
                   swap_temp=B(row_max,:);
                   B(row_max,:)=B(r,:);
                   B(r,:)=swap_temp;

                   swap_temp=P(row_max,:);
                   P(row_max,:)=P(r,:);
                   P(r,:)=swap_temp;
               end
               
               %计算U矩阵的第r行元素，将增广矩阵中的b一起代入运算
               for col_u = r:n+1
                    B(r,col_u)=B(r,col_u)-sum(B(r,1:r-1)*B(1:r-1,col_u));
               end
               
               %计算L矩阵的第r列元素，增加对最后一列的判断以防止数组下标溢出
               if r<n
                   for row_l = r+1:n
                        B(row_l,r)=B(row_l,r)/B(r,r);
                   end
               end
           end
           %三角分解过程完成
           
           %从运算结果中获得U、L、b矩阵
           U=triu(B(1:n,1:n));      %获取上三角矩阵
           L=tril(B(1:n,1:n));        %获取下三角矩阵
           L(logical(eye(n)))=1;   %L矩阵对角线全为1
           b=B(1:n,n+1);

           %回代过程开始
           X(n,1)=b(n,1)/U(n,n);
           for i=n-1:-1:1
               X(i,1)=(b(i,1)-sum(U(i,i+1:n)*X(i+1:n,1)))/U(i,i);
           end
           %回代过程结束
       end
    end