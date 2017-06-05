function [a]=polyLSM(X,Y,n)   %X,Y为按行输入的向量，在函数中进行转置处理
    DESC = 'polyLSM:';   %函数名称
    point_count = length(X);    %获取数据点总个数
    C=zeros(point_count,n+1);   %建立中间矩阵，初始化元素均为0
    
    X=X';   %X、Y矩阵预处理
    Y=Y';
    C(:,1) = 1;     %C矩阵第一列计算
    C(:,2) = X(:,1);    %C矩阵第二列计算
    for i = 3:n+1   %C矩阵后续列迭代计算
         C(:,i) = C(:,i-1).*X(:,1);
    end
    
    A=C'*C;     %计算方程系数矩阵
    b=C'*Y;
    
    a=A^(-1)*b;     %解方程得到待求多项式系数an
    a=flipud(a);        %多项式系数矩阵an上下翻转加转置，为打印出多项式做准备
    a=a';
    poly_res=poly2sym(a);       %打印出多项式
    disp(DESC);
    disp('拟合多项式为：');
    disp(poly_res);
    
    Y1=polyval(a,X);        %计算多项式拟合方法的最大偏差
    delta = (Y1-Y);
    maxdelta = max(delta);
    delta = delta.^2;
    squareDelta = sqrt(sum(delta));     %计算多项式拟合方法的均方误差
    disp('均方误差：');
    disp(squareDelta);
    disp('最大偏差：');
    disp(maxdelta);
