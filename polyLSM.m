function [a]=polyLSM(X,Y,n)   %X,YΪ����������������ں����н���ת�ô���
    DESC = 'polyLSM:';   %��������
    point_count = length(X);    %��ȡ���ݵ��ܸ���
    C=zeros(point_count,n+1);   %�����м���󣬳�ʼ��Ԫ�ؾ�Ϊ0
    
    X=X';   %X��Y����Ԥ����
    Y=Y';
    C(:,1) = 1;     %C�����һ�м���
    C(:,2) = X(:,1);    %C����ڶ��м���
    for i = 3:n+1   %C��������е�������
         C(:,i) = C(:,i-1).*X(:,1);
    end
    
    A=C'*C;     %���㷽��ϵ������
    b=C'*Y;
    
    a=A^(-1)*b;     %�ⷽ�̵õ��������ʽϵ��an
    a=flipud(a);        %����ʽϵ������an���·�ת��ת�ã�Ϊ��ӡ������ʽ��׼��
    a=a';
    poly_res=poly2sym(a);       %��ӡ������ʽ
    disp(DESC);
    disp('��϶���ʽΪ��');
    disp(poly_res);
    
    Y1=polyval(a,X);        %�������ʽ��Ϸ��������ƫ��
    delta = (Y1-Y);
    maxdelta = max(delta);
    delta = delta.^2;
    squareDelta = sqrt(sum(delta));     %�������ʽ��Ϸ����ľ������
    disp('������');
    disp(squareDelta);
    disp('���ƫ�');
    disp(maxdelta);
