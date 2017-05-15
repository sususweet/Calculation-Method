function [RA,RB,n,X]=gauss_solve(A,b)   %bΪ����������������ں����н���ת�ô���
    TAB = 'gauss_solve:';
    B=[A b'];
    n=length(b);
    RA=rank(A);
	RB=rank(B);
    has_solution = (RB == RA);
    if ~has_solution
       sprintf('%s%s',TAB,'�������޽�');
    else
       if (RA~=n)
           sprintf('%s%s',TAB,'��������������');
       else
           sprintf('%s%s',TAB,'��������Ψһ��');
           X=zeros(n,1);
           
           %��ȥ���̿�ʼ
           for step = 1:n-1
               %����Ԫ�Ĳ���
               col_max=find(abs(B(step:end, step))==max(abs(B(step:end, step))))+step-1;
               %����Ԫ�Ľ���
               swap_temp=B(col_max,:);
               B(col_max,:)=B(step,:);
               B(step,:)=swap_temp;
               
               for row = step+1 : n
                   m = B(row,step)/B(step,step);
                   B(row,step:n+1)=B(row,step:n+1)-m.*B(step,step:n+1); %ѡ��n+1��Ŀ����ͬʱ����������е�b����
               end
               
           end
           %��ȥ���̽���
           
           %�ش����̿�ʼ
           %��дA�����b����
           A=B(1:n,1:n);
           b=B(1:n,n+1);
           X(n,1)=b(n,1)/A(n,n);
           for k=n-1:-1:1
               X(k,1)=(b(k,1)-sum(A(k,k+1:n)*X(k+1:n,1)))/A(k,k);
           end
           %�ش����̽���
       end
    end
    disp('������Ľ�Ϊ��');X