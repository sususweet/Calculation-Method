function [RA,RB,n,X]=tri_solve(A,b)   %bΪ����������������ں����н���ת�ô���
    TAB = 'tri_solve:';   %�������ƣ�����Ԫ���Ƿֽⷨ
    B=[A b'];
    n=length(b);
    RA=rank(A);
	RB=rank(B);
    has_solution = (RB == RA);  %�������ȣ��ж��Ƿ����
    if ~has_solution                   %������Ȳ���ȣ��򷽳����޽�
       sprintf('%s%s',TAB,'�������޽�');
    else
       if (RA~=n)
           sprintf('%s%s',TAB,'��������������');
       else
           %����������Ψһ���ʱ��ʼ����
           sprintf('%s%s',TAB,'��������Ψһ��');
           %��ʼ�����������м����
           P=eye(n,n);
           X=zeros(n,1);
           y=zeros(n,1);
           
           %���Ƿֽ���̿�ʼ����r�����㿪ʼ
           for r = 1:n
               %����Ԫ�����ӳ���ֻʣ�����һ�е�ʱ���ý���
               if r<n
                   %����L�����r�еĸ��бȽ���Si
                   for row_l = r:n
                        B(row_l,r)=(B(row_l,r)-sum(B(row_l,1:r-1)*B(1:r-1,r)));
                   end
                   %��ȡ���бȽ���Si���ֵ
                   row_max=find(abs(B(r:end, r))==max(abs(B(r:end, r))))+r-1;

                   %��������Ԫ�Ľ�����ͬʱ����P����
                   swap_temp=B(row_max,:);
                   B(row_max,:)=B(r,:);
                   B(r,:)=swap_temp;

                   swap_temp=P(row_max,:);
                   P(row_max,:)=P(r,:);
                   P(r,:)=swap_temp;
               end
               
               %����U����ĵ�r��Ԫ�أ�����������е�bһ���������
               for col_u = r:n+1
                    B(r,col_u)=B(r,col_u)-sum(B(r,1:r-1)*B(1:r-1,col_u));
               end
               
               %����L����ĵ�r��Ԫ�أ����Ӷ����һ�е��ж��Է�ֹ�����±����
               if r<n
                   for row_l = r+1:n
                        B(row_l,r)=B(row_l,r)/B(r,r);
                   end
               end
           end
           %���Ƿֽ�������
           
           %���������л��U��L��b����
           U=triu(B(1:n,1:n));      %��ȡ�����Ǿ���
           L=tril(B(1:n,1:n));        %��ȡ�����Ǿ���
           L(logical(eye(n)))=1;   %L����Խ���ȫΪ1
           b=B(1:n,n+1);

           %�ش����̿�ʼ
           X(n,1)=b(n,1)/U(n,n);
           for i=n-1:-1:1
               X(i,1)=(b(i,1)-sum(U(i,i+1:n)*X(i+1:n,1)))/U(i,i);
           end
           %�ش����̽���
       end
    end