%integration using Romberg method:�������㷨
function g=Romberg(f,N,a,b,e) %fΪ�����ֵĺ�����NΪ�����ִ�����aΪ�������ޣ�bΪ�������ޣ�eΪ�����
    if N<4
        error('���ִ���̫�٣�����߶��ִ��������ִ�����Ҫ���ڵ���4');
    else
        T=zeros(N+1,4);
        T(1,1) = (b-a)*(f(a)+f(b))/2;

        for k = 1:N
            m=0;
            for i=1:2^(k-1)
                m=m+f(a+(2*i-1)*(b-a)/(2^k));
            end

            T(k+1,1)=0.5*T(k,1)+(b-a)*m/2^k;

            for t=1:3
                if t<k+1
                    T(k+1,t+1)=(4^t*T(k+1,t)-T(k,t))/(4^t-1);
                end
            end

            if abs(T(k+1,4)-T(k,4))<e && k>=4
                break;
            end
        end
        
        if k == N && abs(T(k+1,4)-T(k,4))>=e
            error('���ִ��������򲻿ɻ�')
        else
            disp('����Ļ���ֵΪ��');
            disp(T(k+1,4));
            g = T(k+1,4);
        end
    end
    