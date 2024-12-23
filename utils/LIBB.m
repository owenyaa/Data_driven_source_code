%% Author : TAO ZHANG  * zt1996nic@gmail.com *
% Created Time : 2023-05-11 08:58
% Last Revised : TAO ZHANG ,2023-06-01
% Remark : Library of Parametric oscillator or Lorenz and so on


function [Theta,Sym]=LIBB(X,X_OrderMax,Trig_OrderMax,nonsmooth_OrderMax)
% 数据维度
[DataN,N_X]=size(X);
% 库对应的符号表示  存储符号
Sym_X=sym('x',[N_X,1]);
Sym_sin=sym('sin');
Sym_cos=sym('cos');
Sym_sign=sym('sign');

Theta=[];
%% abs(x+1) abs(x-1)
Index=0;
for i=1:N_X
    Index=Index+1;
    Theta(:,Index)=abs(X(:,i)+1);
    Sym{1,Index}=Sym_X(i,1)+1;
end
for i=1:N_X
    Index=Index+1;
    Theta(:,Index)=abs(X(:,i)-1);
    Sym{1,Index}=Sym_X(i,1)-1;
end
%% 位移代数表示
%一阶
if X_OrderMax>=1
    for i=1:N_X
        Index=Index+1;
        Theta(:,Index)=X(:,i);
        Sym{1,Index}=Sym_X(i,1);
    end
end

%二阶
if X_OrderMax>=2
    for i=1:N_X
        for j=i:N_X
            Index=Index+1;
            Theta(:,Index)=X(:,i).*X(:,j);
            Sym{1,Index}=Sym_X(i,1)*Sym_X(j,1);
        end
    end
end

%三阶
if X_OrderMax>=3
    for i=1:N_X
        for j=i:N_X
            for k=j:N_X
                Index=Index+1;
                Theta(:,Index)=X(:,i).*X(:,j).*X(:,k);
                Sym{1,Index}=Sym_X(i,1)*Sym_X(j,1)*Sym_X(k,1);
            end
        end
    end
end

%四阶
if X_OrderMax>=4
    for i=1:N_X
        for j=i:N_X
            for k=j:N_X
                for m=k:N_X
                    Index=Index+1;
                    Theta(:,Index)=X(:,i).*X(:,j).*X(:,k).*X(:,m);
                    Sym{1,Index}=Sym_X(i,1)*Sym_X(j,1)*Sym_X(k,1)*Sym_X(m,1);
                end
            end
        end
    end
end

%五阶
if X_OrderMax>=5
    for i=1:N_X
        for j=i:N_X
            for k=j:N_X
                for m=k:N_X
                    for p=m:N_X
                        Index=Index+1;
                        Theta(:,Index)=X(:,i).*X(:,j).*X(:,k).*X(:,m).*X(:,p);
                        Sym{1,Index}=Sym_X(i,1)*Sym_X(j,1)*Sym_X(k,1)*Sym_X(m,1)*Sym_X(p,1);
                    end
                end
            end
        end
    end
end


%% 三角函数表示
% sinx cosx
if Trig_OrderMax>=1
    for i=1:N_X
        Index=Index+1;
        Theta(:,Index)=sin(X(:,i));
        Sym{1,Index}=Sym_sin*Sym_X(i,1);
    end
    for i=1:N_X
        Index=Index+1;
        Theta(:,Index)=cos(X(:,i));
        Sym{1,Index}=Sym_cos*Sym_X(i,1);
    end
end

%% 非光滑表示
% sign 
if nonsmooth_OrderMax>=1
    for i=1:N_X
        Index=Index+1;
        Theta(:,Index)=sign(X(:,i));
        Sym{1,Index}=Sym_sign*Sym_X(i,1);
    end
end