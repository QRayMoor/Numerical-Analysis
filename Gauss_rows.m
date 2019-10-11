function x = Gauss_rows(A,b)
% 功能：用Gausss列主元消去法解n阶线性方程组Ax=b。
n = length(b);

for k = 1:n-1 %步骤2 消元
    [max_value,max_index]=max(abs(A(k:n,k)));  %步骤2.1 选主元
    rk = k+ max_index -1;
   
    if max_value == 0  %步骤2.2 
        warning('系数矩阵奇异');
        return;
    end
    
    if rk ~= k %步骤2.3
        t = A(k,:);         %交换矩阵的第k行与第rk行（从第k列开始）
        A(k,:) = A(rk,:); 
        A(rk,:) = t;
        t = b(k);
        b(k) = b(rk);
        b(rk) = t;
    end
     
    if i == k + 1 : n        %步骤2.4
        L(i,k) = A(i,k)/A(k,k);
        A(i,k + 1:n) = A(i,k + 1:n) - L(i,k)*A(k,k+1:n);
        b(i)= b(i) - L(i,k)*b(k);
    end     
end

if A(n,n) == 0
    warning('系数矩阵奇异！');
    return;
end

for k = n:-1:1       %回代求解
    if k == n
        x(n) = b(n)/A(n,n);
    else
        x(k) = (b(k) - sum(A(k,k+1:n) .* x(k+1:n)))/A(k,k);
    end
end