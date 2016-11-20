function ospa = Ospa_Adapted(gmm_u, data, cutoff_c, order_p)
ind4X = {gmm_u.w}.'; 
nonZero = find([ind4X{:}] > 0);
X = {gmm_u.m}.'; 
X = X(nonZero);
%fix coordinates-velocities positions of X
for i = 1 : size(X,1)
    X{i} = [X{i}(1) X{i}(3) X{i}(2) X{i}(4)] ;
end
X = cell2mat(X)';
%add fake velocities for Y    
Y = [data'; zeros(2,size(data,1)) ]   ;
ospa = ospa_dist(X, Y, 300, 2);
end













