function ospa = Ospa_Adapted(gmm_u, data, cutoff_c, order_p)
ind4X = {gmm_u.w}.'; 
% nonZero = find([ind4X{:}] > 0);
X = {gmm_u.m}.'; 
X = X([ind4X{:}] > 0);
%fix coordinates-velocities positions of X
for i = 1 : size(X,1)
    X{i} = [X{i}(1) X{i}(3)] ;
end
X = cell2mat(X)';
Y = data';
ospa = ospa_dist(X, Y, cutoff_c, order_p);
end













