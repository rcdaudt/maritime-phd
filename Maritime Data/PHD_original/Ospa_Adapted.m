function ospa = Ospa_Adapted(data, gt, cutoff_c, order_p)
ospa = zeros(1,size(data,1));
for i = 1:size(data,1)
    Nobj_data = size(data{i},1);
    Nobj_gt = size(gt{i},1);    
    X = [data{i}'; zeros(2,Nobj_data)];
    Y = [gt{i}'; zeros(2,Nobj_gt)];
    ospa(i) = ospa_dist(X, Y, cutoff_c, order_p); 

end
end













