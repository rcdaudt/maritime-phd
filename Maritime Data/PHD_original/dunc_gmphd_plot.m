function [] = dunc_gmphd_plot(gm_prev, gm_current, assoc_thresh)
    figure(420);hold on; box on; grid on;

    % For each previous object, find its closest neighbor
    for i = 1:size(gm_prev, 2)
        mahals = [];
        for j = 1:size(gm_current, 1)
            diff = gm_current(j).m - gm_prev(i).m;
            mahals(j) = diff' \ gm_prev(i).C * diff;
        end
        
        [min, idx] = min(mahals);
        if min < assoc_thresh
            plot([gm_prev(i).m(1), gm_current(idx).m(1)], ...
                 [gm_prev(i).m(3), gm_current(idx).m(3)] );
        end
    end

end