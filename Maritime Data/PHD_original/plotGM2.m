function plotGM2(data,gt,GM,cst,timestamp)

% pathName = '/home/isabel/Documents/MATLAB/PHD_NBclutter/simulation/';

figure(22);
clf
ind = find([GM.i]);

w = zeros(length(ind),1);
mu = zeros(length(ind),2);
sigma = zeros(2,2,length(ind)); %cat(3,[2 0;0 .5],[1 0;0 1]);
for ii=1:length(ind)
    w(ii) = GM(ind(ii)).w;
    if cst.Nx == 4
        mu(ii,:) = GM(ind(ii)).m([1,3])';
        sigma(:,:,ii) = GM(ind(ii)).C([1,3],[1,3]);
    elseif cst.Nx == 2
        mu(ii,:) = GM(ind(ii)).m';
        sigma(:,:,ii) = GM(ind(ii)).C;        
    else
        error('invalid dimension: %d',cst.Nx);
    end
    sigma(2,1,ii) = sigma(1,2,ii); %to make sure they're really symmetrical
    if ~isreal(sigma(:,:,ii))
        keyboard;
    end
end
obj1 = gmdistribution(mu,sigma,w);

scaling = 10;
[yy,xx]=meshgrid(1:1/scaling:cst.xwidth,1:1/scaling:cst.ywidth);
a1 = reshape(pdf(obj1, [xx(:),yy(:)]),size(xx,1),size(xx,2));

hold on
%plot the GMM intensity
imagesc(sqrt(a1)',[0 0.1]) 
% plot the GT
if ~isempty(gt)
    plot(scaling*gt(:,1),scaling*gt(:,2),'go');
end
%trajectories
% for ii=1:length(trajectories)
%     currenttrajectory = trajectories{ii};
%     currenttrajectory(currenttrajectory(:,5)>timestamp,:)=[];
%     if max(currenttrajectory(:,5))<timestamp
%         continue;
%     end
%     if ~isempty(currenttrajectory)
%         plot(scaling*currenttrajectory(:,1),scaling*currenttrajectory(:,3),'g','LineWidth',2);
%     end
% end
% measurements
plot(scaling*data(:,1),scaling*data(:,2),'rx');
axis equal
axis tight
axis ij
axis off

% pause(0.1)
drawnow
%     B = getframe(fig1,winsize);
%     if tt==1
%         imwrite(B.cdata,[pathName, 'movie.tif']);
%     else
%         imwrite(B.cdata,[pathName, 'movie.tif'],'WriteMode','append');
%     end
%     fprintf('t %i -- w1 %f w2 %f -- max1 %f max2 %f \n', tt, max(w1),max(w2), sqrt(max(max(a1))), sqrt(max(max(a2))))
%     ezsurf(@(x,y)pdf(obj,[x y]),windowx,windowy)
