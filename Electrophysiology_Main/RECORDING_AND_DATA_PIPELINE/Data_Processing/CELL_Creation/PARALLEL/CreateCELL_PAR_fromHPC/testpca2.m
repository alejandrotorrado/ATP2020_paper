X = [1*ones(1,10);-40*ones(1,10);3*ones(1,10);-40*ones(1,10);2*ones(1,10);-45*ones(1,10);1*ones(1,10);-40*ones(1,10);1*ones(1,10);-40*ones(1,10);3*ones(1,10);-40*ones(1,10)];

mu = mean(X);

[eigenvectors, scores] = pca(X);

nComp = 2;
Xhat = scores(:,1:nComp) * eigenvectors(:,1:nComp)';
Xhat = bsxfun(@plus, Xhat, mu);

Xhat(1,:)