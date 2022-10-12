Npt = 160; Nh = 2;
WinLen = 240;
dcor = 20;     % Correlation distance in meter

% Add 3D correlated shadow fading
%  - 3D Filter design
stepSizeZ = 3;
dx0 = [0: (1/(Npt - 1) * WinLen): dcor * 2]; dx1 = [-dx0(end:-1:1), dx0(2:end)];
dz0 = [0: stepSizeZ: dcor * 2]; dz1 = [-dz0(end:-1:1), dz0(2:end)];
dX = repmat(dx1(:), 1, length(dx1));
dY = repmat(dx1(:).', length(dx1), 1);
Hcor = zeros(length(dx1), length(dx1), length(dz1));
for i = 1:length(dz1)
    Hcor(:, :, i) = exp(- sqrt(dX.^2 + dY.^2 + dz1(i)^2) / dcor);
end
Hcor = Hcor / norm(Hcor(:));

sigma = 0.1;
SFnoise3 = randn(Npt, Npt, Nh) * sigma;
Shadowing3 = imfilter(SFnoise3, Hcor);

figure(1),
surf(Shadowing3(:,:,1), 'edgecolor', 'none');
view(0, 90);

figure(2),
stem(svd(Shadowing3(:,:,1)));