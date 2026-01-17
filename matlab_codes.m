% Generate Mueller matrix maps for tight focusing (Mf) and combined system (Mf * Mwpc)

% Grid setup
N = 100;
kx = linspace(-1, 1, N);
ky = linspace(-1, 1, N);
[KX, KY] = meshgrid(kx, ky);
K = sqrt(KX.^2 + KY.^2);
PHI = atan2(KY, KX);

% Annular mask for NA range 0.8 to 0.92
mask = (K >= 0.8) & (K <= 0.92);

% Parameters for tight focusing (Mf)
df = 0.1; 
xf = sqrt(1 - df^2);

% Mf computation
M11_f = ones(size(K));
M12_f = df * cos(2*PHI);
M13_f = df * sin(2*PHI);
M14_f = zeros(size(K));
M21_f = M12_f;
M22_f = cos(2*PHI).^2 + xf * sin(2*PHI).^2;
M23_f = (1 - xf) * cos(2*PHI) .* sin(2*PHI);
M24_f = zeros(size(K));
M31_f = M13_f;
M32_f = M23_f;
M33_f = sin(2*PHI).^2 + xf * cos(2*PHI).^2;
M34_f = zeros(size(K));
M41_f = zeros(size(K));
M42_f = zeros(size(K));
M43_f = zeros(size(K));
M44_f = xf * ones(size(K));

% Apply mask
M_f_elements = {M11_f, M12_f, M13_f, M14_f, ...
                M21_f, M22_f, M23_f, M24_f, ...
                M31_f, M32_f, M33_f, M34_f, ...
                M41_f, M42_f, M43_f, M44_f};
for i = 1:16
    M_f_elements{i} = M_f_elements{i} .* mask;
    M_f_elements{i}(~mask) = NaN;
end

% Plot Figure S1
figure('Name', 'Figure S1', 'NumberTitle', 'off');
for i = 1:4
    for j = 1:4
        subplot(4,4,(i-1)*4 + j);
        imagesc(kx, ky, M_f_elements{(i-1)*4 + j});
        axis square; colormap('jet'); colorbar;
        title(sprintf('M_{%d%d}', i, j));
        xlabel('k_x'); ylabel('k_y');
        set(gca, 'YDir', 'normal');
    end
end
sgtitle('Figure S1: Mf for Tight Focusing (NA 0.8-0.92)');

% WPC parameters
dwpc = 0.5; delta_wpc = pi/4;
M_wpc = [1, dwpc, 0, 0;
         dwpc, 1, 0, -sin(delta_wpc);
         0, 0, sqrt(1 - dwpc^2)*cos(delta_wpc), 0;
         0, sin(delta_wpc), 0, cos(delta_wpc)];

% Combined M = Mf * Mwpc
M_elements = cell(4,4);
for i = 1:N
    for j = 1:N
        phi = PHI(i,j);
        Mf = [1, df*cos(2*phi), df*sin(2*phi), 0;
              df*cos(2*phi), cos(2*phi)^2 + xf*sin(2*phi)^2, (1-xf)*cos(2*phi)*sin(2*phi), 0;
              df*sin(2*phi), (1-xf)*cos(2*phi)*sin(2*phi), sin(2*phi)^2 + xf*cos(2*phi)^2, 0;
              0, 0, 0, xf];
        M = Mf * M_wpc;
        for r = 1:4
            for c = 1:4
                if i == 1 && j == 1
                    M_elements{r,c} = zeros(N,N);
                end
                M_elements{r,c}(i,j) = M(r,c);
            end
        end
    end
end

% Apply mask to combined M
for i = 1:4
    for j = 1:4
        M_elements{i,j} = M_elements{i,j} .* mask;
        M_elements{i,j}(~mask) = NaN;
    end
end

% Plot Figure S2
figure('Name', 'Figure S2', 'NumberTitle', 'off');
for i = 1:4
    for j = 1:4
        subplot(4,4,(i-1)*4 + j);
        imagesc(kx, ky, M_elements{i,j});
        axis square; colormap('jet'); colorbar;
        title(sprintf('M_{%d%d}', i, j));
        xlabel('k_x'); ylabel('k_y');
        set(gca, 'YDir', 'normal');
    end
end
sgtitle('Figure S2: Combined Mf × Mwpc (λ = 530 nm)');