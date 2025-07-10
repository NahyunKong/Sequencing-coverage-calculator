%% Determine required supporting read
k0 = 1; % k0 is how many supporting read are required to confirm existance

k0 = input('Enter how many supporting read are required (k0): ');
VAF_input = input('Enter the variant allele frequency (VAF): ');
n_input = input('Enter the total read depth (n): ');
pE = input('Enter the sequencing error rate (pE): ');

if isempty(k0); k0 = 1; end
if isempty(VAF_input); VAF_input = 0.01; end
if isempty(n_input); n_input = 500; end
if isempty(pE); pE = 0.001; end

if ~(isnumeric(k0) && isscalar(k0) && k0 > 0 && mod(k0, 1) == 0)
    warning('K0 must be a positive integer. You entered: %g', k0);
end
if ~(isnumeric(VAF_input) && isscalar(VAF_input) && VAF_input > 0 && VAF_input < 1 )
    warning('VAF must be a value between 0 and 1. You entered: %g', VAF_input);
end
if ~(isnumeric(n_input) && isscalar(n_input) && n_input > 0 && mod(n_input, 1) == 0)
    warning('Depth(n) must be a positive integer. You entered: %g', n_input);
end
if ~(isnumeric(pE) && isscalar(pE) && pE > 0 && pE < 1 )
    warning('Sequencing error must be a value between 0 and 1. You entered: %g', pE);
end

fprintf(['Running model with k0 = %d, VAF = %g, read depth = %d, ' ...
    'and sequencing error = %g\n.\n.\n.\n'],...
    k0, VAF_input, n_input, pE);

%% Figure 4; Heatmap of VAF(%) x N(#Seq)
VAF = [0.1 0.05 0.025 0.01 0.005 0.0025 0.001];
n_seq = [30 60 100 200 300 400 500 600 1000];
% pE = 0.001; % Assuming sequencing error rate is 0.1%
RR = NaN(length(n_seq), length(VAF));
warning('off', 'MATLAB:nchoosek:LargeCoefficient');

for j = 1:length(n_seq)
    n_temp = n_seq(j);
    p_detect = 0;
    for k = k0:n_temp
        p_detect = p_detect + ...
            nchoosek(n_temp,k)*(VAF*(1-pE)+(1-VAF)*pE/3).^k.*...
            (1-(VAF*(1-pE)+(1-VAF)*pE/3)).^(n_temp-k);
    end
    RR(j,:) = p_detect;
end
RR = RR*100; % convert unit to percent(%)

figure,
h = heatmap(RR');
% h.CellLabelFormat = '%.2e';
title(['Variant detection rate (%; required supporting read = ' num2str(k0), ')'] ) 
h.YLabel = 'VAF (%)';
h.YDisplayLabels = {'10', '5', '2.5', '1', '0.5','0.25', '0.1'};
h.XLabel = 'Coverage (X)';
h.XDisplayLabels = ({'30', '60', '100', '200','300' ,'400', '500', '600', '1000'});
J = customcolormap_preset('red-white-blue');
h.Colormap = J;

%%
n_temp = n_input;
p_detect = 0;
for k = k0:n_temp
    p_detect = p_detect + ...
        nchoosek(n_temp,k)*(VAF_input*(1-pE)+(1-VAF_input)*pE/3).^k.*...
        (1-(VAF_input*(1-pE)+(1-VAF_input)*pE/3)).^(n_temp-k);
end

fprintf('======== Variant detection rate ========\n');
fprintf('Variant detection rate will be: %f\n', p_detect)
fprintf('========================================\n');