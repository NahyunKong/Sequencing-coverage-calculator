%% VAF model fitting depending on k0

% Model: Binomial

% Binomial: f_BI(0,n,VAF) = (1-VAF)^n;
% where, VAF: % of variant; n: number of sequencing

% k0 is how many supporting read are required to confirm existance
k0 = 1;

% True Data Site(row) by VAF%(column) array
% Site: WashU, NYGC, BCM, Broad
% VAF%: 10 5 2 1 0.5 0.25 (%)
TRUE_WASHU = [0.9923887323,0.9725584441,0.9882978969,0.9597331419,0.9624433341,0.8511306606];
TRUE_NYGC = [0.989525087	0.9653008343	0.9844185766	0.9401901719	0.8153270904	0.61025842];
TRUE_BCM = [0.9924640914	0.9712167203	0.9877231828	0.9593421159	0.9453691142	0.816663937];
TRUE_Broad = [0.987716469	0.9638964663	0.9752488076	0.8955793971	0.6991826875	0.5281700027];

data = [TRUE_WASHU;TRUE_BCM;TRUE_NYGC;TRUE_Broad];
data_name = {'TRUE_WASHU','TRUE_BCM','TRUE_NYGC','TRUE_Broad'};
data_i = 1;
TRUE_DATA = data;

% Recall Rate: RR
VAF = [0.1	0.05	0.02	0.01	0.005	0.0025];
n_WASHU = 510; n_Broad= 175; n_NYGC = 271; n_BCM = 463;
n_total = [n_WASHU, n_BCM, n_NYGC, n_Broad];
pE = 0.001; % Illumina Nova seq error rate = 0.1%

RR_WASHU_BI = NaN(length(n_total),length(VAF));
warning('off', 'MATLAB:nchoosek:LargeCoefficient');

for j = 1:length(n_total)
    p_nondetect = 0;
    for k = 0:(k0-1)
        p_nondetect = p_nondetect + ...
            nchoosek(n_total(j),k) .* VAF.^k .* (1-VAF).^(n_total(j)-k);
    end
    RR_WASHU_BI(j,:) = 1-p_nondetect;
end

% Recall Rate with sequencing error
pE = 0.001; % Illumina Nova seq error rate = 0.1%
RR_WASHU_BI_ER = NaN(length(n_total),length(VAF));
for j = 1:length(n_total)
    n_temp = n_total(j);
    p_nondetect = 0;
    for k = 0:k0-1 %k0:n_temp
        p_nondetect = p_nondetect + ...
            nchoosek(n_temp,k)*(VAF*(1-pE)+(1-VAF)*pE/3).^k.*...
            (1-(VAF*(1-pE)+(1-VAF)*pE/3)).^(n_temp-k);
    end
    RR_WASHU_BI_ER(j,:) = 1-p_nondetect;
end

% Recall Rate when sequencing error is Poisson dist.
pE = 0.001; % Illumina Nova seq error rate = 0.1%
RR_WASHU_BI_ER_PO = NaN(length(n_total),length(VAF));
for j = 1:length(n_total)
    n_temp = n_total(j);
    p_nondetect = 0;
    for k = 0:(k0-1)
        poi = 0;
        for l = 0:(k0-k-1)
            lambda = (1-VAF)*pE/3;
            poi = poi + exp(-lambda).*lambda.^l./factorial(l);
        end
        q = VAF.*(1-pE);
        p_nondetect = p_nondetect + nchoosek(n_total(j),k) .* q.^k .* (1-q).^(n_temp-k) .* poi;
    end
    RR_WASHU_BI_ER_PO(j,:) = 1-p_nondetect;
end

% Figure 1; Binomial model fitting
figure, hold on; 
b = bar(TRUE_DATA','FaceColor','flat');
xticks([1:length(VAF)]);
xticklabels({'0.1','0.05','0.02','0.01','0.005','0.0025'});
xticklabels({'10','5','2','1','0.5','0.25'});
x = 1:length(VAF);
for i = 1:length(n_total)
    plot(x+(i-2.5)*0.18,RR_WASHU_BI_ER(i,:),'rs','MarkerSize',7);
end
legend({'Binomial w/SQ-ER'})
for i = 1:length(n_total)
    plot(x+(i-2.5)*0.18,RR_WASHU_BI(i,:),'ks','MarkerSize',7);
end
legend({'WASHU','BCM','NYGC','Broad',...
    'Binomial w/SQ-ER','','','','','',...
    'Binomial','','',''},'Location','southwest')
ylim([0 1]);
title('Variant Detect Rate for different VAF(%)')
set(gca, 'XDir','reverse')
hold off;

% Figure 2; comparing between models (VAF=10%,0.25%)
figure, yl = [0 1];

subplot(1,2,1), hold on;
bar(TRUE_DATA(:,1));
plot(RR_WASHU_BI(:,1),'o-');
plot(RR_WASHU_BI_ER(:,1),'o-');
plot(RR_WASHU_BI_ER_PO(:,1),'o-');
xticks([1:length(n_total)]);
xticklabels({'WASHU','BCM','NYGC','Broad'})
legend({'True set','Binomial',...
    'Bin w/SQ-ER','Bin w/Poi-SQ-ER'},'Location','southeast');
title('Recall rate comparison between models (VAF=10%)')
ylim(yl);
hold off;

subplot(1,2,2), hold on;
bar(TRUE_DATA(:,end));
plot(RR_WASHU_BI(:,end),'o-');
plot(RR_WASHU_BI_ER(:,end),'o-');
plot(RR_WASHU_BI_ER_PO(:,end),'o-');
xticks([1:length(n_total)]);
xticklabels({'WASHU','BCM','NYGC','Broad'})
legend({'True set','Binomial',...
    'Bin w/SQ-ER','Bin w/Poi-SQ-ER'},'Location','southeast');
title('Recall rate comparison between models (VAF=0.25%)')
ylim(yl);
hold off;


%% Example data setup (replace these with your actual data):
% Suppose you have 5 sets of data, each with 6 points
N = 4;  % number of data sets
M = 6;  % number of data points in each set

%% Preallocate arrays for SSE and RMSE
SSE_model1 = zeros(N,1);
SSE_model2 = zeros(N,1);
RMSE_model1 = zeros(N,1);
RMSE_model2 = zeros(N,1);
RMSE_model2 = zeros(N,1);
RMSE_model4 = zeros(N,1);
RMSE_model3 = zeros(N,1);
RMSE_model6 = zeros(N,1);

%% Compute SSE and RMSE for each data set
for i = 1:N
    % Extract real data, model1 predictions, model2 predictions
    y_true = TRUE_DATA(i, :);
    y_pred1 = RR_WASHU_BI(i, :);
    y_pred2 = RR_WASHU_BI_ER(i, :);
    y_pred3 = RR_WASHU_BI_ER_PO(i, :);
    
    % Sum of Squared Errors for model1 and model2
    y_mean = mean(y_true);
    SStot = sum((y_true-y_mean).^2);

    SSE_model1(i) = sum((y_true - y_pred1).^2);
    SSE_model2(i) = sum((y_true - y_pred2).^2);
    SSE_model3(i) = sum((y_true - y_pred3).^2);

    R2_model3(i) = 1-(SSE_model2(i)/SStot);
    
    % Root Mean Square Error for model1 and model2
    RMSE_model1(i) = sqrt(mean((y_true - y_pred1).^2));
    RMSE_model2(i) = sqrt(mean((y_true - y_pred2).^2));
    RMSE_model3(i) = sqrt(mean((y_true - y_pred3).^2));
end

%% Aggregate metrics across all data sets
% You can either look at sums/means across all data sets or look at them individually.

total_SSE_model1 = sum(SSE_model1);

mean_SSE_model1  = mean(SSE_model1);

total_RMSE_model1 = sum(RMSE_model1);

mean_RMSE_model1  = mean(RMSE_model1);
mean_RMSE_model2  = mean(RMSE_model2);
mean_RMSE_model3  = mean(RMSE_model3);

mean_RMSE = [mean_RMSE_model1,mean_RMSE_model2,mean_RMSE_model3];


%% Figure 3; RMSE
figure, hold on;

bar(mean_RMSE,'FaceColor',[0.7 0.7 0.7]);
xticks([1:3]);
xticklabels({'Binomial','Bin w/SQ-ER',...
    'Bin w/Poi-SQ-ER'});

% yticks([]);
ylabel('RMSE')
title('RMSE for six models (lower the better)')
R2 = mean(R2_model3);

hold off;


%% Figure 4; Heatmap of VAF(%) x N(#Seq)
VAF = [0.1 0.05 0.025 0.01 0.005 0.0025 0.001];
n_seq = [30 60 100 200 300 400 500 600 1000];
pE = 0.001; % Assuming sequencing error rate is 0.1%
RR = NaN(length(n_seq), length(VAF));

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
title(['Variant detection rate (%; required supporting read = ' num2str(k0), ')'] ) 
h.YLabel = 'VAF (%)';
h.YDisplayLabels = {'10', '5', '2.5', '1', '0.5','0.25', '0.1'};
h.XLabel = 'Coverage (X)';
h.XDisplayLabels = ({'30', '60', '100', '200','300' ,'400', '500', '600', '1000'});
J = customcolormap_preset('red-white-blue');
h.Colormap = J;


function pdf = poisson(k,lambda)
    pdf = lambda^k*exp(-lambda)/factorial(k);
end

