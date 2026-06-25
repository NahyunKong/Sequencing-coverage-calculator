%% VAF model fitting depending on k0

% Model 1: Binomial
% Model 2: Beta-Binomial

% Binomial: f_BI(0,n,VAF) = (1-VAF)^n;
% Beta-Binomial: f_BB(0|n,VAF) = B(a,n-b)/B(a,b) 
% = B(VAF*n+1,n(2-VAF)+1)/B(n*VAF+1,n(1-VAF)+1)
% where, VAF: % of variant; n: number of sequencing

% k0 is how many supporting read are required to confirm existance
k0 = 5;

% True Data Site(row) by VAF%(column) array
% Site: WashU, NYGC, BCM, Broad
% VAF%: 10 5 2 1 0.5 0.25 (%)
TRUE_OMIM = [0.996225	0.994488	0.993828	0.954044	0.967148	0.865503
0.993708	0.991578	0.991689	0.93481	0.82006	0.614988
0.996225	0.994108	0.993319	0.953724	0.952713	0.825766
0.993321	0.989525	0.982921	0.892933	0.714535	0.536606];
TRUE_ACMG = [1	0.992492	0.996849	0.98843	0.99375	0.856882
0.994186	0.991491	0.997899	0.978099	0.7625	0.611668
1	0.992993	0.995798	0.990289	0.96875	0.828624
1	0.988488	0.993697	0.932645	0.7125	0.564266];
TRUE_DGI = [0.996058	0.992612	0.99369	0.987559	0.969027	0.851907
0.994035	0.989241	0.990891	0.972744	0.816926	0.610122
0.996265	0.991956	0.992655	0.987624	0.948071	0.816276
0.99336	0.987616	0.982703	0.927074	0.709499	0.531356];
TRUE_GIAB = [1	0.993637	0.978526	0.98207	0.944954	0.873096
0.992157	0.987509	0.97267	0.961847	0.816514	0.655457
0.996078	0.992458	0.97755	0.979881	0.954128	0.836929
0.988235	0.984209	0.961933	0.916919	0.733945	0.589467];
TRUE_ONCO = [0.997933	0.994379	0.993394	0.988616	0.964852	0.8615
0.9938	0.991389	0.990947	0.974979	0.809097	0.622778
0.997244	0.994628	0.992844	0.989476	0.944866	0.827839
0.996555	0.990337	0.981283	0.928829	0.696072	0.535147];
TRUE_COSMIC = [1	0.993873	0.993029	0.988401	0.97046	0.86256
0.995671	0.990982	0.991434	0.975049	0.805387	0.620625
1	0.993701	0.992189	0.989025	0.94874	0.825412
0.999038	0.990259	0.97976	0.927568	0.717637	0.539121];
TRUE_CHIP = [1	0.991558	0.984166	0.990068	1	0.868852
0.994413	0.988961	0.98173	0.97833	0.672414	0.622951
1	0.992857	0.985384	0.990971	0.948276	0.812221
1	0.988961	0.974421	0.927991	0.724138	0.515648];

data = {TRUE_OMIM,TRUE_ACMG,TRUE_DGI,TRUE_GIAB,TRUE_ONCO,...
    TRUE_COSMIC,TRUE_CHIP};
data_name = {'TRUE_OMIM','TRUE_ACMG','TRUE_DGI','TRUE_GIAB','TRUE_ONCO',...
    'TRUE_COSMIC','TRUE_CHIP'};
data_i = 1;
TRUE_DATA = data{data_i};

% Recall Rate: RR
VAF = [0.1	0.05	0.02	0.01	0.005	0.0025];
n_WASHU = 491; n_Broad= 173; n_NYGC = 260; n_BCM = 447;
n_total = [n_WASHU, n_BCM, n_NYGC, n_Broad];
pE = 0.001; % Illumina Nova seq error rate = 0.1%
RR_weight = NaN(length(n_total),length(VAF));
rmse = NaN(10,1);
for t = 1:10
    for j = 1:length(n_total)
        for i = 1:length(VAF)
            n_t = ceil(n_total(j)*(t+5)/10);
            p_nondetect = 0;
            for k = 0:k0-1 %k0:n_temp
            p_nondetect = p_nondetect + nchoosek(n_t,k) * ...
                betabinomial(k,n_t,(VAF(i)*(1-pE)+(1-VAF(i))*pE/3));
            end
            RR_weight(j,i) = 1-p_nondetect;
        end
    end
    er = TRUE_DATA - RR_weight;
    rmse(t) = sum(sqrt(er.^2),[1,2]);
end
[m,I1] = min(rmse);

RR_weight = NaN(length(n_total),length(VAF));
rmse = NaN(10,1);
for t = 1:10
    for j = 1:length(n_total)
        for i = 1:length(VAF)
            n_t = ceil(n_total(j)*((I1+5)+(t-5)/10)/10);
            p_nondetect = 0;
            for k = 0:k0-1 %k0:n_temp
            p_nondetect = p_nondetect + nchoosek(n_t,k) * ...
                betabinomial(k,n_t,(VAF(i)*(1-pE)+(1-VAF(i))*pE/3));
            end
            RR_weight(j,i) = 1-p_nondetect;
        end
    end
    er = TRUE_DATA - RR_weight;
    rmse(t) = sum(sqrt(er.^2),[1,2]);
end
[m,I2] = min(rmse);
w = (I1+5)/10+(I2-5)/100;
fprintf('Weight(%s)= %.4f\n', data_name{data_i},w);


RR_WASHU_BI = NaN(length(n_total),length(VAF));
RR_WASHU_BB = NaN(length(n_total),length(VAF));
warning('off', 'MATLAB:nchoosek:LargeCoefficient');

for j = 1:length(n_total)
    p_nondetect = 0;
    for k = 0:(k0-1)
        p_nondetect = p_nondetect + ...
            nchoosek(n_total(j),k) .* VAF.^k .* (1-VAF).^(n_total(j)-k);
    end
    RR_WASHU_BI(j,:) = 1-p_nondetect;
    for i = 1:length(VAF)
        p_nondetect = 0;
        for k = 0:(k0-1)
            p_nondetect = p_nondetect + ...
                nchoosek(n_total(j),k) .* betabinomial(k,n_total(j),VAF(i));
        end

        RR_WASHU_BB(j,i) = 1-p_nondetect;
    end
end

% Recall Rate with sequencing error
pE = 0.001; % Illumina Nova seq error rate = 0.1%
RR_WASHU_BI_ER = NaN(length(n_total),length(VAF));
RR_WASHU_BB_ER = NaN(length(n_total),length(VAF));
for j = 1:length(n_total)
    n_temp = n_total(j);
    p_nondetect = 0;
    for k = 0:k0-1 %k0:n_temp
        p_nondetect = p_nondetect + ...
            nchoosek(n_temp,k)*(VAF*(1-pE)+(1-VAF)*pE/3).^k.*...
            (1-(VAF*(1-pE)+(1-VAF)*pE/3)).^(n_temp-k);
    end

    RR_WASHU_BI_ER(j,:) = 1-p_nondetect;

    for i = 1:length(VAF)
        p_nondetect = 0;
        for k = 0:k0-1 %k0:n_temp
        p_nondetect = p_nondetect + nchoosek(n_temp,k) * ...
            betabinomial(k,n_total(j),(VAF(i)*(1-pE)+(1-VAF(i))*pE/3));
        end
        RR_WASHU_BB_ER(j,i) = 1-p_nondetect;
        
    end
end

% Recall Rate when sequencing error is Poisson dist.
pE = 0.001; % Illumina Nova seq error rate = 0.1%
RR_WASHU_BI_ER_PO = NaN(length(n_total),length(VAF));
RR_WASHU_BB_ER_PO = NaN(length(n_total),length(VAF));
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


    for i = 1:length(VAF)
        p_nondetect = 0;
        for k = 0:k0-1 
            poi = 0;
            for l = 0:(k0-k-1)
                lambda = (1-VAF(i))*pE/3;
                poi = poi + exp(-lambda).*lambda.^l./factorial(l);
            end
            p_nondetect = p_nondetect + nchoosek(n_temp,k) * ...
            betabinomial(k,n_total(j),(VAF(i)*(1-pE))) .* poi;
        end

        RR_WASHU_BB_ER_PO(j,i) = 1-p_nondetect;
    end
end

% Figure 1; Binomial model fitting
figure, hold on; 
b = bar(TRUE_DATA','FaceColor','flat');
xticks([1:length(VAF)]);
xticklabels({'0.1','0.05','0.02','0.01','0.005','0.0025'});
xticklabels({'10','5','2','1','0.5','0.25'});
% for i = 1:4
%     b(i).CData = [1 0 0];
% end
% b(5).CData = [0.8 0.8 0.8];
% b(6).CData = [0.8 0.8 0.8];
x = 1:length(VAF);
for i = 1:length(n_total)
    plot(x+(i-2.5)*0.18,RR_WASHU_BI_ER(i,:),'rs','MarkerSize',7);
end
legend({'Binomial w/SQ-ER'})
for i = 1:length(n_total)
    plot(x+(i-2.5)*0.18,RR_WASHU_BI(i,:),'ks','MarkerSize',7);
end
% legend({'WASHU','BCM','NYGC','Broad','SIM485','SIM4000',...
%     'Binomial w/SQ-ER','','','','','',...
%     'Binomial','','',''},'Location','southwest')
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
plot(RR_WASHU_BB(:,1),'o-');
plot(RR_WASHU_BI_ER(:,1),'o-');
plot(RR_WASHU_BB_ER(:,1),'o-');
plot(RR_WASHU_BI_ER_PO(:,1),'o-');
plot(RR_WASHU_BB_ER_PO(:,1),'o-');
xticks([1:length(n_total)]);
% xticklabels({'WASHU','BCM','NYGC','Broad','SIM485','SIM4000'})
xticklabels({'WASHU','BCM','NYGC','Broad'})
legend({'True set','Binomial','Beta-Binomial',...
    'Bin w/SQ-ER','Beta-Bin w/SQ-ER','Bin w/Poi-SQ-ER',...
    'Beta-Bin w/Poi-SQ-ER'},'Location','southeast');
title('Recall rate comparison between models (VAF=10%)')
ylim(yl);
hold off;

subplot(1,2,2), hold on;
bar(TRUE_DATA(:,end));
plot(RR_WASHU_BI(:,end),'o-');
plot(RR_WASHU_BB(:,end),'o-');
plot(RR_WASHU_BI_ER(:,end),'o-');
plot(RR_WASHU_BB_ER(:,end),'o-');
plot(RR_WASHU_BI_ER_PO(:,end),'o-');
plot(RR_WASHU_BB_ER_PO(:,end),'o-');
xticks([1:length(n_total)]);
% xticklabels({'WASHU','BCM','NYGC','Broad','SIM485','SIM4000'})
xticklabels({'WASHU','BCM','NYGC','Broad'})
legend({'True set','Binomial','Beta-Binomial',...
    'Bin w/SQ-ER','Beta-Bin w/SQ-ER','Bin w/Poi-SQ-ER',...
    'Beta-Bin w/Poi-SQ-ER'},'Location','southeast');
title('Recall rate comparison between models (VAF=0.25%)')
ylim(yl);
hold off;


% figure, hold on;
% bar(TRUE_DATA');
% plot(RR_WASHU_BI','o-');
% ylim([0.3 1])
% title('Binomial estimate');
% 
% figure, hold on;
% bar(TRUE_DATA');
% plot(RR_WASHU_BB','o-');
% ylim([0.3 1])
% title('Beta-Binomial estimate');



%% Example data setup (replace these with your actual data):
% Suppose you have 5 sets of data, each with 6 points
N = 4;  % number of data sets
M = 6;  % number of data points in each set

%% Preallocate arrays for SSE and RMSE
SSE_model1 = zeros(N,1);
SSE_model2 = zeros(N,1);
RMSE_model1 = zeros(N,1);
RMSE_model2 = zeros(N,1);
RMSE_model3 = zeros(N,1);
RMSE_model4 = zeros(N,1);
RMSE_model5 = zeros(N,1);
RMSE_model6 = zeros(N,1);

%% Compute SSE and RMSE for each data set
for i = 1:N
    % Extract real data, model1 predictions, model2 predictions
    y_true = TRUE_DATA(i, :);
    y_pred1 = RR_WASHU_BI(i, :);
    y_pred2 = RR_WASHU_BB(i, :);
    y_pred3 = RR_WASHU_BI_ER(i, :);
    y_pred4 = RR_WASHU_BB_ER(i, :);
    y_pred5 = RR_WASHU_BI_ER_PO(i, :);
    y_pred6 = RR_WASHU_BB_ER_PO(i, :);
    
    % Sum of Squared Errors for model1 and model2
    y_mean = mean(y_true);
    SStot = sum((y_true-y_mean).^2);

    SSE_model1(i) = sum((y_true - y_pred1).^2);
    SSE_model2(i) = sum((y_true - y_pred2).^2);
    SSE_model3(i) = sum((y_true - y_pred3).^2);
    SSE_model4(i) = sum((y_true - y_pred4).^2);
    SSE_model5(i) = sum((y_true - y_pred5).^2);
    SSE_model6(i) = sum((y_true - y_pred6).^2);

    R2_model3(i) = 1-(SSE_model3(i)/SStot);
    
    % Root Mean Square Error for model1 and model2
    RMSE_model1(i) = sqrt(mean((y_true - y_pred1).^2));
    RMSE_model2(i) = sqrt(mean((y_true - y_pred2).^2));
    RMSE_model3(i) = sqrt(mean((y_true - y_pred3).^2));
    RMSE_model4(i) = sqrt(mean((y_true - y_pred4).^2));
    RMSE_model5(i) = sqrt(mean((y_true - y_pred5).^2));
    RMSE_model6(i) = sqrt(mean((y_true - y_pred6).^2));
end

%% Aggregate metrics across all data sets
% You can either look at sums/means across all data sets or look at them individually.

total_SSE_model1 = sum(SSE_model1);
total_SSE_model2 = sum(SSE_model2);

mean_SSE_model1  = mean(SSE_model1);
mean_SSE_model2  = mean(SSE_model2);

total_RMSE_model1 = sum(RMSE_model1);
total_RMSE_model2 = sum(RMSE_model2);

mean_RMSE_model1  = mean(RMSE_model1);
mean_RMSE_model2  = mean(RMSE_model2);
mean_RMSE_model3  = mean(RMSE_model3);
mean_RMSE_model4  = mean(RMSE_model4);
mean_RMSE_model5  = mean(RMSE_model5);
mean_RMSE_model6  = mean(RMSE_model6);

mean_RMSE = [mean_RMSE_model1,mean_RMSE_model2,mean_RMSE_model3,...
    mean_RMSE_model4,mean_RMSE_model5,mean_RMSE_model6];


%% Figure 3; RMSE
figure, hold on;
bar(mean_RMSE,'FaceColor',[0.7 0.7 0.7]);
xticks([1:6]);
xticklabels({'Binomial','Beta-Binomial','Bin w/SQ-ER','Beta-Bin w/SQ-ER',...
    'Bin w/Poi-SQ-ER','Beta-Bin w/Poi-SQ-ER'});
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
% h.CellLabelFormat = '%.2e';
title(['Variant detection rate (%; required supporting read = ' num2str(k0), ')'] ) 
h.YLabel = 'VAF (%)';
h.YDisplayLabels = {'10', '5', '2.5', '1', '0.5','0.25', '0.1'};
h.XLabel = 'Coverage (X)';
h.XDisplayLabels = ({'30', '60', '100', '200','300' ,'400', '500', '600', '1000'});
J = customcolormap_preset('red-white-blue');
h.Colormap = J;

%% Display results
fprintf('=== Model Comparison Results ===\n');

fprintf('Model 1: Total SSE  = %.4f, Mean SSE  = %.4f\n', total_SSE_model1, mean_SSE_model1);
fprintf('Model 1: Total RMSE = %.4f, Mean RMSE = %.4f\n', total_RMSE_model1, mean_RMSE_model1);

fprintf('Model 2: Total SSE  = %.4f, Mean SSE  = %.4f\n', total_SSE_model2, mean_SSE_model2);
fprintf('Model 2: Total RMSE = %.4f, Mean RMSE = %.4f\n', total_RMSE_model2, mean_RMSE_model2);

if mean_SSE_model1 < mean_SSE_model2
    disp('=> Model 1 has a lower SSE and may be the better fit.');
else
    disp('=> Model 2 has a lower SSE and may be the better fit.');
end

if mean_RMSE_model1 < mean_RMSE_model2
    disp('=> Model 1 has a lower RMSE and may be the better fit.');
else
    disp('=> Model 2 has a lower RMSE and may be the better fit.');
end





function pdf = betabinomial(k,n,p)
    pdf = beta(k+n*p+1,-k+n*(2-p)+1)/beta(n*p+1,n*(1-p)+1);
end

function pdf = poisson(k,lambda)
    pdf = lambda^k*exp(-lambda)/factorial(k);
end

