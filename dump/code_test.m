clc; clear all;

cd 'C:\Users\Sanghun Jee\Desktop\Misha\GUI_data_examplefish\subject_6'

%Let's get the data
data = open('data_full.mat');
stim = data.data.stim_full; %Stim-> would be a criteria to sort 'data' for each trial
CellResp = h5read('TimeSeries.h5', '/CellResp'); %CellResp matrix(normalized calcium response)
absIX = h5read('TimeSeries.h5', '/absIX'); %constant matrix to match CellResp with X, Y, Z
r = 3; %rank of the cp_als

X = data.data.CellXYZ_norm(:, 1); %X coordinates
Y = data.data.CellXYZ_norm(:, 2); %Y coordinates
Z = data.data.CellXYZ_norm(:, 3); %Z coordinates

X = X(absIX());
Y = Y(absIX());
Z = Z(absIX());
fpsec = data.data.fpsec; %Frame per second, because CellResp(:, i) represents each frame
[Resp_n, Resp_d] = size(CellResp);
% figure(1)
% for i = 1:d
%     scatter3(X, Y, Z, 0.5, CellResp);
%     drawnow
% end

%Reorganize the response matrix by each stim_type
stim_1 = CellResp(:, find(stim==1)); 
stim_2 = CellResp(:, find(stim==2));
stim_3 = CellResp(:, find(stim==3));

%prestep for Reorganization the stim matrix
%Check the stim type and its iterations
current_stim = 1;
stim_num = 1;
frame_num = 1;

for i=1:Resp_d-1
    if stim(i) == stim(i+1)
        ;
    else
        stim_num = stim_num +1;
    end
end
real_stim = zeros(stim_num, 2);
for i =1:Resp_d-1
    if stim(i) == stim(i+1)
        frame_num = frame_num +1;
    else
        real_stim(current_stim, 1) = stim(i);
        real_stim(current_stim, 2) = frame_num;
        frame_num = 1;
        current_stim = current_stim +1;
    end
end
if i==Resp_d-1
    real_stim(current_stim, 1) = stim(i);
    real_stim(current_stim, 2) = frame_num;
end

%Let's make real_stim_N = (CellResp, Time(Frame), Trial)
%real_stim_N is the 'original data'
num_stim1 = size(find(real_stim(:, 1)==1), 1);
num_stim2 = size(find(real_stim(:, 1)==2), 1);
num_stim3 = size(find(real_stim(:, 1)==3), 1);

real_stim_1 = zeros(length(CellResp), real_stim(2, 2), num_stim1);
real_stim_2 = zeros(length(CellResp), real_stim(2, 2), num_stim2);
real_stim_3 = zeros(length(CellResp), real_stim(1, 2), num_stim3);

for i=1:num_stim1 %40을 인수로 바꾸어줘야함
    real_stim_1(:, :, i) = stim_1(:, 1+(i-1)*40:40*i);
    real_stim_2(:, :, i) = stim_2(:, 1+(i-1)*40:40*i);

end

for i=1:27*2
    real_stim_3(:, :, i) = stim_3(:, 1+(i-1)*30:30*i);
end

%%Prestep for TCA

%SVD, Let's get the W and B
[U1, ~, V1] = svd(stim_1-mean(stim_1, 2), 'econ'); %U:W(Neural activity), V:B(Temporal)
[U2, ~, V2] = svd(stim_2-mean(stim_2, 2), 'econ'); %U:W(Neural activity), V:B(Temporal)
[U3, ~, V3] = svd(stim_3-mean(stim_3, 2), 'econ'); %U:W(Neural activity), V:B(Temporal)

temp_stim1 = U1*V1';
temp_stim2 = U2*V2';
temp_stim3 = U3*V3';

%Make matrix reorganized by cp_als
X_hat1 = zeros(length(CellResp), real_stim(2, 2), num_stim1);
X_hat2 = zeros(length(CellResp), real_stim(2, 2), num_stim2);
X_hat3 = zeros(length(CellResp), real_stim(1, 2), num_stim3);

for i=1:num_stim1 %40을 인수로 바꾸어줘야함
    X_hat1(:, :, i) = temp_stim1(:, 1+(i-1)*40:40*i);
    X_hat2(:, :, i) = temp_stim2(:, 1+(i-1)*40:40*i);
end

for i=1:27*2
    X_hat3(:, :, i) = temp_stim3(:, 1+(i-1)*30:30*i);
end

%Convert them to tensor
T_stim1 = tensor(X_hat1);
T_stim2 = tensor(X_hat2);
T_stim3 = tensor(X_hat3);

Ori_stim1 = tensor(real_stim_1);

% temp_err = zeros(10, 3);
temp_sim = zeros(10, 3);
figure(1)
for r=2:100
    %Do cp_als
    TCA_hat1 = cp_als(T_stim1, r-1);
%     TCA_hat2 = cp_als(T_stim2, r-1);
%     TCA_hat3 = cp_als(T_stim3, r-1);
%     TCA_hat11 = cp_als(T_stim1, r);
%     TCA_hat22 = cp_als(T_stim2, r);
%     TCA_hat33 = cp_als(T_stim3, r);

    %Calculate error score
    X_Xhat1 = minus(full(Ori_stim1), full(TCA_hat1));
    
%     X_Xhat2 = minus(full(real_stim_2), full(TCA_hat2));
%     X_Xhat3 = minus(full(real_stim_3), full(TCA_hat3));

    err1 = norm(X_Xhat1(:), 'fro')^2/norm(real_stim_1(:), 'fro')^2;
%     err2 = square(norm(X_Xhat2(:), 'fro'))/square(norm(real_stim_2(:), 'fro'));
%     err3 = square(norm(X_Xhat3(:), 'fro'))/square(norm(real_stim_3(:), 'fro'));
    
    temp_err(r) = [err1;];
%     temp_TCA_hat(r) = [TCA_hat1; TCA_hat2; TCA_hat3;];
    
    %Similarity score
%     sim_1 = score(TCA_hat1, TCA_hat11);
%     sim_2 = score(TCA_hat2, TCA_hat22);
%     sim_3 = score(TCA_hat3, TCA_hat33);
    
%     temp_sim(r) = [sim_1;];
    plot(temp_err)
    drawnow
end





% rank = 4;
% %cp_arls
% arls_stim1 = cp_arls(T_stim1, rank);
% arls_stim2 = cp_arls(T_stim2, rank); 
% arls_stim3 = cp_arls(T_stim3, rank); 
% 
% % Error rate
% Xhat = 

% %Visualization
% % vizopts = {'PlotCommands',{'bar','line','line'},...
% %     'ModeTitles',{'Concentration','Emission','Excitation'},...
% %     'BottomSpace',0.10,'HorzSpace',0.04,'Normalize',0};
% % info1 = viz(M1,'Figure',1,vizopts{:});
% info1 = viz(arls_stim1, 'Figure', 4);
% info2 = viz(arls_stim2, 'Figure', 5);
% info3 = viz(arls_stim3, 'Figure', 6);