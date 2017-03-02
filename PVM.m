%% PVM implementation as seen in http://blog.piekniewski.info/2016/11/04/predictive-vision-in-a-nutshell/
% Started writing by Jeff Rodny 2-28-2017
% Organization:
% will need to implement autoencoders, which take, as their output target
%   the next timestep's input 

%% example code to read a gif
% figure(1)
% [matrix, map] = imread('2.gif','frames',20);
% % plot that matrix into a figure, using default colormap
% image(matrix)
% % update that colormap to the colormap of the gif
% colormap(map)
%%
% clear all


% load gif into memory
[gif_matrix, color_map] = imread('1.gif','frames','all');

%initialize the autoencoder structure
% auto_encoder.params.gif_size = length(gif_matrix)-2;
% auto_encoder.params.reach_of_autoencoder = 1;%1=9area,2=25area,3=49area
auto_encoder.params.auto_encoders_per_layer = [48,24,12,6,3,1];%[94:-2:1];
auto_encoder.params.num_layers = size(auto_encoder.params.auto_encoders_per_layer,2);
% auto_encoder.params.size_1_layer_1 = 9;
% auto_encoder.params.size_3_layer_1 = 9;
auto_encoder.params.size_3 = 40;
auto_encoder.params.size_2 = 10;
auto_encoder.params.size_1 = 10 * auto_encoder.params.size_2;
auto_encoder.params.gif_time_length = size(gif_matrix,4);
auto_encoder.layer_1 = cell(auto_encoder.params.num_layers,1);
auto_encoder.layer_2 = cell(auto_encoder.params.num_layers,1);
auto_encoder.layer_2_old = cell(auto_encoder.params.num_layers,1);
auto_encoder.layer_3 = cell(auto_encoder.params.num_layers,1);
auto_encoder.weight_matrix_1_2 = cell(auto_encoder.params.num_layers,1);
auto_encoder.weight_matrix_2_3 = cell(auto_encoder.params.num_layers,1);
% auto_encoder.weight_matrix_between_layers = cell(auto_encoder.params.num_layers-1,1);
for ii = 1:auto_encoder.params.num_layers
    auto_encoder.layer_1{ii} = cell(auto_encoder.params.auto_encoders_per_layer(ii));
    auto_encoder.layer_2{ii} = cell(auto_encoder.params.auto_encoders_per_layer(ii));
    auto_encoder.layer_2_old{ii} = cell(auto_encoder.params.auto_encoders_per_layer(ii));
    auto_encoder.layer_3{ii} = cell(auto_encoder.params.auto_encoders_per_layer(ii));
    auto_encoder.weight_matrix_1_2{ii} = cell(auto_encoder.params.auto_encoders_per_layer(ii));
    auto_encoder.weight_matrix_2_3{ii} = cell(auto_encoder.params.auto_encoders_per_layer(ii));
%     if ii < auto_encoder.params.num_layers
%         auto_encoder.weight_matrix_between_layers{ii} = rand(auto_encoder.params.auto_encoders_per_layer(ii),auto_encoder.params.auto_encoders_per_layer(ii+1));
%     end
    for jj = 1:length(auto_encoder.layer_1{ii})
        for kk = 1:length(auto_encoder.layer_1{ii})
%             if ii == 1 % set bottom layer to 9, all other to 12
%                 auto_encoder.layer_1{ii}{jj,kk} = zeros(auto_encoder.params.size_1_layer_1,1);
%                 auto_encoder.layer_3{ii}{jj,kk} = zeros(auto_encoder.params.size_3_layer_1,1);
%                 auto_encoder.weight_matrix_2_3{ii}{jj,kk} = rand(auto_encoder.params.size_3_layer_1,auto_encoder.params.size_2);
%                 auto_encoder.weight_matrix_1_2{ii}{jj,kk} = rand(auto_encoder.params.size_2,auto_encoder.params.size_1_layer_1);
%             else
                auto_encoder.layer_1{ii}{jj,kk} = rand(auto_encoder.params.size_1,1);
                auto_encoder.layer_3{ii}{jj,kk} = rand(auto_encoder.params.size_3,1);
                auto_encoder.weight_matrix_2_3{ii}{jj,kk} = randn(auto_encoder.params.size_3,auto_encoder.params.size_2);
                auto_encoder.weight_matrix_1_2{ii}{jj,kk} = randn(auto_encoder.params.size_2,auto_encoder.params.size_1);
%             end
            auto_encoder.layer_2{ii}{jj,kk} = rand(auto_encoder.params.size_2,1);
            auto_encoder.layer_2_old{ii}{jj,kk} = rand(auto_encoder.params.size_2,1);
        end
    end
end


%% train
% provide one round of input
numTotalRuns = 100;
tstart = tic;
for runTimes = 1:numTotalRuns
    disp(['runTime = ',num2str(runTimes)]);
    disp([' - ', num2str(100*(runTimes-1)/numTotalRuns),'%']);
    for num_gifs = 1:6
        [gif_matrix, color_map] = imread([num2str(num_gifs),'.gif'],'frames','all');
        auto_encoder.params.gif_time_length = size(gif_matrix,4);
        disp(['gif_num = ',num2str(num_gifs)]);
        disp([' - - ', num2str(100*(num_gifs-1)/6),'%']);
        for tt = 2:auto_encoder.params.gif_time_length
            if mod(tt,auto_encoder.params.gif_time_length/10)<1
                disp([' - - - ', num2str(100*(tt-1)/auto_encoder.params.gif_time_length),'%']);
            end
            tmp_input = gif_matrix(:,:,1,tt-1);
            tmp_target = gif_matrix(:,:,1,tt);
            is_Learn_On = 1;
            iterateStep
        %     auto_encoder = iterateStep(auto_encoder, gif_matrix(:,:,1,tt),
        %     gif_matrix(:,:,1,tt+1), 1);
        end
    end
    disp([' Time Remaining (HH:MM:SS):  ', ...
                datestr((toc(tstart)/(runTimes/numTotalRuns)-toc(tstart))/24/60/60, 'HH:MM:SS')] )
end
%% test
return;

% keyboard;
% iterateStep(auto_encoder, gif_matrix(:,:,1,tt));

figure(1);
image(gif_matrix(:,:,1,1));
colormap(color_map);
tmp_input = gif_matrix(:,:,1,1);
iterateStep
for tt = 2:10
    disp(['time = ',num2str(tt)]);
    
    tmp_input = zeros(96,96);
    for ii = 1:auto_encoder.params.auto_encoders_per_layer(1)
        for jj = 1:auto_encoder.params.auto_encoders_per_layer(1)
            tmp_input(ii*2-1,jj*2-1) = auto_encoder.layer_3{1}{ii,jj}(1)*256;
            tmp_input(ii*2,jj*2-1) = auto_encoder.layer_3{1}{ii,jj}(2)*256;
            tmp_input(ii*2-1,jj*2) = auto_encoder.layer_3{1}{ii,jj}(4)*256;
            tmp_input(ii*2,jj*2) = auto_encoder.layer_3{1}{ii,jj}(5)*256;
        end
    end
    
    figure(1);
    image(tmp_input);
    colormap(color_map);
    
    is_Learn_On = 0;
    iterateStep
end

