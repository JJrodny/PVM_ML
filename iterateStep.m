% function auto_encoder = iterateStep(auto_encoder, tmp_input, tmp_target, is_Learn_On)

% % provide one round of input
% for tt = 1:auto_encoder.params.gif_time_length
%     disp(['time = ',num2str(tt)]);
    % for each layer
%     tmp_input;
    for ii = 1:auto_encoder.params.num_layers
%         disp(['layer = ',num2str(ii)]);
%         tmp_input2 = zeros(96,96);
%         for kk = 1:auto_encoder.params.auto_encoders_per_layer(1)
%             for jj = 1:auto_encoder.params.auto_encoders_per_layer(1)
%                 tmp_input2(ii*2-1,jj*2-1) = auto_encoder.layer_3{1}{ii,jj}(1)*256;
%                 tmp_input2(ii*2,jj*2-1) = auto_encoder.layer_3{1}{ii,jj}(2)*256;
%                 tmp_input2(ii*2-1,jj*2) = auto_encoder.layer_3{1}{ii,jj}(4)*256;
%                 tmp_input2(ii*2,jj*2) = auto_encoder.layer_3{1}{ii,jj}(5)*256;
%             end
%         end
%         disp([num2str(sum(sum(tmp_input))), ' ', num2str(sum(sum(tmp_input2)))]);
%         if abs(sum(sum(tmp_input))-sum(sum(tmp_input2)))<50
%             is_Learn_On = -1;
%         end
        % for each Xautoencoder on that layer
        for Xcntr = 1:auto_encoder.params.auto_encoders_per_layer(ii)
            % for each Yautoencoder on that layer
            for Ycntr = 1:auto_encoder.params.auto_encoders_per_layer(ii)
            
                %% if it's an input layer, give it the input
                if ii == 1
                    % set the center point
                    tmpCenter_idx = [Xcntr,Ycntr];
                    %% 1,1
                    % 1,1 1,2 1,3
                    % 2,1 2,2 2,3
                    % 3,1 3,2 3,3
                    
                    %% 1,2
                    % 1,3 1,4 1,5
                    % 2,3 2,4 2,5
                    % 3,3 3,4 3,5
                    
                    %% 
                    
                    if Xcntr == 1
                        tmpCenter_idx(1) = Xcntr+1;
                    elseif Xcntr == auto_encoder.params.auto_encoders_per_layer(ii)
                        tmpCenter_idx(1) = Xcntr-1;
                    end
                    if Ycntr == 1
                        tmpCenter_idx(2) = Ycntr+1;
                    elseif Ycntr == auto_encoder.params.auto_encoders_per_layer(ii)
                        tmpCenter_idx(2) = Ycntr-1;
                    end
                    
                    % set the indices for that point and the points around it
                    tmp_indices = [...
                        tmpCenter_idx(1)*2-2,tmpCenter_idx(2)*2-2;
                        tmpCenter_idx(1)*2-1,tmpCenter_idx(2)*2-2;
                        tmpCenter_idx(1)*2,tmpCenter_idx(2)*2-2;
                        tmpCenter_idx(1)*2-2,tmpCenter_idx(2)*2-1;
                        tmpCenter_idx(1)*2-1,tmpCenter_idx(2)*2-1;
                        tmpCenter_idx(1)*2,tmpCenter_idx(2)*2-1;
                        tmpCenter_idx(1)*2-2,tmpCenter_idx(2)*2;
                        tmpCenter_idx(1)*2-1,tmpCenter_idx(2)*2;
                        tmpCenter_idx(1)*2,tmpCenter_idx(2)*2];

                    % convert it via
                    % http://stackoverflow.com/questions/792683/compact-matlab-matrix-indexing-notation
                    % to make our array of indices actually index gif_matrix
                    tmp_indexCell = num2cell(tmp_indices,1);
                    tmp_linearIndexMatrix = sub2ind(size(tmp_input),tmp_indexCell{:});
                    tmpArr = tmp_input(tmp_linearIndexMatrix);
                    tmpArr_Input = [double(tmpArr)/256;zeros(auto_encoder.params.size_1-9,1)];
                    
                    if is_Learn_On
                        tmp_indexCell = num2cell(tmp_indices,1);
                        tmp_linearIndexMatrix = sub2ind(size(tmp_target),tmp_indexCell{:});
                        tmpArr = tmp_target(tmp_linearIndexMatrix);
                        tmpArr_Target = [double(tmpArr)/256;zeros(auto_encoder.params.size_3-9,1)];
                    end
                    
                    clear tmpArr tmp_indexCell tmp_linearIndexMatrix tmp_indices tmpCenter_idx
                else
                    % layer_2 from the 4 autoencoders below it
                    
                    % give it the previous time step, target is the current
                    % timestep
%                     tmpArr_Input = [(diff(auto_encoder.layer_2_old{ii-1}{Xcntr*2-1,Ycntr*2-1})./2+.5);
%                                     diff(auto_encoder.layer_2_old{ii-1}{Xcntr*2,Ycntr*2-1}./2+.5);
%                                     diff(auto_encoder.layer_2_old{ii-1}{Xcntr*2-1,Ycntr*2}./2+.5);
%                                     diff(auto_encoder.layer_2_old{ii-1}{Xcntr*2,Ycntr*2}./2+.5)];
%                                     
%                     if is_Learn_On
%                         tmpArr_Target = [diff(auto_encoder.layer_2{ii-1}{Xcntr*2-1,Ycntr*2-1}./2+.5);
%                                         diff(auto_encoder.layer_2{ii-1}{Xcntr*2,Ycntr*2-1}./2+.5);
%                                         diff(auto_encoder.layer_2{ii-1}{Xcntr*2-1,Ycntr*2}./2+.5);
%                                         diff(auto_encoder.layer_2{ii-1}{Xcntr*2,Ycntr*2}./2+.5)];
%                     end
                    tmpArr_Input = [((auto_encoder.layer_2_old{ii-1}{Xcntr*2-1,Ycntr*2-1}));
                                    (auto_encoder.layer_2_old{ii-1}{Xcntr*2,Ycntr*2-1});
                                    (auto_encoder.layer_2_old{ii-1}{Xcntr*2-1,Ycntr*2});
                                    (auto_encoder.layer_2_old{ii-1}{Xcntr*2,Ycntr*2});zeros(auto_encoder.params.size_1-(auto_encoder.params.size_2*4),1)];
                                    
                    if is_Learn_On
                        tmpArr_Target = [(auto_encoder.layer_2{ii-1}{Xcntr*2-1,Ycntr*2-1});
                                        (auto_encoder.layer_2{ii-1}{Xcntr*2,Ycntr*2-1});
                                        (auto_encoder.layer_2{ii-1}{Xcntr*2-1,Ycntr*2});
                                        (auto_encoder.layer_2{ii-1}{Xcntr*2,Ycntr*2});];
                    end
                    %% else its not an input layer
                    % here, we'll use the lower level's activation as input
                    % as well use integral of it, derivative, neighboring
                    % autoencoders too etc.
                    
                end
                
                %% Add to tmpArr_Input the left lateral, right lateral,
                %% upleft feedback and upright feedback context units
                tmp40 = auto_encoder.params.size_3;
                tmp10 = auto_encoder.params.size_2;
                if Xcntr > 1
                    tmpArr_Input(tmp40+1:tmp40+tmp10) = auto_encoder.layer_2{ii}{Xcntr-1,Ycntr};
                elseif Xcntr < auto_encoder.params.auto_encoders_per_layer(ii)
                    tmpArr_Input(tmp40+1+tmp10:tmp40+tmp10*2) = auto_encoder.layer_2{ii}{Xcntr+1,Ycntr};
                end
                
                if Ycntr > 1
                    tmpArr_Input(tmp40+1+tmp10*2:tmp40+tmp10*3) = auto_encoder.layer_2{ii}{Xcntr,Ycntr-1};
                elseif Ycntr < auto_encoder.params.auto_encoders_per_layer(ii)
                    tmpArr_Input(tmp40+1+tmp10*3:tmp40+tmp10*4) = auto_encoder.layer_2{ii}{Xcntr,Ycntr+1};
                end
                    
                if ii < length(auto_encoder.params.auto_encoders_per_layer)-1
                    tmpArr_Input(tmp40+1+tmp10*4:tmp40+tmp10*5) = auto_encoder.layer_2{ii+1}{ceil(Xcntr/2),ceil(Ycntr/2)};
                end
                
%                 if ii < length(auto_encoder.params.auto_encoders_per_layer)
%                     tmpArr_Input(tmp40+1+tmp10*5:tmp40+tmp10*6) = auto_encoder.layer_2{length(auto_encoder.params.auto_encoders_per_layer)}{ceil(Xcntr/2),ceil(Ycntr/2)};
%                 end
                
%                 tmpArr_Input = [tmpArr_Input;...
%                 05-07    auto_encoder.layer_2{ii}{Xcntr+1,Ycntr};...
%                 08-10    auto_encoder.layer_2{ii}{Xcntr,Ycntr+1};...
%                 11-13    auto_encoder.layer_2{ii}{Xcntr-1,Ycntr};...
%                 14-16    auto_encoder.layer_2{ii}{Xcntr,Ycntr-1};...
%                 17-19    auto_encoder.layer_2{ii+1}{Xcntr,Ycntr};...
                
                
                
                %% Forward propogation
                auto_encoder.layer_1{ii}{Xcntr,Ycntr} = tmpArr_Input;
                % remember layer 2 old output
                auto_encoder.layer_2_old{ii}{Xcntr,Ycntr} = auto_encoder.layer_2{ii}{Xcntr,Ycntr};
                
                tmp_netInput_into_2 = auto_encoder.weight_matrix_1_2{ii}{Xcntr,Ycntr} * auto_encoder.layer_1{ii}{Xcntr,Ycntr};
                auto_encoder.layer_2{ii}{Xcntr,Ycntr} = 1.0 ./ (1.0 + exp(-tmp_netInput_into_2));
                
                tmp_netInput_into_3 = auto_encoder.weight_matrix_2_3{ii}{Xcntr,Ycntr} * auto_encoder.layer_2{ii}{Xcntr,Ycntr};
                auto_encoder.layer_3{ii}{Xcntr,Ycntr} = 1.0 ./ (1.0 + exp(- tmp_netInput_into_3));
                
                if is_Learn_On
                    %% Backward propogation
                    error = tmpArr_Target - auto_encoder.layer_3{ii}{Xcntr,Ycntr};
                    gPrime = max(0,min(1,(1.0 ./ (1.0 + exp(- tmp_netInput_into_3)))));
                    gPrime = (gPrime .* (1-gPrime));
                    tmp_unit_delta_3 = error .* gPrime;

                    error = auto_encoder.weight_matrix_2_3{ii}{Xcntr,Ycntr}' * tmp_unit_delta_3;
                    gPrime = max(0,min(1,(1.0 ./ (1.0 + exp(- tmp_netInput_into_2)))));
                    gPrime = (gPrime .* (1-gPrime));
                    tmp_unit_delta_2 = error .* gPrime;

                    %% Change the weights
                    auto_encoder.weight_matrix_2_3{ii}{Xcntr,Ycntr} = auto_encoder.weight_matrix_2_3{ii}{Xcntr,Ycntr} + 0.1 * (tmp_unit_delta_3' * tmp_netInput_into_3);
                    auto_encoder.weight_matrix_1_2{ii}{Xcntr,Ycntr} = auto_encoder.weight_matrix_1_2{ii}{Xcntr,Ycntr} + 0.1 * (tmp_unit_delta_2' * tmp_netInput_into_2);
                else
%                     
%                     if is_Learn_On == -1
%                         %% Backward propogation
%                         error = tmpArr_Target - auto_encoder.layer_3{ii}{Xcntr,Ycntr};
%                         gPrime = max(0,min(1,(1.0 ./ (1.0 + exp(- tmp_netInput_into_3)))));
%                         gPrime = (gPrime .* (1-gPrime));
%                         tmp_unit_delta_3 = -error .* gPrime;
% 
%                         error = auto_encoder.weight_matrix_2_3{ii}{Xcntr,Ycntr}' * tmp_unit_delta_3;
%                         gPrime = max(0,min(1,(1.0 ./ (1.0 + exp(- tmp_netInput_into_2)))));
%                         gPrime = (gPrime .* (1-gPrime));
%                         tmp_unit_delta_2 = -error .* gPrime;
% 
%                         %% Change the weights
%                         auto_encoder.weight_matrix_2_3{ii}{Xcntr,Ycntr} = auto_encoder.weight_matrix_2_3{ii}{Xcntr,Ycntr} + 0.1 * (tmp_unit_delta_3' * tmp_netInput_into_3);
%                         auto_encoder.weight_matrix_1_2{ii}{Xcntr,Ycntr} = auto_encoder.weight_matrix_1_2{ii}{Xcntr,Ycntr} + 0.1 * (tmp_unit_delta_2' * tmp_netInput_into_2);
%                         is_Learn_On = 0;
%                     end
                end
                    
                clear tmp_netInput_into_3 tmp_netInput_into_2 tmp_unit_delta_3 tmp_unit_delta_2 tmpArr_Input tmpArr_Target tmp40 tmp10
            end
        end
    end
% end

% end