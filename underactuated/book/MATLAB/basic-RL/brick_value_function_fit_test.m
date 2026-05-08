function brick_value_function_fit_test
close all
clear all
 q_bins = [-3:.2:3];   
 qdot_bins = [-4:.2:4];
Value = load('BrickValueFunction.mat');

%specify these values:
LayerOputSize1 = 20;
LayerOputSize2 = 20;
MaxEpochs = 200;
MiniBatchSize = 30;

%% Define a layer architecture
layers = [ ...
    featureInputLayer(2, "Name", "myFeatureInputLayer")
    fullyConnectedLayer(LayerOputSize1, "Name", "myFullyConnectedLayer1")
    reluLayer("Name", "myReLu1")
    fullyConnectedLayer(LayerOputSize2, "Name", "myFullyConnectedLayer2")
    reluLayer("Name", "myReLu2")
    fullyConnectedLayer(1, "Name", "myFullyConnectedLayer3")
    regressionLayer("Name", "myRegressionLayer")
];
%% Define options for the training
opts = trainingOptions('adam', ...
    'MaxEpochs',MaxEpochs, ...
    'InitialLearnRate',0.01,...
    'MiniBatchSize',MiniBatchSize, ...
    'Verbose',false, ...
    Plots="training-progress");
%% Train the network
[trainedNet, info] = trainNetwork(Value.s', Value.J, layers, opts);

 q_bins2 = [-3:.02:3];   
 qdot_bins2 = [-4:.02:4];
[q, qdot] = ndgrid(q_bins2,qdot_bins2);
s = [reshape(q,1,numel(q)); reshape(qdot,1,numel(qdot))];


y = predict(trainedNet, s')

save ('bricknet.mat', 'trainedNet')


%  [x,t] = simplefit_dataset;
%  net = fitnet([10,10]);
%  net = train(net,Value.s,Value.J);
%  view(net)
%  y = net(x);

vi_plot(Value.J,[],q_bins,qdot_bins,10)

vi_plot(y,[],q_bins2,qdot_bins2,11)

end

% ===============================================================
% This function plots the value function and policy
%================================================================
function vi_plot(J,PI,q_bins,qdot_bins, fig_num)
figure(fig_num); n1 = size(q_bins,2); n2 = size(qdot_bins,2);
%subplot(2,1,1);imagesc(q_bins,qdot_bins,reshape(PI,n1,n2)'); axis xy;
imagesc(q_bins,qdot_bins,reshape(J,n1,n2)'); axis xy;
drawnow; 
end