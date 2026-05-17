function pend_value_function_fit_test
close all
clear all
q_bins = linspace(0,2*pi,25);   
qdot_bins = linspace(-10,10,25);
Value = load('PendValueFunction.mat');

%specify these values:
LayerOputSize1 = 64;
LayerOputSize2 = 64;
MaxEpochs = 300;
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

q_bins2 = linspace(0,2*pi,50);   
qdot_bins2 = linspace(-10,10,50);
[q, qdot] = ndgrid(q_bins2,qdot_bins2);
s = [reshape(q,1,numel(q)); reshape(qdot,1,numel(qdot))];

Jest = predict(trainedNet, s');

save ('pendnet.mat', 'trainedNet')

vi_plot(Value.J,q_bins,qdot_bins,10)
vi_plot(Jest,q_bins2,qdot_bins2,11)

end

% ===============================================================
% This function plots the value function and policy
%================================================================
function vi_plot(J,q_bins,qdot_bins, fig_num)
figure(fig_num); n1 = size(q_bins,2); n2 = size(qdot_bins,2);
imagesc(q_bins,qdot_bins,reshape(J,n1,n2)'); axis xy;
drawnow; 
end