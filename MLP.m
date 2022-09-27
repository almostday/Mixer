clc;
clear;
close all;

% Multilayer Perceptron with Backpropagation in Regression Problem.
% Author: Mr.Tanadon Ratsameephen

N = 100;
x = linspace(0, 1, N)';
t = x.^2 + unifrnd(-0.1, 0.1, N, 1);

x = [ones(size(x, 1),1) x];
% t = [1 1 1 -1 -1 -1];

feature_number = size(x,2);

w_hidden_node_number = 10;
w_output_node_number= 1; % regression problem

w_hidden = rand(feature_number + 1, w_hidden_node_number);
w_output = rand(w_hidden_node_number + 1, w_output_node_number);

w = rand(2, 1);
n = 0.05;

E = [];
for i = 1:500
    % Feed forward through neural nets
    y_hidden = tanh([ones(size(x, 1), 1) x] * w_hidden);
    y_output = logsig([ones(size(y_hidden, 1), 1) y_hidden] * w_output);
    
    y = y_output;
    
    e = t - y;
    
    % Backpropagation hidden layer
    delta_w_output = ((e .*(y_output.*(1-y_output))) .* [ones(size(y_hidden, 1), 1) y_hidden]);
    delta_w_hidden = (delta_w_output(:, 2:end) .* ((1-y_hidden.^2)))' * [ones(size(x, 1), 1) x];
    
    w_output = w_output + n * sum(delta_w_output)';
    w_hidden = w_hidden + n * delta_w_hidden';
    
    E(i) = mse(e);
    
    subplot(2, 1, 1)
    plot(x(:, 2), t, '.b', x(:, 2), logsig([ones(N, 1) tanh([ones(N, 1) x] * w_hidden)] * w_output), 'r');
    title('Actual vs. Prediction');
    xlabel('x'); ylabel('y'); legend('Actual', 'Predicted');
    
    subplot(2,1,2);
    plot(E); title('Error vs. Iteration');
    xlabel('Iteration (n)'); ylabel('MSE');

    drawnow;
end