function [J, grad] = costFunction(theta, X, y)
%COSTFUNCTION Compute cost and gradient for logistic regression
%   J = COSTFUNCTION(theta, X, y) computes the cost of using theta as the
%   parameter for logistic regression and the gradient of the cost
%   w.r.t. to the parameters.

% Initialize some useful values
m = length(y); % number of training examples

% You need to return the following variables correctly 
J = 0;
grad = zeros(size(theta));

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost of a particular choice of theta.
%               You should set J to the cost.
%               Compute the partial derivatives and set grad to the partial
%               derivatives of the cost w.r.t. each parameter in theta
%
% Note: grad should have the same dimensions as theta
%

% for i = 1:m
%     g = sigmoid(X(i,:) * theta);
%     J = J + (-y(i) * log(g)) - (1 - y(i)) * log(1 - g);
%     for j = 1:size(X,2)
%         grad(j) = grad(j) + (g - y(i)) * X(i,j);
%     end
% end
% 
% J = J / m;
% grad = grad / m;

% =============================================================
g = sigmoid(X * theta);
tmp_J = -y .* log(g) - (1 - y) .* log(1 - g);
J = sum(tmp_J);
grad = X' * (g - y);

% for i = 1:m
%   %  J = J + (-y(i) * log(g(i)) - (1 - y(i)) * log(1 - g(i)));
%     
%     for j = 1:size(X,2)
%         grad(j) = grad(j) + (g(i) - y(i)) * X(i,j);
%     end
% end
J = J / m;
grad = grad / m;


 end
