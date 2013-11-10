function [J, grad] = linearRegCostFunction(X, y, theta, lambda)
%LINEARREGCOSTFUNCTION Compute cost and gradient for regularized linear 
%regression with multiple variables
%   [J, grad] = LINEARREGCOSTFUNCTION(X, y, theta, lambda) computes the 
%   cost of using theta as the parameter for linear regression to fit the 
%   data points in X and y. Returns the cost in J and the gradient in grad

% Initialize some useful values
m = length(y); % number of training examples

% You need to return the following variables correctly 
J = 0;
grad = zeros(size(theta));

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost and gradient of regularized linear 
%               regression for a particular choice of theta.
%
%               You should set J to the cost and grad to the gradient.
%

n = length(theta);

delta = X * theta - y;

J = sum(delta.^ 2);

reg = 0; 
for i = 2:n
    reg = reg + theta(i) ^2;
end

reg = lambda * reg;

J = (J + reg) / (2 * m);

for i = 1:m
    for j = 1:n
        grad(j) = grad(j) + delta(i) * X(i,j);
    end
end

for i = 2:n
    grad(i) = grad(i) + lambda * theta(i);
end

grad = grad / m;


% =========================================================================

grad = grad(:);

end
