
#Answer to Question 1. goes in this Markdown cell.

res = get("http://pages.cs.wisc.edu/~okan/cs412/assignment3/X.txt")
X = readdlm(IOBuffer(res.data))
X = X[:, 2:end]
display(X)
d = size(X, 2)
n = size(X,1)
res = get("http://pages.cs.wisc.edu/~okan/cs412/assignment3/Y.txt")
Y = readdlm(IOBuffer(res.data))
b = zeros(2*n,1)
for i= 1:n
b[i,1] = 1
end
display(b)
A = zeros(2*n,d+1+n)
for i= 1:n
for j = 1:d
A[i,j] = Y[i].*X[i,j]
end
A[i,d+1] = Y[i]
for k = d+2: d+1+n
if(k-(d+1) == i)
A[i,k]= 1
else
A[i,k]= 0
end
end
end
for g = n+1:2*n
for p = d+2: d+1+n
if(p-(d+1) == g-n)
A[g,p]= 1
else
A[g,p]= 0
end 
end
end

Answer to Question 3.a goes in this Markdown cell.
##a.
a=zeros(2n,1)
lamda

using PyPlot

nIter = 200

# Generate random data
A = rand(10, 10)

y = rand(10)

b = A * y

# Pick initial point
# When solving a linear system with gradient descent
# you can start from any point
# Convergence is guaranteed as long as alpha > cond(A'*A)
u=[zeros(n,1) ones(n+1,1)]
alpha = 20

# Best parameter possible
# alpha = norm(A'*A) # as challenging as solving the linear system

print("Alpha = $alpha\n")
x_history = zeros(10, nIter)
x_err_history = zeros(nIter, 1)
b_err_history = zeros(nIter, 1)

for i = 1:nIter
    # Store interesting data for plots
    x_history[:, i] = x
    x_err_history[i] = norm(x-y)
    b_err_history[i]= norm(b-A*x)
    
    # This is another possible update rule, it is slower but converges
    # regardless of the value of beta
    # Numerical instabilities are still possible
    # Picking a good beta value increases computational efficiency
    # beta = 1 
    # x = x - (beta/i) * A' * (A * x - b)
    x = x - (1/alpha) * A' * (A * x - b)
end

err_norm = norm(y - x)
y_norm = norm(y)
b_norm = norm(b)
b_err = norm(A*x - b)
cond_A = cond(A' * A)

subplot(2,1,1)
plot(1:nIter, x_err_history)
title("Plot of \$ \\|x_k - y\\|_2\$")

subplot(2,1,2)
plot(1:nIter, b_err_history)
title("Plot of \$\|Ax_k -b\|_2\$")

print("cond(A'A): $cond_A\nNorm of solution: $y_norm\nNorm of Error (x-y): " *
    "$err_norm\nNorm of b: $b_norm\nNorm of b - Ax: $b_err")

Answer to Question 4.a goes in this Markdown cell.

# You can use this outline as the basis for your implementation
# Make sure that your implementation is readible
# You should also order matrix multiplication to minimize the computational complexity
# eg. Choose (A(Bx)) over (AB)x since the first one requires 2n^2 multiplication whereas
# the second one requires approximately n^3 + n^2 multiplications

using Requests

# If verbose is 1, your function should create a 2x1 subplot
# plot at 1,1 should show the optimal objective value for every t
# plot 1,2 should show the accuracy on X for every t 
# You should also output the total time spent in the function
function createDiag(x)
sz = size(x, 1) 
diag = zeros(sz, sz)
for i = 1:sz
diag[i, i] = x[i]
end
return diag
end
function solveSVM(X, Y, λ; newtonMethod=0, verbose=0)
# Number of data points
n = size(X, 1)
# Dimension of data points
d = size(X, 2)
#20 Set the initial point
u = 2*ones(n+1+d,1);
‪#‎for‬ i = d+2:n+d+1
# u[i,1]=2;
‪#‎end‬
‪#‎H‬ = inv(laplacianF(u))
# Set a value for the previous point
u_pr = zeros(n+1+d,1);
# Construct the constraint matrix A and vector b
A = [diagm(Y[:,1])*X Y eye(n);zeros(n,d+1) eye(n)]
‪#‎println‬(A);
b = [ones(n,1) ; zeros(n, 1)]
#println(b);
# Start with 1/t = 16
t = 1/20;
# Set the current solution $x^*_{t_k}$ and the previous
# solution $x^*_{t_{k-1}}$ to their initial values
sol = u;
sol_pr = u_pr;
firstRun = 1;
#40 
println(u);
while (norm(sol-sol_pr)/norm(sol)>0.0004) # convergence criterion 
#println("outer loop, convergence not yet achieved");
# Current solution becomes the old solution
sol_pr = sol;
u_pr = u;
# start the algorithm from the previous solution
# (warm-starting)
u = sol
# Solve the barrier penalized minimization problem
# with Newton's method
while( (norm(u-u_pr)/norm(u)>0.0004) || (firstRun == 1) ) # convergence criterion
firstRun = 0;
#println("inner loop, convergence not yet achieved error: ");
#println(norm(u-u_pr)/norm(u));
u_pr = u;
# Construct s and S
s = A * u - b;
s_inv = s.^-1;
s_inv_2 = s.^-2; 
#println(s_inv_2);
S = createDiag(s_inv_2);
omega = [u[1:d]; 0; λ*ones(n, 1)]
deltaF = -1/t*A'*s_inv+omega;
#println(size(deltaF));
laplacianF = 1/t*A'*S*A + [eye(d) zeros(d,n+1); zeros(n+1, d) zeros(n+1,n+1)];
#println(size(laplacianF));
# Solve the Newton equation with backslash operator
if(newtonMethod == 0)
y = deltaF\laplacianF;
#println(size(y));
#println("In backslash");
# Solve the Newton method with GD
elseif(newtonMethod == 1) 
println("In GD");
# while(norm(deltaF-(laplacianF*y))/norm(deltaF)>.001)
# y = y - laplacianF*(laplacianF*y-deltaF);
# end
# Solve the Newton equation with SR1
return sol[1:(d+1)];
elseif(newtonMethod == 2)
println("In SR1");
‪#‎s‬ = u-u_prev;
‪#‎y‬ = deltaF(u)-deltaF(u_prev);
#H = H +(s-H*y)(s-H*y)'./(s-H*y)'y;
#y = u - H*deltaF;
return sol[1:(d+1)];
end
u = u - y';
end
sol = u;
# Set the next value of t
t = 2 * t;
end
# First d+1 values correspond to w and c
return sol[1:(d+1)];
end
# Set the random seed
# Do not change the following line
srand(1)
# Load the datasets
# Every row of X is a datapoint 
# computed from a breast mass image
res = get("http://pages.cs.wisc.edu/~okan/cs412/assignment3/X.txt")
X = readdlm(IOBuffer(res.data))
# Skip the column of user IDs
X = X[:, 2:end]
# Every entry of Y corresponds to the class of corresponding
# row in X
# We have two classes: benign and cancerous tumors
# represented by 1 and -1 values
res = get("http://pages.cs.wisc.edu/~okan/cs412/assignment3/Y.txt")
Y = readdlm(IOBuffer(res.data))
# Number of data points
n = size(X, 1)
# Dimension of data points
d = size(X, 2)
# Size of training, validation and test sets
n_tr = 400
n_val = 50 
n_test = 233
# Partition the dataset into training, validation and test sets
X_tr = X[1:n_tr, smile emoticon
X_val = X[(n_tr+1):(n_tr+n_val), smile emoticon
X_test = X[(n_tr+n_val+1):(n_tr+n_val+n_test), smile emoticon
Y_tr = Y[1:n_tr, smile emoticon
Y_val = Y[(n_tr+1):(n_tr+n_val), smile emoticon
Y_test = Y[(n_tr+n_val+1):(n_tr+n_val+n_test), smile emoticon
# The vector of lambda values you want to try
#0.1 1 10 100 1000
lambdaVals = [0.1];#, 1,10,100,1000]
# Solve the SVM problem over a set of lambda values and 
# pick the solution that gives the best prediction on the
# validation set
best_accuracy = 0 
for lambda in lambdaVals
r = solveSVM(X_tr, Y_tr, 1, newtonMethod=0)
println(r);
w = r[1:d]
c = r[d+1]
Y_est = X_val * w + c
Y_est[Y_est .> 0] = 1
Y_est[Y_est .<= 0] = -1
println(Y_est);
current_accuracy = sum(Y_val .== Y_est) /size(Y_val, 1)
println(current_accuracy);
if current_accuracy > best_accuracy
best_w = w
best_accuracy = current_accuracy
end
end

# Report your results on the test set
...

# Include comments about the performance of individual algorithms
# How did you pick parameters like lambda, alpha, t?
# How did you pick the initial points?
# Which method seems to be the best for Wisconsin Breast Cancer dataset? Why?
