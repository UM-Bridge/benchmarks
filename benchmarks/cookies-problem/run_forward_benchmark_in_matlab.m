%% ----------------------- forward UQ benchmark --------------------------
% 
% -------------------------------------------------------------------------


%% preliminar operations:

% 1) add sparse grids matlab kit and um-bridge matlab client to path
addpath(genpath('~/GIT_projects/Github/sparse-grids-matlab-kit/')) % this is version 23.5 Robert
addpath(genpath('~/GIT_projects/Github/umbridge/matlab/'))

% 2) run docker container with a command like
% sudo docker run -it -p 4242:4242 <image name>

clear 
 
%% setup model 

% specify model port
uri = 'http://0.0.0.0:4242';
model = HTTPModel(uri,'forward');

% config cookie solver
config = struct('NumThreads',4,'BasisDegree',3,'Fidelity',1);


% wrap model in an @-function too
Psi_fun = @(y) model.evaluate(y',config);


%% create sparse grid

% setup sparse grid ingredients. We will build increasingly large grids using a for loop

N = 8;
knots = @(n) knots_CC(n,-0.99,-0.2);
lev2knots = @lev2knots_doubling;
idxset_rule = @(i) sum(i-1);
idxset_level_max = 2; % <--- controls max size of sparse grid


% With the choice above, the sparse grids will be nested. To recycle evaluations from one grid 
% to the next, we need some containers

S_old         = []; % the previous sparse grid in extended format
Sr_old        = []; % the previous sparse grid in reduced format
Psi_evals_old = []; % the evaluations of the model output on the previous sparse grid

% more containers, to save results at each iteration

nb_pts = []; % nb of points of each sparse grid
Psi_EV = []; % the expected value of Psi computed on each grid

% here we go with the loop
for w = 0:idxset_level_max

    % create sparse grid in extended format. Recycle some work from previous sparse grid
    S = create_sparse_grid(N,w,knots,lev2knots,idxset_rule,S_old);

    % reduce sparse grid
    Sr = reduce_sparse_grid(S);
    
    % eval Psi on each point of sparse grid. Recycle available evaluations whenever possible
    Psi_evals = evaluate_on_sparse_grid(Psi_fun,S,Sr,Psi_evals_old,S_old,Sr_old);
    
    % compute expected value and append the new value to previous container. Same for nb_pts
    Psi_EV(end+1) = quadrature_on_sparse_grid(Psi_evals,Sr);    %#ok<SAGROW>
    nb_pts(end+1) = Sr.size;                                    %#ok<SAGROW>
    
    % update containers
    S_old = S;
    Sr_old = Sr;
    Psi_evals_old = Psi_evals;
    
end

%% plot convergence

figure
semilogx(nb_pts,Psi_EV,'-ok','LineWidth',2,'MarkerFaceColor','k')

