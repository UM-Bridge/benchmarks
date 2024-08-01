%% ----------------------- forward UQ benchmark --------------------------
% 
% -------------------------------------------------------------------------


%% preliminar operations:

clear 

% 1) run docker container with a command like
% sudo docker run -it -p 4242:4242 <image name>


% 2) add sparse grids matlab kit and um-bridge matlab client to path
% addpath(genpath('~/GIT_projects/Github/sparse-grids-matlab-kit/')) % this is version 23.5 Robert
% addpath(genpath('~/GIT_projects/Github/umbridge/matlab/'))
addpath(genpath('../../../sparse-grids-matlab-kit/'))
addpath(genpath('../../../umbridge'))

% 3) to save results on file (sparse grids in .txt / .mat, results in .txt / .mat), set the flag below to true
saving_stuff = true;


% 4) to plot results, set the flag below to true
plotting = true;


%% setup model 

% specify model port
uri = 'http://0.0.0.0:4242';
model = HTTPModel(uri,'benchmark');

% config cookie solver with num threades
config = struct();


% wrap model in an @-function too
Psi_fun = @(y) model.evaluate(y',config);

% a simple call to test that things are working fine
% y_test = [-0.5; -0.5; -0.5; -0.5; -0.5; -0.5; -0.5; -0.5;];
% Psi_fun(y_test)

%% create sparse grid

% setup sparse grid ingredients. Instead of building the final grid in one go, we build it gradaully using a for
% loop that creates a sequence of nested grids

N = 8;
knots = @(n) knots_CC(n,-0.99,-0.2);
lev2knots = @lev2knots_doubling;
idxset_rule = @(i) sum(i-1);
idxset_level_max = 5; % <--- controls max size of sparse grid


% To recycle evaluations from one grid to the next, we need some containers

S_old         = []; % the previous sparse grid in extended format
Sr_old        = []; % the previous sparse grid in reduced format
Psi_evals_old = []; % the evaluations of the model output on the previous sparse grid

% more containers, to save results at each iteration

nb_pts = []; % nb of points of each sparse grid
Psi_EV = []; % the expected value of Psi computed on each grid

% here we go with the loop
for w = 0:idxset_level_max

    disp('========================================')
    disp(strcat('w=',num2str(w)))
    disp('========================================')

    % create sparse grid in extended format. Recycle some work from previous sparse grid
    S = create_sparse_grid(N,w,knots,lev2knots,idxset_rule,S_old);

    % reduce sparse grid
    Sr = reduce_sparse_grid(S);
        
    % eval Psi on each point of sparse grid. Recycle available evaluations whenever possible
    Psi_evals = evaluate_on_sparse_grid(Psi_fun,S,Sr,Psi_evals_old,S_old,Sr_old);
    
    if saving_stuff
        % save it to file. Type help export_sparse_grid_to_file for info on the saving format
        grid_filename = strcat('sparse_grid_w=',num2str(w),'.txt');
        export_sparse_grid_to_file(Sr,grid_filename,'with_weights');
        % save also some minimal information in .mat format
        grid_filename_mat = grid_filename(1:end-4); % i.e., remove .txt
        save(grid_filename_mat,'S','Sr','Psi_evals')
    end
    
    % compute expected value and append the new value to previous container. Same for nb_pts
    Psi_EV(end+1) = quadrature_on_sparse_grid(Psi_evals,Sr);    %#ok<SAGROW>
    nb_pts(end+1) = Sr.size;                                    %#ok<SAGROW>
    
    % update containers
    S_old = S;
    Sr_old = Sr;
    Psi_evals_old = Psi_evals;
    
end


if saving_stuff
    % save values to be plotted on file, in txt format. The data saved consists of the two vectors nb_pts and Psi_EV,
    % they will be stored as rows of the txt file
    save('cookies-benchmark-output.txt','nb_pts','Psi_EV','-ascii', '-double')
    % save also on .mat
    save('cookies-benchmark-output','nb_pts','Psi_EV')
end


if plotting
    % in principle, only nb-pts and Psi_EV are needed, so if you have them on file already you can just load results 
    % load('cookies-benchmark-output')
    figure
    semilogx(nb_pts,Psi_EV,'-ok','LineWidth',2,'MarkerFaceColor','k','DisplayName','sparse grid approx. of $\mathbf{E}[\Psi]$')
    grid on
    xlabel('sparse grids points')
    ylabel('$\mathbf{E}[\Psi]$','interpreter','latex','rotation',0)
    legend show
    set(legend,'interpreter','latex','location','northwest','fontsize',14)
end