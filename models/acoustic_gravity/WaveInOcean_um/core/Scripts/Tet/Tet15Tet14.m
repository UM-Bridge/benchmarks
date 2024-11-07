%%
syms u v w;

%syms z;
z = 1 - u - v - w;

Nnodes = 15;
Nquads = 14;

% Volume.
phi_15 = 256 * z * u * v * w;

% Faces.
phi_14 = 27 * u * v * w * (1 - 4 * z);
phi_13 = 27 * z * v * w * (1 - 4 * u);
phi_12 = 27 * z * u * w * (1 - 4 * v);
phi_11 = 27 * z * u * v * (1 - 4 * w);

% Edges.
phi_10 = 4 * v * w * (1 - 3 * (z * (1 - 4 * u) + u * (1 - 4 * z)) - 16 * z * u);
phi_9 = 4 * u * w * (1 - 3 * (z * (1 - 4 * v) + v * (1 - 4 * z)) - 16 * z * v);
phi_8 = 4 * z * w * (1 - 3 * (u * (1 - 4 * v) + v * (1 - 4 * u)) - 16 * u * v);
phi_7 = 4 * z * v * (1 - 3 * (u * (1 - 4 * w) + w * (1 - 4 * u)) - 16 * u * w);
phi_6 = 4 * u * v * (1 - 3 * (z * (1 - 4 * w) + w * (1 - 4 * z)) - 16 * z * w);
phi_5 = 4 * z * u * (1 - 3 * (v * (1 - 4 * w) + w * (1 - 4 * v)) - 16 * v * w);

% Vertices.
phi_4 = w - 0.5 * (phi_8 + phi_9 + phi_10) - (1 / 3) * (phi_12 + phi_13 + phi_14) - 0.25 * phi_15;
phi_3 = v - 0.5 * (phi_6 + phi_7 + phi_10) - (1 / 3) * (phi_11 + phi_13 + phi_14) - 0.25 * phi_15;
phi_2 = u - 0.5 * (phi_5 + phi_6 + phi_9) - (1 / 3) * (phi_11 + phi_12 + phi_14) - 0.25 * phi_15;
phi_1 = z - 0.5 * (phi_5 + phi_7 + phi_8) - (1 / 3) * (phi_11 + phi_12 + phi_13) - 0.25 * phi_15;

% Node positions.
positions = [[0,0,0];[1,0,0]; [0,1,0]; [0,0,1]; ...
             [0.5,0,0]; [0.5,0.5,0]; [0,0.5,0]; [0,0,0.5]; [0.5,0,0.5]; [0,0.5,0.5]; ...
             [1/3, 1/3, 0]; [1/3, 0, 1/3]; [0, 1/3, 1/3]; [1/3, 1/3, 1/3]; ...
             [0.25, 0.25, 0.25]];
         
shape_functions = [phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, phi_8, phi_9, phi_10, phi_11, phi_12, phi_13, phi_14, phi_15];         


 

 %%
% %Verify that phi_i(x_j) == 1 if i == j, 0 otherwise.
% verif_array = zeros(Nnodes);
%          
% for i = 1:Nnodes
%     for j = 1:Nnodes
%         verif_array(i,j) = subs(shape_functions(i), [u,v,w], positions(j, :));
%     end
% end
% assert(norm(verif_array - eye(Nnodes)) == 0)
% 




syms c1 c2 d;

positions_tet14 = [[c1,c1,c1]; [1-3*c1,c1,c1]; [c1,1-3*c1,c1]; [c1,c1,1-3*c1]; ...
                   [c2,c2,c2]; [1-3*c2,c2,c2]; [c2,1-3*c2,c2]; [c2,c2,1-3*c2]; ...
                   [0.5-d,d,d]; [d,0.5-d,d]; [d,d,0.5-d]; [0.5-d,0.5-d,d]; [d,0.5-d,0.5-d]; [0.5-d,d,0.5-d];];


for j = 1:Nquads
    positions_tet14_bar(j,:) = [positions_tet14(j,:) 1-positions_tet14(j,1)-positions_tet14(j,2)-positions_tet14(j,3)];
end


% Compute derivatives
 
for i = 1:Nnodes

    du_phi = diff(shape_functions(i),u);
    dv_phi = diff(shape_functions(i),v);
    dw_phi = diff(shape_functions(i),w);

    for j = 1:Nquads
        Ru(j,i) = simplify(expand(subs(du_phi, [u,v,w], positions_tet14(j, :))));
        Rv(j,i) = simplify(expand(subs(dv_phi, [u,v,w], positions_tet14(j, :))));
        Rw(j,i) = simplify(expand(subs(dw_phi, [u,v,w], positions_tet14(j, :))));
    end
end



%%
 fileID = fopen('gradientinterpolation1514.txt', 'w');

 
 c1 = 0.09273525031089123;
 c2 = 0.3108859192633006;
 d  = 0.04550370412564965;
 
 for j = 1:15   
     for i = 1:14
        if (Ru(i,j) ~= 0)
            %s = replace(string(expand(Ru(i,j))),'^','_');;
            s = num2str(eval(Ru(i,j)),'%.16f');
            fprintf(fileID,'Buffer = (%s)*in[%i]; \n', s , j-1); 
            fprintf(fileID,'out[%i] += Buffer; \n \n', 3*(i-1));  
        end
        if (Rv(i,j) ~= 0)
            %s = replace(string(expand(Rv(i,j))),'^','_');;
            s = num2str(eval(Rv(i,j)),'%.16f');
            fprintf(fileID,'Buffer = (%s)*in[%i]; \n', s, j-1); 
            fprintf(fileID,'out[%i] += Buffer; \n \n', 3*(i-1)+1);  
        end
        if (Rw(i,j) ~= 0)
            % s = replace(string(expand(Rw(i,j))),'^','_');
            s = num2str(eval(Rw(i,j)),'%.16f');
            fprintf(fileID,'Buffer = (%s)*in[%i]; \n', s , j-1); 
            fprintf(fileID,'out[%i] += Buffer; \n \n', 3*(i-1)+2);  
        end
     end
 end

 %%
 fileID = fopen('transposegradientinterpolation1514.txt', 'w');

 for i = 1:14   

     for j = 1:15
        if (Ru(i,j) ~= 0)
            % s = replace(string(expand(Ru(i,j))),'^','_');
            s = num2str(eval(Ru(i,j)),'%.16f');
            fprintf(fileID,'Buffer = (%s)*in[%i]; \n',s, 3*(i-1)); 
            fprintf(fileID,'out[%i] += Buffer; \n \n', (j-1));  
        end
     end

     for j = 1:15
        if (Rv(i,j) ~= 0)
            % s = replace(string(expand(Rv(i,j))),'^','_');
            s = num2str(eval(Rv(i,j)),'%.16f');
            fprintf(fileID,'Buffer = (%s)*in[%i]; \n', s, 3*(i-1)+1); 
            fprintf(fileID,'out[%i] += Buffer; \n \n', (j-1));  
        end
     end

     for j = 1:15
        if (Rw(i,j) ~= 0)
            % s = replace(string(expand(Rw(i,j))),'^','_');
            s = num2str(eval(Rw(i,j)),'%.16f');
            fprintf(fileID,'Buffer = (%s)*in[%i]; \n', s, 3*(i-1)+2); 
            fprintf(fileID,'out[%i] += Buffer; \n \n', (j-1));  
        end
     end
     
 end




%  %% Generate a more efficient way to compute the matrix vector product
% 
%  c1_value = 0.09273525031089123;
%  c2_value = 0.3108859192633006;
%  d_value  = 0.04550370412564965;
%  
% num_values = [c1_value c2_value d_value];
% 
%  fileID = fopen('gradientinterpolation1514_opt.txt', 'w');
% 
%  for j = 1:15
% 
%      
%      for i = 1:14
%          Value(3*(i-1)+1) = Ru(i,j);
%          Value(3*(i-1)+2) = Rv(i,j);
%          Value(3*(i-1)+3) = Rw(i,j);
%      end
% 
%      Indices_to_treat = 1:(3*14);
% 
% 
%      while (size(Indices_to_treat,2) > 0)
% 
%          curr_index = Indices_to_treat(1);
% 
%          if (Value(curr_index) == 0)
%              Indices_to_treat(1) = [];
%          else 
%          
%             s = num2str(eval(subs(Value(curr_index),[c1 c2 d], num_values)),'%.16f');
%             fprintf(fileID,'\nBuffer = (%s)*in[%i]; \n', s , j-1); 
%             fprintf(fileID,'out[%i] += Buffer; \n', curr_index-1);  
%                 
%             Indices_to_treat(1) = [];
% 
%             for k=1:size(Indices_to_treat,2)
% 
%                 curr_index_copy = Indices_to_treat(k);
% 
%                 if (Value(curr_index) == Value(curr_index_copy))
% 
%                     fprintf(fileID,'out[%i] += Buffer; \n', curr_index_copy-1);  
%                     Indices_to_treat(k) = -1;
% 
%                 else
%                 
%                     if (Value(curr_index) == -Value(curr_index_copy))
% 
%                         fprintf(fileID,'out[%i] -= Buffer; \n', curr_index_copy-1);  
%                         Indices_to_treat(k) = -1;
%                 
%                     end
%                 end
%             end
% 
%             Indices_to_treat(Indices_to_treat == -1) = [];
% 
%          end
% 
%      end
% 
%   
%      
%      end
%  







 






