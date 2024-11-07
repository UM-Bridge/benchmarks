%% Enriched P3
syms u v w z;


 

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
Nnodes = size(positions,1);


% Print shape functions.
fileID = fopen('tetra15.txt','w');
for i = 1:Nnodes
    fprintf(fileID, 'Phi %i: %s\n', i, string(simplify(shape_functions(i))));
end
fclose(fileID);


% Print derivative shape functions.
fileID = fopen('tetra15_d1.txt','w');
for i = 1:Nnodes
    fprintf(fileID, 'dPhi_1 %i: %s\n', i, string(simplify(diff(shape_functions(i),u)-diff(shape_functions(i),z))));
end
fclose(fileID);

% Print derivative shape functions.
fileID = fopen('tetra15_d2.txt','w');
for i = 1:Nnodes
    fprintf(fileID, 'dPhi_2 %i: %s\n', i, string(simplify(diff(shape_functions(i),v)-diff(shape_functions(i),z))));
end
fclose(fileID);

% Print derivative shape functions.
fileID = fopen('tetra15_d3.txt','w');
for i = 1:Nnodes
    fprintf(fileID, 'dPhi_3 %i: %s\n', i, string(simplify(diff(shape_functions(i),w)-diff(shape_functions(i),z))));
end
fclose(fileID);

% Goes to uvw coordinate only
for i = 1:Nnodes
    shape_functions(i) =    subs(shape_functions(i),z,1 - u - v - w);
end


 %%
%Verify that phi_i(x_j) == 1 if i == j, 0 otherwise.
verif_array = zeros(Nnodes);
         
for i = 1:Nnodes
    for j = 1:Nnodes
        verif_array(i,j) = subs(shape_functions(i), [u,v,w], positions(j, :));
    end
end
assert(norm(verif_array - eye(Nnodes)) == 0)


%%
% Compute derivatives
 

for i = 1:Nnodes

    du_phi = diff(shape_functions(i),u);
    dv_phi = diff(shape_functions(i),v);
    dw_phi = diff(shape_functions(i),w);

    for j = 1:Nnodes
        Rs(j,i) = subs(du_phi, [u,v,w], positions(j, :));
        Rr(j,i) = subs(dv_phi, [u,v,w], positions(j, :));
        Rt(j,i) = subs(dw_phi, [u,v,w], positions(j, :));
    end
end

nnz(Rs)/(size(Rs,1)*size(Rs,2))

%  RsScaled = sym(diag([1 1 1 1 1 4 4 4 4 4 3 3 27 27 16]))*Rs*sym(diag([1 1 1/3 1/3 1/4 1/4 1/4 1/4 1/4 1/32 1/27 1/27 1/27 1/27 1/256])) 
%  RrScaled = sym(diag([1 1 1 1 4 4 1 4 4 4 3 27 3 27 16]))*Rr*sym(diag([1 1/3 1 1/3 1/4 1/4 1/4 1/4 1/32 1/4 1/27 1/27 1/27 1/27 1/256])) 
%  RtScaled = sym(diag([1 1 1 1 4 4 4 1 4 4 27 3 3 27 16]))*Rt*sym(diag([1 1/3 1/3 1 1/4 1/32 1/4 1/4 1/4 1/4 1/27 1/27 1/27 1/27 1/256])) 
% % 
% %
% size(nonzeros(RsScaled(abs(RsScaled)~=1)))
% size(nonzeros(RrScaled(abs(RrScaled)~=1)))
% size(nonzeros(RtScaled(abs(RtScaled)~=1)))


%New method : same scaling is applied

 RsScaled = sym(diag([1 1 1 1 4 4 4 4 4 4 27 27 27 27 16]))*Rs*sym(diag([1 1 1 1 1/4 1/4 1/4 1/4 1/4 1/4 1/27 1/27 1/27 1/27 1/16])) 
 RrScaled = sym(diag([1 1 1 1 4 4 4 4 4 4 27 27 27 27 16]))*Rr*sym(diag([1 1 1 1 1/4 1/4 1/4 1/4 1/4 1/4 1/27 1/27 1/27 1/27 1/16])) 
 RtScaled = sym(diag([1 1 1 1 4 4 4 4 4 4 27 27 27 27 16]))*Rt*sym(diag([1 1 1 1 1/4 1/4 1/4 1/4 1/4 1/4 1/27 1/27 1/27 1/27 1/16])) 
 
PrintScaling('RsScaled', RsScaled);
PrintScaling('RrScaled', RrScaled);
PrintScaling('RtScaled', RtScaled);
% 
%
size(nonzeros(RsScaled(abs(RsScaled)~=1)))
size(nonzeros(RrScaled(abs(RrScaled)~=1)))
size(nonzeros(RtScaled(abs(RtScaled)~=1)))

%%Check
assert(norm(Rs -  sym(diag([1 1 1 1 1/4 1/4 1/4 1/4 1/4 1/4 1/27 1/27 1/27 1/27 1/16])) *RsScaled*sym(diag([1 1 1 1 4 4 4 4 4 4 27 27 27 27 16])))==0)
assert(norm(Rr -  sym(diag([1 1 1 1 1/4 1/4 1/4 1/4 1/4 1/4 1/27 1/27 1/27 1/27 1/16])) *RrScaled*sym(diag([1 1 1 1 4 4 4 4 4 4 27 27 27 27 16])))==0)
assert(norm(Rt -  sym(diag([1 1 1 1 1/4 1/4 1/4 1/4 1/4 1/4 1/27 1/27 1/27 1/27 1/16])) *RtScaled*sym(diag([1 1 1 1 4 4 4 4 4 4 27 27 27 27 16])))==0)


function PrintScaling(filename, scaling_matrix)
    fileID = fopen([filename '.txt'], 'w');
    for i = 1:size(scaling_matrix, 1)
        fprintf(fileID,'out[%i] =', i-1);    
        for j = 1:size(scaling_matrix, 2)
            value = GetPrintValue(scaling_matrix(i,j));
            if (scaling_matrix(i,j) ~= 0)
                if (j == 1 && value == '+')
                    fprintf(fileID,' in[%i]', j-1);
                else
                    fprintf(fileID,' %s in[%i]', value, j-1);
                end
            end
        end
        fprintf(fileID, ';\n');  
    end
end

function ret = GetPrintValue(value)
    if (value == 0)
        ret = '';
    end
    if (abs(value) ~= 1)
        ret = string(value) + " *";
        if (value > 0)
            ret = "+ " + ret;
        end
    end
    if (value == 1)
        ret = '+';
    end
    if (value == -1)
        ret = '-';
    end
end


