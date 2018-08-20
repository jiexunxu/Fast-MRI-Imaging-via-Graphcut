function num_edges = write_graph(grstr, conn_mtx, prune_flg, nodetypes)

if nargin < 3,
    prune_flg = 0;
end
szc = size(conn_mtx);
nnodes = szc(1);
if nnodes ~= szc(2)
    disp('connectivity matrix is not square - check!!');
end
if nargin < 4 || isempty(nodetypes),
    nodetypes = ones(nnodes,1);
end

if length(szc) > 2
    nsubject = szc(3);
else
    nsubject = 1;
end

num_edges = zeros(nsubject, 1);

%added new to prune matrix to remove zero rows/cols
conn_mtx_pruned = cell(nsubject);

for ind = 1:nsubject
    if prune_flg
        I = find(sum(abs(conn_mtx(:,:,ind))) > 0);
        conn_mtx_pruned{ind} = conn_mtx(I,I,ind);
    else
        conn_mtx_pruned{ind} = conn_mtx(:,:,ind);
    end
end

% szc = size(conn_mtx),
% szc_pruned = size(conn_mtx_pruned{ind}),
nnodes = szc(1);
fgr = fopen(grstr,'wt');
if fgr == -1
    printf('Cannot open graph file %s for writing\n', grstr);
else
    for ind = 1:nsubject
        fprintf(fgr,'#label\n%s_subj%d\n', grstr, ind);
        %pruned_size = size(conn_mtx_pruned{ind}, 1),
        nnodes = size(conn_mtx_pruned{ind}, 1);
        if nnodes>1
            fprintf(fgr,'#types\n');
            for i=1:nnodes
                fprintf(fgr,'\t%d', nodetypes(i));	% all nodes get type 1 (for now...)
            end
            fprintf(fgr,'\n');
            fprintf(fgr,'#edges\n');
            for i=1:nnodes-1
              for j = i+1:nnodes
%                    szij = [ind, nnodes, i,j],
                    if conn_mtx_pruned{ind}(i,j) > 0
                        %disp('found at least one edge!');
                        num_edges(ind) = num_edges(ind) + 1;
                        %disp(sprintf('%d %d %d\n', i-1, j-1, round(conn_mtx_pruned{ind}(i,j)))); 
                        fprintf(fgr,'%d %d %d\n', i-1, j-1, round(conn_mtx_pruned{ind}(i,j))); 
                    end
                end
            end
            if num_edges(ind) == 0
                fprintf(fgr,'%d %d %d\n', 0, 1, 1);
            end
        else
            disp(sprintf('found empty graph at %d-th subject', ind));
            % write empty but valid graph file
            fprintf(fgr,'#types\n');
            fprintf(fgr,'\t%d\t%d\n', 1, 1);	% all nodes get type 1 (for now...)
            fprintf(fgr,'#edges\n');
            %disp(sprintf('%d %d %d\n', 0, 1, 0)); 
            fprintf(fgr,'%d %d %d\n', 0, 1, 1); 
        end
    end
    fclose( fgr );
end
