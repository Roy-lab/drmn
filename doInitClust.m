% get initial cluster assignments. input: one column gene expression file; list of k for kmeans
% will add k_clusterassign.txt to the end of the file prefix.
function doInitClust(infile, klist, outpref)

    for k=klist
        exprdata=readtable(infile,'Delimiter','\t','ReadVariableNames',true);
        exprmat=table2array(exprdata(:,2:end));
        
        [ cid centers ] = kmeans(exprmat, k, 'Replicates', 10);
        
        % sort centers by expression, low to high
        [ newC ix ] = sort(centers, 'ascend');
        
        % repopulate according to the new IDs
        Genes={};
        ClusterIDs=[];
        % get the genes in order of new IDs -- low to high
        for i=1:numel(ix)
            % i is the new cluster ID, c is the original cluster ID.
            c=ix(i);
            
            % get the rows for this cluster from the ORIGINAL data
            myrows=find(cid==c);
            % rename to new clusters and subtract 1 to start from 0
            ClusterIDs=[ClusterIDs; repmat(i-1, size(myrows))];
            newgenes=exprdata.Gene(myrows);
            Genes=[Genes; newgenes];
        end
        
        size(ClusterIDs)
        size(Genes)
        
        newTable=table(Genes, ClusterIDs);
        
        % print for this cell type
        newFN=sprintf('%s_k%d_clusterassign.txt', outpref, k);
        writetable(newTable, newFN, 'WriteVariableNames',false,'Delimiter','\t');
      
    end % end loop over k
end % end function
