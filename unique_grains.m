function [grains_unique,indx] = unique_grains(grains,criangle)
mis = zeros(size(grains,2));
grains_unique_ind = ones(size(grains));
indx = [];
for i = 1:size(grains,2)
    for j = i+1:size(grains,2)
        mis(i,j)=misori2(grains(i).refined_ori_matrix',grains(j).refined_ori_matrix');
    end
    ind = find(mis(i,:)<criangle);
    
    ind(ind<=i) = [];
    if ~isempty(ind)
        grains_unique_ind(ind) = 0;
        indx = [indx,ind];
    end
end
grains_unique = grains(find(grains_unique_ind==1));