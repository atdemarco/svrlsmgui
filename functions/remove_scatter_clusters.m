function [out_map] = remove_scatter_clusters(in_map, ClusterThreshold)
    out_map = zeros(size(in_map));
    CC = bwconncomp(in_map, 6);
    count = 0;
    for ti=1:CC.NumObjects
        ObjIdx = CC.PixelIdxList{1,ti};
        if(size(ObjIdx,1) < ClusterThreshold)
            continue;
        else
            out_map(CC.PixelIdxList{1,ti}) = in_map(CC.PixelIdxList{1,ti});
        end
    end
end