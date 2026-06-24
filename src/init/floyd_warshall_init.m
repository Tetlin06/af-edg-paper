function [Dist_fw,pt_fw,notinf_idx] = floyd_warshall_init(Dist_input,weight_input,num_eigen)
    % number of points
    num_pts = length(Dist_input);
    Dist_sq= sqrt(Dist_input);
    Dist = Dist_sq.*Dist_sq;
    Distw= sqrt(Dist.*weight_input);
    Distw(Distw==0)=Inf;
    Dist_fw = r8mat_floyd(Distw);
    % Find rows with Inf and remove them
    inf_idx = [];
    i_list = 1:num_pts;
    j_list = 1:num_pts;
    for i = 1:length(i_list)
        for j = 1:length(j_list)
           if Dist_fw(i,j)==Inf
               inf_idx = [inf_idx j];              
           end
           j_list = sort(setdiff(1:num_pts,inf_idx));
        end
        i_list = sort(j_list);
    end
    % indices not infinity
    notinf_idx = setdiff(1:num_pts,inf_idx);
    %notinf_idx = 1:num_pts;
    ni = length(notinf_idx);
    Dist_fw = Dist_fw(notinf_idx(1:ni),notinf_idx(1:ni));
    Dist_fw = Dist_fw.*Dist_fw;
    Dist_fw = Dist_fw - diag(diag(Dist_fw));
    IPM_est = Dist_fw - mean(Dist_fw,2)*ones(1,ni);
    IPM_est =  - 1/2*(IPM_est - ones(ni,1)*mean(IPM_est,1));
    [V, D] = eigs(IPM_est,num_eigen,'lr');
    D = diag(D);
    [D,IJ] = sort(D,'descend');
    V = V(:,IJ);
    pt_fw = real(V(:,1:num_eigen)*diag(sqrt(D(1:num_eigen))));
end




