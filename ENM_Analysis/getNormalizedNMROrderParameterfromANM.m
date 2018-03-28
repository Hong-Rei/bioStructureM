function [order_parameter] = getNormalizedNMROrderParameterfromANM(PDB_Structure,eigvalues,mode_selection)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the ANM theroretical predicted NMR NOE order parameter profile of a specified bond vector.
% Here you can specify how many modes you want to include in this computation, normally the slowest modes are excluded one step each.
% 
% The individual terms in the order parameter before squaring them were normalized by the sum of xijxij, yijyij and zijzij.
% 
% input:
%   eigvalues are ANM eigenvalues
%   PDB_Structure are Structure obtained from ANM.m with ANM fields.
%   IMPORTANT: The input eigvalues and eigvectors should have smallest variance mode at top (MATLAB default!!) 
%               and the first six zero-modes should be already eliminated.
%   mode_selection is an array that specify the index of modes you want to use to reform covariance.
%   
% return:
%   order_parameter: order parameter(S^2) profile calculated theoretically from ANM for each residue.
%   It is a Nx1 column vector. 
%
% Editor: Hong-RUi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    atom_num = length(PDB_Structure);
    SpringConst = 100;

    if ~exist('mode_selection','var')
        mode_selection = 1:(atom_num*3-6);
    end

    eigvalues = eigvalues(mode_selection);
    eigvalues = diag(eigvalues);
    eigvectors = getANM(PDB_Structure,mode_selection);

    reduced_hessian = eigvectors*inv(eigvalues)*(eigvectors');
    reduced_covariance = 1.38064852*10^(-23)*300/SpringConst*10^(20)*8*pi^(2)/3*reduced_hessian;

    res_num = atom_num/2;
    order_parameter = zeros(res_num,1);
    displacement_product = zeros(6,6);

    for i = 0:(res_num-1)
        covariance_res = reduced_covariance(6*i+1:6*i+6,6*i+1:6*i+6);
        mean_struct_vector = [PDB_Structure(2*i+1).coord',PDB_Structure(2*i+2).coord'];
        displacement_product(:,:) = mean_struct_vector'*mean_struct_vector + covariance_res;

        xijxij = displacement_product(1,1) + displacement_product(4,4) - 2*displacement_product(1,4);
        yijyij = displacement_product(2,2) + displacement_product(5,5) - 2*displacement_product(2,5);
        zijzij = displacement_product(3,3) + displacement_product(6,6) - 2*displacement_product(3,6);
        xijyij = displacement_product(1,2) + displacement_product(4,5) - displacement_product(1,5) - displacement_product(4,2);
        xijzij = displacement_product(1,3) + displacement_product(4,6) - displacement_product(1,6) - displacement_product(4,3);
        yijzij = displacement_product(2,3) + displacement_product(5,6) - displacement_product(2,6) - displacement_product(5,3);

        correction_factor = (xijxij + yijyij + zijzij);

        xijxij = xijxij/correction_factor;
        yijyij = yijyij/correction_factor;
        zijzij = zijzij/correction_factor;
        xijyij = xijyij/correction_factor;
        xijzij = xijzij/correction_factor;
        yijzij = yijzij/correction_factor;

        order_parameter(i+1) = (3/2)*(xijxij^(2) + yijyij^(2) + zijzij^(2) + 2*xijyij^(2) + 2*xijzij^(2) + 2*yijzij^(2)) - (1/2);
    end
end

