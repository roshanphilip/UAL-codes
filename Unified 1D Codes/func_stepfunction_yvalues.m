function [b] = func_stepfunction_yvalues(a)
% This function calcualtes the y values of the step function used in
% plotting damage/strain against the gauss points

% a matrix 
    size_a=size(a,1);
    [max_a,pos_max_a]=max(a);

    a_flagged_values=[a(1,1);a(pos_max_a,1);a(size_a,1)];

    a_flagged_positions=[1;pos_max_a;size_a];

    j=1;
    for i=1:1:size_a
        if ismember(i,a_flagged_positions)~=1       
            a_nonflagged_positions(j,1)=i;
            j=j+1;
        end
    end
    clear j
    a_nonflagged_positions;

    for i=1:1:size(a_nonflagged_positions,1)
        a_nonflagged_values(i,1)=a(a_nonflagged_positions(i,1),1);
    end
    a_nonflagged_values;

% k matrix

    size_k=(size_a*2)-1;
    k=zeros(size_k,1);
    mod_pos_max=(pos_max_a*2)-1;

    k_1=1;
    k_max_damage=mod_pos_max;
    k_max_adjacent_plus=mod_pos_max+1;
    k_max_adjacent_minus=mod_pos_max-1;
    k_end=(size_a*2)-1;

    k_flagged_positions=[k_1;k_max_adjacent_minus;k_max_damage;k_max_adjacent_plus;k_end];

    j=1;
    for i=1:1:size_k
        if ismember(i,k_flagged_positions)~=1
            k_nonflagged_positions(j,1)=i;
            j=j+1;
        end
    end
    clear j
    k_nonflagged_positions;

    for i=1:1:size(k_flagged_positions,1)
        if i==1 
        k_flagged_values(i,1)=a_flagged_values(i,1);
        elseif i==size(k_flagged_positions,1)
        k_flagged_values(i,1)=a_flagged_values(end,1);    
        else 
        k_flagged_values(i,1)=a_flagged_values(2,1);   
        end
    end
    
    k_flagged_values;
    
    j=1;
    for i=1:1:(size(a_nonflagged_values,1))
        k_nonflagged_values(j,1)=a_nonflagged_values(i,1);
        k_nonflagged_values(j+1,1)=a_nonflagged_values(i,1);
        j=j+2;
    end
    clear j
    k_nonflagged_values;  
        
    % Calculating values of b matrix
    
    for i=1:1:size_k
        if ismember(i,k_flagged_positions)==1
            o=find(k_flagged_positions==i);
            b(i,1)=k_flagged_values(o,1);        
        elseif ismember(i,k_nonflagged_positions)==1
            p=find(k_nonflagged_positions==i);
            b(i,1)=k_nonflagged_values(p,1);
        end      
    end
    
    b;

end

