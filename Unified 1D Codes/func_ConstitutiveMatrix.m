function [ConstitutiveMatrix,YoungsModulus] = func_ConstitutiveMatrix(n_TotalElements,E,Area,l_EachElement,damaged_length,l_Total,notch_type,notch)

% // This function calculates the constitutive matrix of the 1D bar based on the material properties defined in the input file //

ConstitutiveMatrix=zeros(1,n_TotalElements);

middle_element=round(n_TotalElements/2);

% Add a notch to the middle of the bar 

% 1. No notch

if notch_type==1
    for i=1:1:n_TotalElements
        ConstitutiveMatrix(1,i)=(E*Area);
        YoungsModulus(i)=E;
    end

% 2. Notch in single element

elseif notch_type==2
    for i=1:1:n_TotalElements
        if i==middle_element
            ConstitutiveMatrix(1,i)=(notch*E*Area);
            YoungsModulus(i)=notch*E;
        else
            ConstitutiveMatrix(1,i)=(E*Area);
            YoungsModulus(i)=E;
        end

    end

% 3. Notch across a fixed length in the middle of the bar

elseif notch_type==3
    middle_element=round(n_TotalElements/2);

    % Check if there are enough elements to apply damage properly across the length of the bar
    if damaged_length<l_EachElement
        disp("Warning: Too few elements! Use atleast "+round(l_Total/damaged_length)+" elements or don't apply a notch")
        pause
    end

    n_damaged_elements=round(damaged_length/l_EachElement);

    % Identify the elements where the notch should be applied

    zero_plus_elements=round(n_damaged_elements/2-0.1);
    if rem(n_damaged_elements,2)~=0
        zero_minus_elements=zero_plus_elements;
    elseif rem(n_damaged_elements,2)==0
        zero_minus_elements=zero_plus_elements-1;
    end
   
    add=(round(-zero_minus_elements:zero_plus_elements));
    add=add';
    
    c_notch=zeros(n_damaged_elements,1);    

    for i=1:1:n_damaged_elements
        c_notch(i,1)=middle_element+add(i,1);
    end
   
    counter=0;

    for i=1:1:n_TotalElements
        [check1,~]=ismember(i,c_notch);
        if check1==1
            ConstitutiveMatrix(1,i)=(notch*E*Area);
            YoungsModulus(i)=notch*E;
            counter=counter+1;
        else
            ConstitutiveMatrix(1,i)=(E*Area);
            YoungsModulus(i)=E;
        end
    end

end

end

