function [GlobalCoordinates_CurrentElement] = func_GlobalCoordinates_CurrentElement(n_NodesPerElement,el,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode)

% // This function calculates the Global Coordinates of the current element under consideration //

    GlobalCoordinates_CurrentElement=zeros(1,n_NodesPerElement);

    for i=1:n_NodesPerElement    
        GlobalCoordinates_CurrentElement(i)=Globalcoordinates_EachNode(connect_LocaltoGlobal_Nodenumbers(i,el));    
    end

end

