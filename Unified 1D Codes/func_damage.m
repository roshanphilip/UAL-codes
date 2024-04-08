function [Damage,dDdestar] = func_damage(DamageThreshholdStrain,a,b,Se)

 % This function calculates the Damage induced in each element 
                    
    if (Se>=DamageThreshholdStrain)
        
        Damage=1-(((DamageThreshholdStrain*(1-a))/Se)+(a/(exp(b*(Se-DamageThreshholdStrain)))));
        
        % Partial derivative of damage w.r.t. equivalent strain     
        dDdestar=(a*b*(exp(b*(DamageThreshholdStrain-Se))))+(((1-a)*DamageThreshholdStrain)/(Se^2));
               
    else
        
        Damage=0;
        dDdestar=0;

    end

end

