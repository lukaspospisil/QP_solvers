classdef SmalseOptions
    
    properties
        rho % rho > 0
        eta % eta > 0
        beta % 0 < beta < 1
        M % M0 > 0

        myeps
        maxit
        
        normA
        normBB
        
        dispdebug       % plot info about progress or be quite (true/false)
    end
    
    methods
        function obj = SmalseOptions()
            % constructor - set default values
            obj.rho = 1e3;
            obj.eta = 1e-4;
            obj.beta = 0.5;
            obj.M = 1;
            
            obj.myeps = 1e-6;
            obj.maxit = 1e2;
            
            obj.dispdebug = true;
            
        end
        
    end
end    
