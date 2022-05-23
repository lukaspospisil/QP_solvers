classdef MprgpOptions
    
    properties
        gamma

        maxit
        myeps
        
        smalsem
        smalsem_obj
        
        dispdebug       % plot info about progress or be quite (true/false)
    end
    
    methods
        function obj = MprgpOptions()
            % constructor - set default values
            obj.gamma = 0.75;
            
            obj.maxit = 1e3;
            obj.myeps = 1e-6;
            
            obj.smalsem = false;
            obj.smalsem_obj = struct;
            
            obj.dispdebug = true;
        end
        
    end
end    
