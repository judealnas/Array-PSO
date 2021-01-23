classdef (Abstract) Problem < handle
   properties (Abstract)
       dimensions
   end
   
   methods (Abstract)
       evaluate(obj, input)
   end
end