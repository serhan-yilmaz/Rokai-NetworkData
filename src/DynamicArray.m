classdef DynamicArray < handle
    
    properties (SetAccess = private)
        array;
        capacity;
        currentSize = 0;
        isarraycell = false;
    end
    
    methods
        function obj = DynamicArray(capacity, type)
            if(nargin >= 1)
                obj.capacity = capacity;
            else
                obj.capacity = 100;
            end
            if(nargin < 2)
               type = 'numeric'; 
            end
            switch(type)
                case 'numeric'
                    obj.array = zeros(obj.capacity, 1);
                case 'cell'
                    obj.array = cell(obj.capacity, 1);
                    obj.isarraycell = true;
                otherwise
                    error('Invalid type argument.');
            end
        end
        
        function [obj] = enlarge(obj)
            obj.array = [obj.array; zeros(obj.capacity, 1)];
            obj.capacity = obj.capacity * 2;
        end
        
        function [obj] = add(obj, list)
            obj = insert(obj, list);
        end
        
        function [obj] = insert(obj, list)
            if(~isvector(list))
               error('List must be a vector.'); 
            end
%             list = squeeze(list);
            if(obj.isarraycell && ~iscell(list))
               list = {list}; 
            end
            nElement = length(list);
            while(obj.capacity < obj.currentSize + nElement)
               obj = enlarge(obj);
            end
            indexStart = 1 + obj.currentSize;
            indexEnd = obj.currentSize + nElement;
            update_indices = indexStart:indexEnd;
            obj.array(update_indices) = list;
            obj.currentSize = obj.currentSize + nElement;
        end
        
        function [val] = get(obj, index)
            if(index > obj.currentSize)
               error('Array index is out of bounds.');
            end
            if(obj.isarraycell)
               val = obj.array{index};
            else
               val = obj.array(index);
            end
        end
        
        function [obj] = set(obj, index, value)
            if(index > obj.currentSize)
               error('Array index is out of bounds.');
            end
            if(obj.isarraycell)
               obj.array{index} = value;
            else
               obj.array(index) = value;
            end
        end
        
        function [currentSize] = size(obj)
            currentSize = obj.currentSize;
        end
        
        function [array] = finalize(obj)
            array = obj.array(1:obj.currentSize);
        end
    end
    
end

