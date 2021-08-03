ME1 = MException('X', 'Z');
condition = true;

try
    if condition
        throw(ME1);
    end
catch
    
end
