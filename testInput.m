function testInput()
    p = inputParser;
    defaultDamping = false;
    addOptional(p, 'Damped', defaultDamping, @islogical);
    parse(p, 'Damped', true);
    disp('Damped is:');
    disp(p.Results.Damped);
end
