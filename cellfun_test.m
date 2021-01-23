func1 = @(x) x.^2;
func2 = @(y) y.^3 - y;
func3 = @(z) 4;

fm = {func1 func2 func3}

x = [2 0 1];
results = cellfun(@(f) f(x),fm, 'UniformOutput', false)