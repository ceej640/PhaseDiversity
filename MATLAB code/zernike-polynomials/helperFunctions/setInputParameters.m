function obj = setInputParameters(className, obj, par)
    nSetParameters = 0;
    propertyNames = properties(className);
    for k = 1:size(propertyNames,1)
        if isfield(par, propertyNames{k})
            obj.(propertyNames{k}) = par.(propertyNames{k});
            nSetParameters = nSetParameters + 1;
        end
    end
    if nSetParameters ~= length(fieldnames(par))
        warning('Parameter structure contains invalid fields! Not all parameters were set!')
        parFieldnames = fieldnames(par);
        for k = 1:length(parFieldnames)
            if ~any(strcmp(propertyNames,parFieldnames{k}))
                warning(strcat(parFieldnames{k},' is not a valid property!'))
            end
        end
    end
end