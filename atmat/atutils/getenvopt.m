function value = getenvopt(name,default)
%   GETENVOPT(NAME, DEFAULTVALUE)
%   Looks for an environment variable and return a default value if absent

value=getenv(name);
if isempty(value)
    value=default;
end
end

