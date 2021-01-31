function result = endsWith(str,pat)
%%atEndsWith    Check if string ends with pattern
%%
%%RESULT=ATENDSWITH(STR,PAT)
%%
%%   STR:     Target string
%%   PAT:     Pattern to look for
  result = 0;

  strLength = length(str);
  patLength = length(pat);

  if strLength < patLength
    return;
  endif

  str = tolower(str);
  pat = tolower(pat);

  result = strcmp(substr(str, -patLength), pat);
end
