v=runtests('attests',ProcedureName="x_*");
success=all(cat(1,v.Passed));
exit(~success);