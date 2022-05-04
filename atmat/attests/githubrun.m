v=runtests('attests',ProcedureName="x_*");
success=all(v.Passed);
exit(~success);