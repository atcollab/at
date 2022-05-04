v=runtests('attests',ProcedureName=["x_*","linopt");
success=all(cat(1,v.Passed));
exit(~success);