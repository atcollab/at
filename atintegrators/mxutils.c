#include <string.h>

static char* ElemName(const mxArray *ElemData)
{
    static char buffer[50];
    mxArray *tmpmxptr=mxGetField(ElemData,0,"FamName");
    if (tmpmxptr==NULL)
        strcpy(buffer,"Invalid");
    else
        mxGetString(tmpmxptr,buffer,50);
    return buffer;
}

static void RequiredFieldError(const mxArray *ElemData, const char *fieldname)
{
    mexErrMsgIdAndTxt("AT:MissingField",
    "Element '%s': Required field '%s' was not found in the element data structure",
    ElemName(ElemData), fieldname);
}

static int GetRequiredFieldNumber(const mxArray *ElemData, const char *fieldname)
{
    int fnum=mxGetFieldNumber(ElemData, fieldname);
    if (fnum<0) RequiredFieldError(ElemData, fieldname);
    return fnum;
}

static mxArray *GetRequiredField(const mxArray *ElemData, const char *fieldname)
{
    mxArray *tmpmxptr=mxGetField(ElemData,0,fieldname);
    if (tmpmxptr==NULL) RequiredFieldError(ElemData, fieldname);
    return tmpmxptr;
}

