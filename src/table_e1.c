#include <R.h>
#include <Rdefines.h>

SEXP table_e1(SEXP a, SEXP b)
{
    PROTECT(a = AS_INTEGER(a));
    PROTECT(b = AS_INTEGER(b));

    int i, index = 0, n = GET_LENGTH(a);
    if (n <= 0) {
        UNPROTECT(2);
        return R_NilValue;
    }

    SEXP x;
    PROTECT(x = NEW_INTEGER(6));

    const int *pa = INTEGER_POINTER(a);
    const int *pb = INTEGER_POINTER(b);

    int *px = INTEGER_POINTER(x);
    for(i = 0; i <= 5; i++){
      *(px + i) = 0;
    }

    for (i = 0; i < n; ++i){
      index = *pb * 3 + *pa;
      pa++;
      pb++;
      *(px + index) = *(px + index) + 1;
    }

    UNPROTECT(3);

    return x;
}

