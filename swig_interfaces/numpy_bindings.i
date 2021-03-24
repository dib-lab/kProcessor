
%{
    #define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"


%init %{
    import_array();
%}



%apply (unsigned long long* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(unsigned long long* rangevec, int rows, int cols)}
