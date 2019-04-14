%module kProcessor

%{
/* Everything in this block will be copied in the wrapper file. We include the C header file necessary to compile the interface */
#include "kDataFrame.hpp" // including kDataframe
%}

using namespace std; // Extremly important

%include "std_vector.i"     /*Using std::vector typemaps*/
%template(IntVector) vector<int>; /*vector to tuple conversion*/

%import stdint.i            /*This mainly used for converting python int to C++ uint64_t*/
%include std_string.i       /*And this for converting python str to C++ std::string*/

/*Just copy/paste the snippet I'm interested in to be wrapped!*/

// no need to include subclasses if all their methods (that we are interested in) are defined in the superclass

/******** kDataFrame Interface ************/
%include "swig_interfaces/kDataFrame/kDataFrameIterator.i"
%include "swig_interfaces/kDataFrame/kmerRow.i"
%include "swig_interfaces/kDataFrame/kDataFrame.i"
%include "swig_interfaces/kDataFrame/kDataFrameMQF.i"
%include "swig_interfaces/kDataFrame/kDataFrameMAP.i"

%{
#include "algorithms.hpp" // including algorithms
%}

%template(kFramesVector) vector<kDataFrame*>; /*vector to tuple conversion*/

%include "swig_interfaces/algorithms/algorithms.i"
/******** END kDataFrame ************/