%module kProcessor

%{
/* Everything in this block will be copied in the wrapper file. We include the C header file necessary to compile the interface */
#include "kDataFrame.hpp" // including kDataframe
#include "defaultColumn.hpp"
%}

using namespace std; // Extremly important

%include "std_vector.i"     /*Using std::vector typemaps*/
%include std_map.i

%template(IntVector) vector<int>; /*vector to tuple conversion*/
%template(FloatVector) vector<float>; /*vector to tuple conversion*/



//typedef long int 		    int64_t;
typedef unsigned long int 	uint64_t;
//%apply long int  { int64_t };
%apply unsigned long int { uint64_t };

%include "stdint.i"            /*This mainly used for converting python int to C++ uint64_t*/
%include std_string.i       /*And this for converting python str to C++ std::string*/

/*  ~~~~ENABLED FOR NOW~~~~~~  */
/*  ~~~~COMMENT TO DISABLE namesMap interface~~~~  */
/*  ~~~~WORKS ONLY IN SWIG4.0 ~~~~  */
%include std_unordered_map.i

%template(MAPsi) unordered_map<int, std::string>;
%template(MAPis) unordered_map<std::string, int>;

%include "swig_interfaces/modules/custom_typemaps.i"

%template(colorsList) vector<uint32_t>;
%template(kmerDecoderParams) map<string, int>;
%template(DoubleVector) vector<double>; /*vector to tuple conversion*/
// %template(map_str_dbl_bec) std::unordered_map<std::string, std::vector<double>>;

/*Just copy/paste the snippet I'm interested in to be wrapped!*/

// no need to include subclasses if all their methods (that we are interested in) are defined in the superclass

%{
#include "algorithms.hpp" // including algorithms
#include "kDataframes/kDataFrameBlight.hpp"
#include "kDataframes/kDataFrameMQF.hpp"
#include "kDataframes/kDataFrameSTL.hpp"
#include "kDataframes/kDataFrameBMQF.hpp"
#include "defaultColumn.hpp"
%}

%include "swig_interfaces/kDataFrame/kmerRow.i"
%include "swig_interfaces/algorithms/algorithms.i"

/******** kDataFrame Interface ************/
%include "swig_interfaces/kDataFrame/kDataFrameIterator.i"
%include "swig_interfaces/kDataFrame/kDataFrame.i"
%include "swig_interfaces/kDataFrame/kDataFrameFactory.i"
%include "swig_interfaces/kDataFrame/kDataFrameMQF.i"
%include "swig_interfaces/kDataFrame/kDataFrameSTL.i"
%include "swig_interfaces/kDataFrame/kDataFrameBlight.i"
%include "swig_interfaces/kDataFrame/kDataFrameBMQF.i"
%include "swig_interfaces/kDataFrame/dbgIterator.i"
%include "swig_interfaces/kDataFrame/defaultColumn.i"

/******** kDataFrame Interface ************/

%template(kFramesVector) vector<kDataFrame*>; /*vector to tuple conversion*/


// extend_algorithms

// %{
// #include "extend_algorithms.hpp" // including algorithms
// %}

// %include "swig_interfaces/algorithms/extend_algorithms.i"

/******** HashUtils Interface ************/

/******** END kDataFrame ************/