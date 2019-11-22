%module kProcessor

%{
/* Everything in this block will be copied in the wrapper file. We include the C header file necessary to compile the interface */
#include "kDataFrame.hpp" // including kDataframe
%}

using namespace std; // Extremly important

%include "std_vector.i"     /*Using std::vector typemaps*/
%include std_map.i

%template(IntVector) vector<int>; /*vector to tuple conversion*/

%import stdint.i            /*This mainly used for converting python int to C++ uint64_t*/
%include std_string.i       /*And this for converting python str to C++ std::string*/

/*  ~~~~ENABLED FOR NOW~~~~~~  */
/*  ~~~~COMMENT TO DISABLE namesMap interface~~~~  */
/*  ~~~~WORKS ONLY IN SWIG4.0 ~~~~  */
%include std_unordered_map.i
/*
%template(MAPsi) unordered_map<int, std::string>;
%template(MAPis) unordered_map<std::string, int>;
*/
%template(kmerDecoderParams) map<string, int>;
%include "swig_interfaces/modules/custom_typemaps.i"

/*Just copy/paste the snippet I'm interested in to be wrapped!*/

// no need to include subclasses if all their methods (that we are interested in) are defined in the superclass

/******** kDataFrame Interface ************/
%include "swig_interfaces/kDataFrame/kDataFrameIterator.i"
%include "swig_interfaces/kDataFrame/kmerRow.i"
%include "swig_interfaces/kDataFrame/kDataFrame.i"
%include "swig_interfaces/kDataFrame/kDataFrameMQF.i"
%include "swig_interfaces/kDataFrame/kDataFrameMAP.i"
%include "swig_interfaces/kDataFrame/kDataFramePHMAP.i"
%include "swig_interfaces/kDataFrame/kDataFrameBMQF.i"
/******** kDataFrame Interface ************/


/******** colored_kDataFrame Interface ************/

%{
#include "colored_kDataFrame.hpp"
%}

%include "swig_interfaces/colored_kDataFrame.i"

%template(colorsList) vector<uint32_t>;

/******** colored_kDataFrame Interface ************/

/******** colorTable Interface ************/

%{
#include "colorTable.hpp"
%}
%include "swig_interfaces/colorTable.i"

/******** colorTable Interface ************/

%{
#include "algorithms.hpp" // including algorithms
%}

%template(kFramesVector) vector<kDataFrame*>; /*vector to tuple conversion*/

%include "swig_interfaces/algorithms/algorithms.i"


/******** HashUtils Interface ************/

%{
#include "HashUtils/hashutil.hpp" // including HashUtils
%}

%include "swig_interfaces/HashUtils/hashutil.i"

%{
#include "Utils/kmer.h" // including Kmer
%}

%include "swig_interfaces/Utils/kmer.i"

/******** END kDataFrame ************/