
// /* unordered_map<int, std::string> */
// %typemap(out) std::unordered_map<int, std::string> {
//     $result = PyDict_New();
//     std::unordered_map<int, string>::iterator iter;

//     for (iter = $1.begin(); iter != $1.end(); ++iter) {
//         PyDict_SetItem($result, PyInt_FromLong(iter->first), PyString_FromString(iter->second.c_str()));
//     }
// }

// /* unordered_map<std::string, int> */
// %typemap(out) std::unordered_map<std::string, int> {
//     $result = PyDict_New();
//     std::unordered_map<string, int>::iterator iter;
//     for (iter = $1.begin(); iter != $1.end(); ++iter) {
//         PyDict_SetItem($result, PyString_FromString(iter->first.c_str()) ,PyInt_FromLong(iter->second));
//     }
// }

// %typemap(out) std::vector<uint32_t> {
//     int ln = $1.size();
//     $result = PyList_New(0);
//     std::vector<uint32_t>::iterator iter;
//     for (int i = 0; i < ln ; i++) {
//         PyList_Append($result, PyInt_FromLong($1[i]));
//     }
// }

// /*----------------------------------------------*/


// //%typemap(out) std::unordered_map<std::string, std::vector<std::vector<uint32_t>>> {
// //
// //    $result = PyDict_New();
// //    std::unordered_map<string, std::vector<std::vector<uint32_t>>> iter;
// //    for (iter = $1.begin(); iter != $1.end(); ++iter) {
// //        std::cout << "Hello Fuckin World!\n";
// //        PyDict_SetItem($result, PyString_FromString(iter->first.c_str()) , iter->second);
// //    }
// //
// //}

// // Typemap to convert kmerSources to Python dictionary
// // Thanks to https://stackoverflow.com/a/27507865
// %typemap(out) std::unordered_map<std::string, std::vector<std::vector<uint32_t>>> &{
//     $result = PyDict_New();
//     std::unordered_map<string, std::vector<std::vector<uint32_t>>>::iterator map_iter;
//     for (map_iter = $1.begin(); map_iter != $1.end(); ++map_iter) {
//         for(int i = 0; i < $1->second.size(); ++i){
//             int subLength = $1->second->data()[i].size();
//             npy_intp dims[] = { subLength };
//             PyObject* temp = PyArray_SimpleNewFromData(1, dims, NPY_INT, $1->second->data()[i].data());
//             $result = PyDict_SetItem($result, PyString_FromString(iter->first.c_str()), temp);
//         }
//     }
// }


// // Typemap to convert kmerSources to Python dictionary
// %typemap(out) std::unordered_map<std::string, std::vector<uint32_t>> &{
//     $result = PyDict_New();

//     std::unordered_map<string, std::vector<std::vector<uint32_t>>>::iterator map_iter;

//     for (map_iter = $1.begin(); map_iter != $1.end(); ++map_iter) {
//             int len = map_iter->second.size();

//             $vec = PyList_New(len);

//             for(int i = 0; i < len; i++){
//                 PyList_Append($vec, PyInt_FromLong(map_iter->second[i]));
//             }

//             $result = PyDict_SetItem($result, PyString_FromString(map_iter->first.c_str()), $vec);

//     }

// }

%include <std_vector.i>
%template() std::vector<double>;


%typemap(out) std::unordered_map<std::string, std::vector<double>> &{

    $result = PyDict_New();

    std::unordered_map<string, std::vector<double>>::iterator map_iter;

    for (map_iter = $1.begin(); map_iter != $1.end(); ++map_iter) {

            int len = map_iter->second.size();

            $vec = PyList_New(len);

            for(int i = 0; i < len; i++){
                PyList_Append($vec, PyFloat_FromDouble(map_iter->second[i]));
            }

            $result = PyDict_SetItem($result, PyString_FromString(map_iter->first.c_str()), $vec);
    }
}