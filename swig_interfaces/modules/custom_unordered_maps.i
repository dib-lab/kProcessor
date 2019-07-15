
/* unordered_map<int, std::string> */
%typemap(out) std::unordered_map<int, std::string> {
    $result = PyDict_New();
    std::unordered_map<int, string>::iterator iter;

    for (iter = $1.begin(); iter != $1.end(); ++iter) {
        PyDict_SetItem($result, PyInt_FromLong(iter->first), PyString_FromString(iter->second.c_str()));
    }
}

/* unordered_map<std::string, int> */
%typemap(out) std::unordered_map<std::string, int> {
    $result = PyDict_New();
    std::unordered_map<string, int>::iterator iter;
    for (iter = $1.begin(); iter != $1.end(); ++iter) {
        PyDict_SetItem($result, PyString_FromString(iter->first.c_str()) ,PyInt_FromLong(iter->second));
    }
}

/*----------------------------------------------*/