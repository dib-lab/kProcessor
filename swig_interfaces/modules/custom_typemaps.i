
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

%typemap(out) std::vector<uint32_t> {
    int ln = $1.size();
    $result = PyList_New(0);
    std::vector<uint32_t>::iterator iter;
    for (int i = 0; i < ln ; i++) {
        PyList_Append($result, PyInt_FromLong($1[i]));
    }
}

/*----------------------------------------------*/