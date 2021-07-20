#ifndef _RESTRICTION_TAGS_H_
#define _RESTRICTION_TAGS_H_

#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>

using std::map;
using std::string;
using std::vector;
using std::to_string;
using std::stoi;
using std::cout;
using std::cerr;
using std::endl;


class tag {

    public:
        vector<string> active_tags;
        map<string, int> restrictions = {
            {"min_kSize", 7},
            {"max_kSize", 31},
            {"sorted", false},
        };
        
        tag(){}
        tag(map<string, int> tags);

        void add_restriction(string tag_name, int value);
        void check_restrictions();

        void tag_min_kSize(int value);
        void tag_max_kSize(int value);

        ~tag() {}
    };


typedef void (tag::*intFunc)(int);
typedef map<string, intFunc> intFuncMap;


#endif
