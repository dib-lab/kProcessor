#include "restriction_tags.hpp"

intFuncMap tagToFunc;

tag::tag(map <string, int> tags) {
    for (auto & tag_info: tags) {
        string tag_name = tag_info.first;
        int tag_value = tag_info.second;
        add_restriction(tag_name, tag_value);
    }

    tagToFunc["min_kSize"] = tag_min_kSize;
    tagToFunc["max_kSize"] = tag_max_kSize;
}

void tag::add_restriction(string tag_name, int value) {
    bool valid_restriction = (restrictions.find(tag_name) != restrictions.end());
    if (valid_restriction) {
        restrictions[tag_name] = value;
        active_tags.push_back(tag_name);
    } else throw std::invalid_argument("(" + tag_name + ") is not recognized.");
}

// Checks for rules violation in active_tags
void tag::check_restrictions() {
    // Check for missing tags
    for (auto & active_tag: active_tags) {
        int tag_value = restrictions[active_tag];
        std::string
        function = std::string(active_tag);

        // Call the validation function
        intFuncMap::iterator tagFunc = tagToFunc.find(function);
        if (tagFunc != tagToFunc.end()) {
            tag m;
            (m.*(tagFunc -> second))(tag_value);
        }
    }
}

// Check functions

void tag::tag_min_kSize(int value) {
    if (value < restrictions["min_kSize"]) {
        throw std::logic_error("kSize must: " + to_string(restrictions["min_kSize"]) + " < kSize > " + to_string(restrictions["max_kSize"]));
    }
}

void tag::tag_max_kSize(int value) {
    if (value > restrictions["max_kSize"]) {
        throw std::logic_error("kSize must: " + to_string(restrictions["min_kSize"]) + " < kSize > " + to_string(restrictions["max_kSize"]));
    }
}