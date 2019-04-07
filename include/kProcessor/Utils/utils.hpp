#ifndef UTILS_HPP
#define UTILS_HPP
#include <stdint.h>
#include <string>
#include <vector>


namespace kProcessor::utils{
bool has_suffix(const std::string& s, const std::string& suffix);
// Taken from
// https://stackoverflow.com/questions/19189014/how-do-i-find-files-with-a-specific-extension-in-a-directory-that-is-provided-by
std::vector<std::string> GetFilesExt(const char *dir, const char *ext);
std::string last_part(std::string str, char c);
std::string first_part(std::string str, char c);
bool FileExists(std::string);
}

#endif
