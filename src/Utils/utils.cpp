#include <stdint.h>
#include <string>
#include <gqf.h>
#include <seqan/seq_io.h>
#include <sys/stat.h>
#include <vector>


namespace kProcessor::utils{
// void removeReadsWithN(std::string inputFilename,std::string outputFilename)
// {
//   SeqFileIn seqFileIn(inputFilename.c_str());
//   SeqFileOut seqFileOut(outputFilename.c_str());
//   CharString id;
//   CharString quals;
//   std::string readT;
//
//   while(!atEnd(seqFileIn))
//   {
//     bool hasN=false;
//     readRecord(id, readT,quals, seqFileIn);
//     for(uint64_t i=0;i<readT.size();i++)
//     {
//       if(readT[i]=='N')
//       {
//         hasN=true;
//         break;
//       }
//     }
//     if(!hasN)
//     {
//       writeRecord(seqFileOut,id,readT,quals);
//     }
//   }
//
//
// }
bool has_suffix(const std::string& s, const std::string& suffix)
		{
			return (s.size() >= suffix.size()) && equal(suffix.rbegin(),
																									suffix.rend(), s.rbegin());
		}
// Taken from
// https://stackoverflow.com/questions/19189014/how-do-i-find-files-with-a-specific-extension-in-a-directory-that-is-provided-by
std::vector<std::string> GetFilesExt(const char *dir, const char *ext)
  {
    DIR *folder = opendir(dir);
		if (!folder) {
			std::cerr << "Directory doesn't exist " << dir << std::endl;
			exit(1);
		}

		std::vector<std::string> ret;
		dirent *entry;
		while((entry = readdir(folder)) != NULL)
		{
			if(has_suffix(entry->d_name, ext))
			{
				std::string filename(entry->d_name);
				std::string dirname(dir);
				ret.push_back(std::string(dirname + filename));
			}
		}

		return ret;
}

std::string last_part(std::string str, char c) {
	uint64_t found = str.find_last_of(c);
	return str.substr(found + 1);
}

std::string first_part(std::string str, char c) {
	uint64_t found = str.find_first_of(c);
	return str.substr(0, found);
}

// copied from mantis
// Taken from
// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
bool FileExists(std::string Spath) {
	const char* path=Spath.c_str();
	struct stat fileStat;
	if (stat(path, &fileStat)) {
		return false;
	}
	if (!S_ISREG(fileStat.st_mode)) {
		return false;
	}
	return true;
}

}
