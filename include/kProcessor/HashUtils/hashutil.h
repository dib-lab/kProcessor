/*
* =====================================================================================
*
*       Filename:  hashutil.h
*
*    Description:
*
*        Version:  1.0
*        Created:  04/18/2016 04:49:32 PM
*       Revision:  none
*       Compiler:  gcc
*
*         Author:  Prashant Pandey (ppandey@cs.stonybrook.edu)
*                  Rob Patro (rob.patro@cs.stonybrook.edu)
*                  Rob Johnson (rob@cs.stonybrook.edu)
*   Organization:  Stony Brook University
*
* =====================================================================================
*/

#ifndef _HASHUTIL_H_
#define _HASHUTIL_H_

#include <sys/types.h>
#include <string>
#include <stdlib.h>
#include <stdint.h>

using namespace std;
class Hasher{
public:
	virtual uint64_t hash(string key){return 0;};
	virtual uint64_t hash(uint64_t key){return 0;};
	virtual Hasher* clone(){return this;};
};
class MumurHasher: public Hasher{
	using Hasher::hash;
private:
	uint64_t seed;
public:
	MumurHasher(uint64_t Iseed){seed=Iseed;}
	Hasher* clone() { return new MumurHasher(seed);}
	uint64_t hash(string kmer);
};

class IntegerHasher: public Hasher{
private:
	uint64_t mask;
	uint64_t kSize;
public:
	IntegerHasher(uint64_t kSize);
	Hasher* clone() { return new IntegerHasher(kSize);}
	uint64_t hash(string key);
	uint64_t hash(uint64_t key);
	string Ihash(uint64_t key);
};



#endif  // #ifndef _HASHUTIL_H_
