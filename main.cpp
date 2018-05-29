/*  bamtk.c -- main samtools command front-end.

    Copyright (C) 2008-2017 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */





#include <iostream>
#include <string>

using namespace std;

#define KPROCESSOR_VERSION "0.1"

//int KmerCounter_main(int argc, char *argv[]);

string KprocessorVersion()
{
    return KPROCESSOR_VERSION;
}

static void usage()
{
    /* Please improve the grouping */

    cout<<
"\n"
"This is Kprocessor version "<<KprocessorVersion()<<" developed by Mostafa Shokrof <mostafa.shokrof@gmail.com>\n"
<<"Usage:   kprocessor <command> [options]\n"<<
"\n"<<
"Commands:\n"<<
"  -- Counter\n"<<
"     count          count kmers in seqeunces file \n";
}

int main(int argc, char *argv[])
{
    if (argc < 2) { usage(); return 1; }

    if (string(argv[1])== "help"  || string(argv[1])== "--help" ) {
        if (argc == 2) { usage(); return 0; }

        // display help of the specific tool
        argv++;
        argc = 2;
    }

    int ret = 0;
    if (string(argv[1])== "view")           ret = 1;
    else if (string(argv[1])== "--version" ) {
      cout<<"This is Kprocessor version "<<KprocessorVersion()<<" developed by Mostafa Shokrof <mostafa.shokrof@gmail.com>\n";
    }
    else {
      cout<<"[main] unrecognized command "<< argv[1]<<endl;
        return 1;
    }
    return ret;
}
