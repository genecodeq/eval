/*
  Copyright (c) 2015, Fonleap Ltd
  All rights reserved.

  Redistribution and use in source and binary forms, with or without modification,
  are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice, this
     list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice, this
     list of conditions and the following disclaimer in the documentation and/or
     other materials provided with the distribution.

  3. Neither the name of the copyright holder nor the names of its contributors may
     be used to endorse or promote products derived from this software without
     specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
  NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
  OF SUCH DAMAGE.


 *
 * mergeq.cpp - Processes fastq or sam file and merges quality scores (separated by newline)
 * from stdin. The new fastq or sam file will be output in stdout.
 */

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <string>

enum SeqFileType { SAM, FASTQ };

using namespace std;

bool endsWith(const string& base, const string& pattern) {
    if (base.length() < pattern.length())
        return false;
    return (base.substr(base.length() - pattern.length()).compare(pattern) == 0);
}


int main(int argc, char *argv[]) {
    if (argc <2) {
        fprintf(stderr, "Usage: %s [filename] (quality scores per line from stdin) (merged result in stdout)\n", argv[0]);
        return 0;
    }
    std::string inputfilepath = argv[1];
    SeqFileType inputfiletype = endsWith(inputfilepath, ".fastq") ? FASTQ : SAM;

    FILE *in = fopen(argv[1], "r");
    if (in == NULL) {
        fprintf(stderr, "Unable to open input file: %s [%s]\n", argv[1], strerror(errno));
        return -1;
    }
    char line[65536], line2[65536];
    while (fgets(line, sizeof(line), in)) {
        if (inputfiletype == FASTQ) {
            // write the first line to the output
            fputs(line, stdout);
            // write the other two lines to the output
            for (int i = 0; i < 2; ++i) {
                if (!fgets(line, sizeof(line), in)) {
                    printf("Failed to read fastq entry: %s\n", line);
                    return -1;
                }
                fputs(line, stdout);
            }
            // retrieve the qscore line
            if (!fgets(line, sizeof(line), in)) {
                printf("Failed to read fastq entry: %s\n", line);
                return -1;
            }
            if (!fgets(line, sizeof(line), stdin)) {
                printf("Failed to read stdin fastq quality score entry\n");
            }
            // write out
            fputs(line, stdout);
            continue;
        }

        if (line[0] == '@') {
            fputs(line, stdout);
            continue;
        }
        bool wasWhitespace = true;
        unsigned int colNum=0, colStartPos=0, pos=0;
        while (line[pos]!=0) {
            bool whitespace = (line[pos] == ' ' or line[pos] == '\t');
            if (whitespace && !wasWhitespace) {
                wasWhitespace = true;
                if (colNum == 11) {
                    if (!fgets(line2, sizeof(line2), stdin)) {
                        printf("Failed to read stdin fastq quality score entry\n");
                        return -1;
                    }
                    for (unsigned int idx=colStartPos; idx<pos; ++idx)
                        line[idx] = line2[idx-colStartPos];
                }
            }
            else if (!whitespace && wasWhitespace) {
                colNum++;
                colStartPos = pos;
                wasWhitespace = false;
            }
            pos++;
        }
        fputs(line, stdout);
    }
    return 0;
}
