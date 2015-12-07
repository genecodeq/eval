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
 * qsxtract.cpp - A tool for extracting quality scores from SAM or FASTQ files.
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

void reportHelp(const char* argv0) {
    fprintf(stderr, "Usage: %s /path/to/filename|- [-o /path/to/output/filename]\n", argv0);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        reportHelp(argv[0]);
        return -1;
    }
    std::string inputfilepath = argv[1];

    SeqFileType inputfiletype = endsWith(inputfilepath, ".fastq") ? FASTQ : SAM;

    FILE *in = NULL;
    if (inputfilepath.compare("-") == 0) {
        inputfiletype = FASTQ;
        in = stdin;
    }
    else {
        in = fopen(inputfilepath.c_str(), "r");
    }
    if (in == NULL) {
        fprintf(stderr, "Unable to open input file: %s - [%s]\n", inputfilepath.c_str(), strerror(errno));
        return -1;
    }
    FILE *out = stdout;
    if (argc >= 3) {
        std::string cmdopt = argv[2];
        if (cmdopt.compare("-o") != 0) {
            fprintf(stderr, "Invalid command option: %s\n", cmdopt.c_str());
            reportHelp(argv[0]);
            return -1;
        }
        std::string outputfilepath = argv[3];
        out = fopen(outputfilepath.c_str(), "w");
        if (out == NULL) {
            fprintf(stderr, "Unable to open output file: %s - [%s]\n", outputfilepath.c_str(), strerror(errno));
            return -1;
        }
    }
    char line[65536];
    while (fgets(line, sizeof(line), in)) {
        if (inputfiletype == FASTQ) {
            // read the other two lines
            for (int i = 0; i < 2; ++i) {
                if (!fgets(line, sizeof(line), in)) {
                    fprintf(stderr, "Failed to read fastq entry: %s\n", line);
                    return -1;
                }
            }
            // retrieve the qscore line
            if (!fgets(line, sizeof(line), in)) {
                fprintf(stderr, "Failed to read fastq entry: %s\n", line);
                return -1;
            }
            unsigned int idx=0;
            while (!(line[idx] == '\0' || line[idx] == '\n')) {
                fputc(line[idx], out);
                ++idx;
            }
            continue;
        }

        if (line[0] == '@') {
            continue;
        }
        bool wasWhitespace = true;
        unsigned int colNum=0, colStartPos=0, pos=0;
        while (line[pos]!=0) {
            bool whitespace = (line[pos] == ' ' or line[pos] == '\t');
            if (whitespace && !wasWhitespace) {
                wasWhitespace = true;
                if (colNum == 11) {
                    for (unsigned int idx=colStartPos; idx<pos; ++idx)
                        fputc(line[idx], out);
                }
            }
            else if (!whitespace && wasWhitespace) {
                colNum++;
                colStartPos = pos;
                wasWhitespace = false;
            }
            pos++;
        }
    }
    fclose(in);
    if (out != stdout)
        fclose(out);
    return 0;
}
