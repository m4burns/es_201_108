/*
 *  g++ -std=c++14 -o mfcc mfcc.cc -lfftw3 -lm
 *  sox in.wav -c 1 -r 8000 -e float -b 64 -t raw - | ./mfcc
 *
 */

#include "mfcc.h"
#include <iostream>

int main() {
  MFCCExtractor<16000> mfcc([](const std::array<double, 14> & features){
      std::cout << features[0];
      for(int i = 1; i < features.size(); i++) {
        std::cout << " " << features[i];
      }
      std::cout << "\n";
    });

  double buf[1024];
  for(int i = 0; i < 1024; i++) { buf[i] = 0.0; }

  while(!std::cin.eof()) {
    std::cin.read((char*)buf, 1024 * sizeof(double));
    mfcc.processSamples(buf, 1024);
  }
}

