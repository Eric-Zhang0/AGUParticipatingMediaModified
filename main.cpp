//
//  main.cpp
//  AGUParticipatingMedia
//
//  Created by 張皓珂 on 2021/11/29.
//

#include <iostream>
#include "render.h"
#include "image.h"
#include "imageiohdr.h"

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Path tracing renderer: edupt" << std::endl
                  << std::endl;
    int samplesPerPixel = atoi(argv[1]);
    int n_mutations = atoi(argv[2]);
    std::cout << "mutation number: " << n_mutations << std::endl;
    
        // Image<double, 3> *obj;
        // loadImageHdr<double, 3>("rnl_probe.hdr", *obj);
        // 640x480の画像、(2x2) * 4 sample / pixel
    edupt::render(640, 480, samplesPerPixel, 10, n_mutations);
}
