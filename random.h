//
//  random.h
//  AGUParticipatingMedia
//
//  Created by 張皓珂 on 2021/11/29.
//

#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <climits>

namespace edupt {

// Xor-Shift�ɂ�闐���W�F�l���[�^
class XorShift {
    unsigned int seed_[4];
public:
    unsigned int next(void) {
        const unsigned int t = seed_[0] ^ (seed_[0] << 11);
        seed_[0] = seed_[1];
        seed_[1] = seed_[2];
        seed_[2] = seed_[3];
        return seed_[3] = (seed_[3] ^ (seed_[3] >> 19)) ^ (t ^ (t >> 8));
    }

    double next01(void) {
        return (double)next() / UINT_MAX;
    }

    XorShift(const unsigned int initial_seed) {
        unsigned int s = initial_seed;
        for (int i = 1; i <= 4; i++){
            seed_[i-1] = s = 1812433253U * (s^(s>>30)) + i;
        }
    }
};

typedef XorShift Random;

};

#endif
