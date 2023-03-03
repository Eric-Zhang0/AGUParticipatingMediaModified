//
//  ray.h
//  AGUParticipatingMedia
//
//  Created by 張皓珂 on 2021/11/29.
//

#ifndef _RAY_H_
#define _RAY_H_

#include "vec.h"

namespace edupt {

struct Ray {
    Vec org, dir;
    // Ray(const Vec &org, const Vec &dir) : org(org), dir(dir) {}
    Ray(const Vec &_org, const Vec &_dir, const int &_object_id) : org(_org), dir(_dir), object_id(_object_id) {}
    int object_id;
};

};

#endif /* ray_h */
