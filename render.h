//
//  render.h
//  AGUParticipatingMedia
//
//  Created by 張皓珂 on 2021/11/29.
//

#ifndef _RENDER_H_
#define _RENDER_H_

#include <iostream>

#include "radiance.h"
#include "ppm.h"
#include "random.h"

#include "image.h"
#include "imageiohdr.h"

namespace edupt
{

    int render(const int width, const int height, const int samples, const int supersamples, const int n_mutations)
    {
        // カメラ位置
        const Vec camera_position = Vec(50.0, 52.0, 220.0);
        const Vec camera_dir = normalize(Vec(0.0, -0.04, -1.0));
        const Vec camera_up = Vec(0.0, 1.0, 0.0);

        // ワールド座標系でのスクリーンの大きさ
        const double screen_width = 30.0 * width / height;
        const double screen_height = 30.0;
        // スクリーンまでの距離
        const double screen_dist = 40.0;
        // スクリーンを張るベクトル
        const Vec screen_x = normalize(cross(camera_dir, camera_up)) * screen_width;
        const Vec screen_y = normalize(cross(screen_x, camera_dir)) * screen_height;
        const Vec screen_center = camera_position + camera_dir * screen_dist;

        Color *image = new Color[width * height];

        std::cout << width << "x" << height << " " << samples * (supersamples * supersamples) << " spp" << std::endl;

        Image<double, 3> obj;
        loadImageHdr<double, 3>("rnl_probe.hdr", obj);

        // OpenMP
        #pragma omp parallel for schedule(dynamic, 1)
        for (int y = 0; y < height; y++)
        {
            std::cerr << "Rendering (y = " << y << ") " << (100.0 * y / (height - 1)) << "%" << std::endl;

            Random rnd(y + 1);
            for (int x = 0; x < width; x++)
            {
                const int image_index = (height - y - 1) * width + x;
                // supersamples x supersamples のスーパーサンプリング
                for (int sy = 0; sy < supersamples; sy++)
                {
                    for (int sx = 0; sx < supersamples; sx++)
                    {
                        Color accumulated_radiance = Color();
                        // 一つのサブピクセルあたりsamples回サンプリングする
                        for (int s = 0; s < samples; s++)
                        {
                            const double rate = (1.0 / supersamples);
                            const double r1 = sx * rate + rate / 2.0;
                            const double r2 = sy * rate + rate / 2.0;
                            // スクリーン上の位置
                            const Vec screen_position =
                                screen_center +
                                screen_x * ((r1 + x) / width - 0.5) +
                                screen_y * ((r2 + y) / height - 0.5);
                            // レイを飛ばす方向
                            const Vec dir = normalize(screen_position - camera_position);
                            //const double r = (1 / kPI) * acos(dir.z) / sqrt(dir.x * dir.x + dir.y * dir.y);
                            //double uvR = obj.atUVNearest((dir.x * r + 1)*0.5, (dir.y * r + 1)*0.5, 0);
                            //double uvG = obj.atUVNearest((dir.x * r + 1) * 0.5, (dir.y * r + 1) * 0.5, 1);
                            //double uvB = obj.atUVNearest((dir.x * r + 1) * 0.5, (dir.y * r + 1) * 0.5, 2);
                            //accumulated_radiance = accumulated_radiance + Color(uvR, uvG, uvB) / samples / (supersamples * supersamples);

                             accumulated_radiance = accumulated_radiance +
                                                    radiance(Ray(camera_position, dir, -1), &rnd, 0, n_mutations) / samples / (supersamples * supersamples);
                        }
                        image[image_index] = image[image_index] + accumulated_radiance;
                    }
                }
            }
        }

        // 出力
        int spp = samples * (supersamples * supersamples);
        std::string samplesPerPixel = std::to_string(spp);
        std::string fileName1("image_mixed_");
        std::string fileName2(".ppm");
        std::string fileName = fileName1 + samplesPerPixel + fileName2;
        save_ppm_file(fileName, image, width, height);
        return 0;
    }
};

#endif
