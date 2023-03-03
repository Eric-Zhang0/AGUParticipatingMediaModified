//
//  Header.h
//  AGUParticipatingMedia
//
//  Created by 張皓珂 on 2021/11/29.
//

#ifndef _RADIANCE_H_
#define _RADIANCE_H_

#include <algorithm>
#include <math.h>
#include <random>
#include <map>

#include "ray.h"
#include "scene.h"
#include "sphere.h"
#include "intersection.h"
#include "random.h"
#include "image.h"
#include "imageiohdr.h"

const double pi = 3.1415926535897932385;

inline double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}
inline double getRadius(const edupt::Sphere& s){
    return s.radius;
}
inline edupt::Vec getCenter(const edupt::Sphere& s){
    return s.position;
}
inline double getMin(const edupt::Ray& r, const edupt::Sphere& s){
    double radius = getRadius(s);
    edupt::Vec center = getCenter(s);
    double direction_length = r.dir.length();
    
    return (center.length() - radius) / direction_length;
}
inline double getMax(const edupt::Ray& r, const edupt::Hitpoint& rec, const edupt::Sphere& s){
    double radius = getRadius(s);
    double direction_length = r.dir.length();
    double t_in = rec.distance;
    
    return (2*radius/direction_length)+t_in;;
}
inline double getObj(const edupt::Ray& r, const edupt::Hitpoint& rec, const edupt::Sphere& s){
    double t_in = rec.distance;
    double radius = getRadius(s);
    edupt::Vec center = getCenter(s);
    edupt::Vec pc_in = (r.org + r.dir * t_in) - center;
    double pc_in_length = pc_in.length();
    double direction_length = r.dir.length();
    
    double cos_angle = dot(r.dir, pc_in) / (pc_in_length*direction_length);
    double t_obj = (2*radius*cos_angle)/direction_length + t_in;
    
    return t_obj;
}
inline double getGrayValue(const edupt::Color& c){
    double grayValue = c.x*0.3 + c.y*0.59 + c.z*0.11;
    return grayValue;
}
inline double transmittance(double t, double sigma){
    return exp(-sigma*t);
}
inline double luminosity(double t, double sigma, const edupt::Color& c){
    return getGrayValue(c)*transmittance(t, sigma)*(1/(4*pi));
}
inline double accept(double t1, double t2, double sigma, const edupt::Color& c){
    double grayValueT1 = luminosity(t1, sigma, c);
    double grayValueT2 = luminosity(t2, sigma, c);
    
    double a = std::min(1.0, grayValueT2/grayValueT1);
    return a;
}
inline double mutate(const edupt::Ray& r, const edupt::Hitpoint& rec, const edupt::Sphere& s, double t){
    auto t_max = getMax(r, rec, s);
    auto t_min = getMin(r, s);
    double epsilon = 0.1;
    double length = t_max - t_min;
    
    t += epsilon * (2 * length * random_double() - length);
    
    if (t > t_max) {
        t -= length;
    }else if(t < t_min){
        t += length;
    }
    return t;
}
namespace edupt {

const Color kBackgroundColor = Color(1.0, 1.0, 1.0);
const int kDepth = 5; // ロシアンルーレットで打ち切らない最大深度
const int kDepthLimit = 2800;

Image<double, 3> bgImage;
bool alreadyLoadHdr = false; // すでにhdrを読み込んでいるか

Color radiance(const Ray& ray, Random* rnd, const int depth, const int n_mutations)
{
    Intersection intersection;
    // ★シーンオブジェクトと交差判定
    if (!intersect_scene(ray, &intersection))
    {
        // シーンの背景にぶつかったとき

        // hdr画像を読み込んでいなければ読み込む
        if (!alreadyLoadHdr)
        {
            loadImageHdr<double, 3>("/Users/eric-zhang/Desktop/c++_source_code/AGParticipatingMedia-master/rnl_probe.hdr", bgImage);

            alreadyLoadHdr = true;
        }
        const double r = (1 / kPI) * acos(ray.dir.z) / sqrt(ray.dir.x * ray.dir.x + ray.dir.y * ray.dir.y);
        double uvR = bgImage.atUVInterpolation((ray.dir.x * r + 1) * 0.5, 1 - (ray.dir.y * r + 1) * 0.5, 0);
        double uvG = bgImage.atUVInterpolation((ray.dir.x * r + 1) * 0.5, 1 - (ray.dir.y * r + 1) * 0.5, 1);
        double uvB = bgImage.atUVInterpolation((ray.dir.x * r + 1) * 0.5, 1 - (ray.dir.y * r + 1) * 0.5, 2);
        return Color(uvR, uvG, uvB);
    }
    const Sphere& next_object = spheres[intersection.object_id];
    const Hitpoint& hitpoint = intersection.hitpoint;
    const Vec orienting_normal = dot(hitpoint.normal, ray.dir) < 0.0 ? hitpoint.normal : (-1.0 * hitpoint.normal); // 交差位置の法線（物体からのレイの入出を考慮）


    // ★ロシアンルーレット関連
    // 色の反射率最大のものを得る。ロシアンルーレットで使う。
    // ロシアンルーレットの閾値は任意だが色の反射率等を使うとより良い。
    double russian_roulette_probability = std::max(next_object.color.x, std::max(next_object.color.y, next_object.color.z));

    // 反射回数が一定以上になったらロシアンルーレットの確率を急上昇させる。（スタックオーバーフロー対策）
    //if (depth > kDepthLimit)
    //    russian_roulette_probability *= pow(0.5, depth - kDepthLimit);

    // ロシアンルーレットを実行し追跡を打ち切るかどうかを判断する。
    // ただしDepth回の追跡は保障する。
    //if (depth > kDepth) {
    if (rnd->next01() >= russian_roulette_probability)
        return next_object.emission;
    //}
    //else
    //    russian_roulette_probability = 1.0; // ロシアンルーレット実行しなかった

    // ★デルタトラッキング距離と面イベント距離の短い方を選択する
    double d_delta_tracking = INFINITY;
    double d_surface = hitpoint.distance;

    double syousan = 0.01;

    
    if (ray.object_id >= 0)
    {
        const Sphere& current_object = spheres[ray.object_id];
        if ( current_object.reflection_type == REFLECTION_TYPE_VOLUME)
        {

        double u0 = rnd->next01();
        d_delta_tracking = -(log(1 - u0) / syousan);
        }
    }

    // ★イベント分岐
    Color incoming_radiance;
    Color weight = 1.0;

    if (d_delta_tracking < d_surface)
    {
        // ★デルタトラッキングイベント
        
        // 新しい座標を決定
        /*Vec newPos = ray.org + d_delta_tracking * ray.dir;

        // ランダムな方向を選ぶ
        Vec newDir = normalize(Vec(rnd->next01() * 2 - 1, rnd->next01() * 2 - 1, rnd->next01() * 2 - 1));
        
        incoming_radiance = radiance(Ray(newPos,newDir, ray.object_id), rnd, depth + 1);
        weight = next_object.color / russian_roulette_probability;*/
        
        // Implementation for PSSMLT:
        const double r = (1 / kPI) * acos(ray.dir.z) / sqrt(ray.dir.x * ray.dir.x + ray.dir.y * ray.dir.y);
        double uvR = bgImage.atUVInterpolation((ray.dir.x * r + 1) * 0.5, 1 - (ray.dir.y * r + 1) * 0.5, 0);
        double uvG = bgImage.atUVInterpolation((ray.dir.x * r + 1) * 0.5, 1 - (ray.dir.y * r + 1) * 0.5, 1);
        double uvB = bgImage.atUVInterpolation((ray.dir.x * r + 1) * 0.5, 1 - (ray.dir.y * r + 1) * 0.5, 2);
        Color initial_color(uvR, uvG, uvB);
        Vec updatePos = ray.org + ray.dir * d_delta_tracking;
        Vec updateDir = ray.dir;
        double t_obj = getObj(ray, hitpoint, next_object);
        
        auto accumulated_delta = d_delta_tracking;

        for (int i = 0; i < n_mutations; i++) {
            double t_c = mutate(ray, hitpoint, next_object, d_delta_tracking);
            double a = accept(d_delta_tracking, t_c, syousan, initial_color);
            double rand = random_double();
            
            if (rand < a){
                if(d_delta_tracking > t_obj) break;
                d_delta_tracking = t_c;
                accumulated_delta += d_delta_tracking;
            }else{
                continue;
            }
        }
        if(accumulated_delta > t_obj) accumulated_delta = t_obj;
        updatePos = ray.org + ray.dir * accumulated_delta;
        updateDir = normalize(Vec(rnd->next01() * 2 - 1, rnd->next01() * 2 - 1, rnd->next01() * 2 - 1));
        incoming_radiance = radiance(Ray(updatePos,updateDir, ray.object_id), rnd, depth + 1, n_mutations);
        weight = next_object.color / russian_roulette_probability;
        return next_object.emission + multiply(weight, incoming_radiance);
    }
    else {
        // ★面イベント

        switch (next_object.reflection_type) {
            // 完全拡散面
        case REFLECTION_TYPE_DIFFUSE: {
            // orienting_normalの方向を基準とした正規直交基底(w, u, v)を作る。この基底に対する半球内で次のレイを飛ばす。
            Vec w, u, v;
            w = orienting_normal;
            if (fabs(w.x) > kEPS) // ベクトルwと直交するベクトルを作る。w.xが0に近い場合とそうでない場合とで使うベクトルを変える。
                u = normalize(cross(Vec(0.0, 1.0, 0.0), w));
            else
                u = normalize(cross(Vec(1.0, 0.0, 0.0), w));
            v = cross(w, u);
            // コサイン項を使った重点的サンプリング
            const double r1 = 2 * kPI * rnd->next01();
            const double r2 = rnd->next01(), r2s = sqrt(r2);
            Vec dir = normalize((
                u * cos(r1) * r2s +
                v * sin(r1) * r2s +
                w * sqrt(1.0 - r2)));

            incoming_radiance = radiance(Ray(hitpoint.position, dir, intersection.object_id), rnd, depth + 1, n_mutations);
            // レンダリング方程式に対するモンテカルロ積分を考えると、outgoing_radiance = weight * incoming_radiance。
            // ここで、weight = (ρ/π) * cosθ / pdf(ω) / R になる。
            // ρ/πは完全拡散面のBRDFでρは反射率、cosθはレンダリング方程式におけるコサイン項、pdf(ω)はサンプリング方向についての確率密度関数。
            // Rはロシアンルーレットの確率。
            // 今、コサイン項に比例した確率密度関数によるサンプリングを行っているため、pdf(ω) = cosθ/π
            // よって、weight = ρ/ R。
            weight = next_object.color / russian_roulette_probability;
        } break;

            // 完全鏡面
        case REFLECTION_TYPE_SPECULAR: {
            // 完全鏡面なのでレイの反射方向は決定的。
            // ロシアンルーレットの確率で除算するのは上と同じ。
            incoming_radiance = radiance(Ray(hitpoint.position, ray.dir - hitpoint.normal * 2.0 * dot(hitpoint.normal, ray.dir),intersection.object_id), rnd, depth + 1, n_mutations);
            weight = next_object.color / russian_roulette_probability;
        } break;

            // 屈折率kIorのガラス
        case REFLECTION_TYPE_REFRACTION: {
            const Ray reflection_ray = Ray(hitpoint.position, ray.dir - hitpoint.normal * 2.0 * dot(hitpoint.normal, ray.dir), intersection.object_id);
            const bool into = dot(hitpoint.normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか

            // Snellの法則
            const double nc = 1.0; // 真空の屈折率
            const double nt = kIor; // オブジェクトの屈折率
            const double nnt = into ? nc / nt : nt / nc;
            const double ddn = dot(ray.dir, orienting_normal);
            const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

            if (cos2t < 0.0) { // 全反射
                incoming_radiance = radiance(reflection_ray, rnd, depth + 1, n_mutations);
                weight = next_object.color / russian_roulette_probability;
                break;
            }
            // 屈折の方向
            const Ray refraction_ray = Ray(hitpoint.position,
                normalize(ray.dir * nnt - hitpoint.normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t))), intersection.object_id);

            // SchlickによるFresnelの反射係数の近似を使う
            const double a = nt - nc, b = nt + nc;
            const double R0 = (a * a) / (b * b);

            const double c = 1.0 - (into ? -ddn : dot(refraction_ray.dir, -1.0 * orienting_normal));
            const double Re = R0 + (1.0 - R0) * pow(c, 5.0); // 反射方向の光が反射してray.dirの方向に運ぶ割合。同時に屈折方向の光が反射する方向に運ぶ割合。
            const double nnt2 = pow(into ? nc / nt : nt / nc, 2.0); // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
            const double Tr = (1.0 - Re) * nnt2; // 屈折方向の光が屈折してray.dirの方向に運ぶ割合

            // 一定以上レイを追跡したら屈折と反射のどちらか一方を追跡する。（さもないと指数的にレイが増える）
            // ロシアンルーレットで決定する。
            const double probability = 0.25 + 0.5 * Re;
            if (depth > 2) {
                if (rnd->next01() < probability) { // 反射
                    incoming_radiance = radiance(reflection_ray, rnd, depth + 1, n_mutations) * Re;
                    weight = next_object.color / (probability * russian_roulette_probability);
                }
                else { // 屈折
                    incoming_radiance = radiance(refraction_ray, rnd, depth + 1, n_mutations) * Tr;
                    weight = next_object.color / ((1.0 - probability) * russian_roulette_probability);
                }
            }
            else { // 屈折と反射の両方を追跡
                incoming_radiance =
                    radiance(reflection_ray, rnd, depth + 1, n_mutations) * Re +
                    radiance(refraction_ray, rnd, depth + 1, n_mutations) * Tr;
                weight = next_object.color / (russian_roulette_probability);
            }
        } break;

        case REFLECTION_TYPE_VOLUME:
        {
            // 新しい座標を決定
            Vec newPos = ray.org + d_surface * ray.dir;

            // ランダムな方向を選ぶ
            Vec newDir = ray.dir;

            incoming_radiance = radiance(Ray(newPos, newDir, intersection.object_id), rnd, depth + 1, n_mutations);
            weight = next_object.color / russian_roulette_probability;

        } break;

        }


        if (!alreadyLoadHdr)
        {
            loadImageHdr<double, 3>("rnl_probe.hdr", bgImage);

            alreadyLoadHdr = true;
        }
        const double r = (1 / kPI) * acos(ray.dir.z) / sqrt(ray.dir.x * ray.dir.x + ray.dir.y * ray.dir.y);
        double uvR = bgImage.atUVNearest((ray.dir.x * r + 1) * 0.5, (ray.dir.y * r + 1) * 0.5, 0);
        double uvG = bgImage.atUVNearest((ray.dir.x * r + 1) * 0.5, (ray.dir.y * r + 1) * 0.5, 1);
        double uvB = bgImage.atUVNearest((ray.dir.x * r + 1) * 0.5, (ray.dir.y * r + 1) * 0.5, 2);


        return next_object.emission + multiply(weight, incoming_radiance);

    }
}
};

#endif
