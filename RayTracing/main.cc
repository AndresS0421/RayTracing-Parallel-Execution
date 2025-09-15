//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "rtweekend.h"

#include "camera.h"
#include "hittable.h"
#include "hittable_list.h"
#include "material.h"
#include "sphere.h"
#include <chrono>
#include <iomanip>
#include <thread>
#include <fstream>


int main() {
    hittable_list world;


    auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, ground_material));


    // Create a simple but artistic arrangement - just 8 spheres
    // A small cluster of colorful spheres (moved lower)
    auto red_material = make_shared<lambertian>(color(0.8, 0.2, 0.2));
    world.add(make_shared<sphere>(point3(-1.5, -0.5, 0), 0.5, red_material));

    auto blue_material = make_shared<lambertian>(color(0.2, 0.4, 0.8));
    world.add(make_shared<sphere>(point3(1.5, -0.5, 0), 0.5, blue_material));

    auto green_material = make_shared<lambertian>(color(0.2, 0.8, 0.4));
    world.add(make_shared<sphere>(point3(0, -0.5, -1.5), 0.5, green_material));

    auto yellow_material = make_shared<lambertian>(color(0.8, 0.8, 0.2));
    world.add(make_shared<sphere>(point3(0, -0.5, 1.5), 0.5, yellow_material));


    // Just 3 main spheres - simple but elegant (moved lower)
    auto glass_material = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 0, 0), 1.0, glass_material));

    auto gold_material = make_shared<metal>(color(1.0, 0.8, 0.2), 0.0);
    world.add(make_shared<sphere>(point3(-3, 0, 0), 1.0, gold_material));

    auto silver_material = make_shared<metal>(color(0.7, 0.7, 0.7), 0.1);
    world.add(make_shared<sphere>(point3(3, 0, 0), 1.0, silver_material));

    // Add glowing spheres with the new emissive material (moved lower)
    auto red_glow = make_shared<emissive>(color(1.0, 0.3, 0.3), 2.0);
    world.add(make_shared<sphere>(point3(0, 2, 0), 0.5, red_glow));

    auto blue_glow = make_shared<emissive>(color(0.3, 0.3, 1.0), 1.5);
    world.add(make_shared<sphere>(point3(-2, 1.5, 0), 0.3, blue_glow));

    auto green_glow = make_shared<emissive>(color(0.3, 1.0, 0.3), 1.8);
    world.add(make_shared<sphere>(point3(2, 1.5, 0), 0.3, green_glow));


    camera cam;

    cam.aspect_ratio = 16.0 / 9.0;
    cam.image_width = 1200;
    cam.samples_per_pixel = 25;
    cam.max_depth = 20;

    cam.vfov = 30;
    cam.lookfrom = point3(6, 2, 4);
    cam.lookat = point3(0, 0, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0.0;
    cam.focus_dist = 10.0;


    // Compare sequential vs parallel rendering times
    std::cout << "Starting performance comparison...\n" << std::endl;

    // Sequential rendering
    std::cout << "=== SEQUENTIAL RENDERING ===" << std::endl;
    auto start_seq = std::chrono::high_resolution_clock::now();
    cam.render(world, "output/sequential_render.ppm");
    auto end_seq = std::chrono::high_resolution_clock::now();
    auto duration_seq = std::chrono::duration_cast<std::chrono::seconds>(end_seq - start_seq);
    std::cout << "Sequential rendering time: " << duration_seq.count() << " seconds" << std::endl;

    std::cout << "\n=== PARALLEL RENDERING ===" << std::endl;
    auto start_par = std::chrono::high_resolution_clock::now();
    cam.render_parallel(world, std::thread::hardware_concurrency(), "output/parallel_render.ppm");
    auto end_par = std::chrono::high_resolution_clock::now();
    auto duration_par = std::chrono::duration_cast<std::chrono::seconds>(end_par - start_par);
    std::cout << "Parallel rendering time: " << duration_par.count() << " seconds" << std::endl;

    // Calculate speedup
    double speedup = (double)duration_seq.count() / duration_par.count();
    std::cout << "\n=== PERFORMANCE SUMMARY ===" << std::endl;
    std::cout << "Speedup: " << std::fixed << std::setprecision(2) << speedup << "x faster" << std::endl;
    std::cout << "Time saved: " << (duration_seq.count() - duration_par.count()) << " seconds" << std::endl;
    std::cout << "CPU threads used: " << std::thread::hardware_concurrency() << std::endl;
}