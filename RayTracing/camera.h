#ifndef CAMERA_H
#define CAMERA_H
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

#include "hittable.h"
#include "material.h"
#include "vec3.h"
#include "ray.h"
#include "color.h"
#include <thread>
#include <mutex>
#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <atomic>


class camera {
public:
    double aspect_ratio = 1.0;  // Ratio of image width over height
    int    image_width = 100;  // Rendered image width in pixel count
    int    samples_per_pixel = 10;   // Count of random samples for each pixel
    int    max_depth = 10;   // Maximum number of ray bounces into scene

    double vfov = 90;              // Vertical view angle (field of view)
    point3 lookfrom = point3(0, 0, 0);   // Point camera is looking from
    point3 lookat = point3(0, 0, -1);  // Point camera is looking at
    vec3   vup = vec3(0, 1, 0);     // Camera-relative "up" direction

    double defocus_angle = 0;  // Variation angle of rays through each pixel
    double focus_dist = 10;    // Distance from camera lookfrom point to plane of perfect focus

    void render(const hittable& world, const std::string& filename = "output.ppm") {
        initialize();

        std::ofstream file(filename);
        file << "P3\n" << image_width << ' ' << image_height << "\n255\n";

        for (int j = 0; j < image_height; j++) {
            std::clog << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
            for (int i = 0; i < image_width; i++) {
                color pixel_color(0, 0, 0);
                for (int sample = 0; sample < samples_per_pixel; sample++) {
                    ray r = get_ray(i, j);
                    pixel_color += ray_color(r, max_depth, world);
                }
                write_color(file, pixel_samples_scale * pixel_color);
            }
        }

        file.close();
        std::clog << "\rDone.                 \n";
    }

    void render_parallel(const hittable& world, int num_threads = std::thread::hardware_concurrency(), const std::string& filename = "output.ppm") {
        initialize();

        std::ofstream file(filename);
        file << "P3\n" << image_width << ' ' << image_height << "\n255\n";

        // Create a vector to store pixel data for each row
		std::vector<std::vector<color>> pixel_data(image_height, std::vector<color>(image_width));

        // Queue-based work distribution
		std::queue<int> row_queue;
		std::mutex queue_mutex;
		std::mutex progress_mutex;
		std::atomic<int> rows_completed(0);

        // Initialize work queue with all row indices
        for (int j = 0; j < image_height; j++) {
            row_queue.push(j);
        }

        // Worker function that pulls work from the queue
        auto worker = [&]() {
            while (true) {
                int row = -1;

                // Get work from the queue
                {
                    std::lock_guard<std::mutex> lock(queue_mutex);
                    if (row_queue.empty()) {
                        return; // No more work
                    }
                    row = row_queue.front();
                    row_queue.pop();
                }

                try {
                    // Render the row
                    for (int i = 0; i < image_width; i++) {
                        color pixel_color(0, 0, 0);

                        for (int sample = 0; sample < samples_per_pixel; sample++) {
                            ray r = get_ray(i, row);
                            pixel_color += ray_color(r, max_depth, world);
                        }

                        pixel_data[row][i] = pixel_samples_scale * pixel_color;
                    }

                    int completed = rows_completed.fetch_add(1) + 1;
                    {
                        std::lock_guard<std::mutex> lock(progress_mutex);
                        std::clog << "\rScanlines remaining: " << (image_height - completed) << ' ' << std::flush;
                    }
                }
                catch (...) {
                    // If rendering fails, put the row back in queue for retry
                    {
                        std::lock_guard<std::mutex> lock(progress_mutex);
                        row_queue.push(row);
                    }
                    // Log error and continue with next work item
                    {
                        std::lock_guard<std::mutex> lock(progress_mutex);
                        std::clog << "\rError processing row " << row << ", retrying... " << std::flush;
                    }
                    // Work with next line
                    continue;
                }
            }
        };

		// Create and launch threads
        std::vector<std::thread> threads;
        for (int t = 0; t < num_threads; t++) {
            threads.emplace_back(worker);
		}

        // Wait for all threads to be completed
        for (auto& thread : threads) {
            thread.join();
        }

        // Output the pixel data in correct order
        for (int j = 0; j < image_height; j++) {
            for (int i = 0; i < image_width; i++) {
                write_color(file, pixel_data[j][i]);
            }
		}

        file.close();

		std::clog << "\rDone.                 \n";
    }

private:
    int    image_height;         // Rendered image height
    double pixel_samples_scale;  // Color scale factor for a sum of pixel samples
    point3 center;               // Camera center
    point3 pixel00_loc;          // Location of pixel 0, 0
    vec3   pixel_delta_u;        // Offset to pixel to the right
    vec3   pixel_delta_v;        // Offset to pixel below
    vec3   u, v, w;              // Camera frame basis vectors
    vec3   defocus_disk_u;       // Defocus disk horizontal radius
    vec3   defocus_disk_v;       // Defocus disk vertical radius

    void initialize() {
        image_height = int(image_width / aspect_ratio);
        image_height = (image_height < 1) ? 1 : image_height;

        pixel_samples_scale = 1.0 / samples_per_pixel;

        center = lookfrom;

        // Determine viewport dimensions.
        auto theta = degrees_to_radians(vfov);
        auto h = std::tan(theta / 2);
        auto viewport_height = 2 * h * focus_dist;
        auto viewport_width = viewport_height * (double(image_width) / image_height);

        // Calculate the u,v,w unit basis vectors for the camera coordinate frame.
        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);

        // Calculate the vectors across the horizontal and down the vertical viewport edges.
        vec3 viewport_u = viewport_width * u;    // Vector across viewport horizontal edge
        vec3 viewport_v = viewport_height * -v;  // Vector down viewport vertical edge

        // Calculate the horizontal and vertical delta vectors from pixel to pixel.
        pixel_delta_u = viewport_u / image_width;
        pixel_delta_v = viewport_v / image_height;

        // Calculate the location of the upper left pixel.
        auto viewport_upper_left = center - (focus_dist * w) - viewport_u / 2 - viewport_v / 2;
        pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

        // Calculate the camera defocus disk basis vectors.
        auto defocus_radius = focus_dist * std::tan(degrees_to_radians(defocus_angle / 2));
        defocus_disk_u = u * defocus_radius;
        defocus_disk_v = v * defocus_radius;
    }

    ray get_ray(int i, int j) const {
        // Construct a camera ray originating from the defocus disk and directed at a randomly
        // sampled point around the pixel location i, j.

        auto offset = sample_square();
        auto pixel_sample = pixel00_loc
            + ((i + offset.x()) * pixel_delta_u)
            + ((j + offset.y()) * pixel_delta_v);

        auto ray_origin = (defocus_angle <= 0) ? center : defocus_disk_sample();
        auto ray_direction = pixel_sample - ray_origin;

        return ray(ray_origin, ray_direction);
    }

    vec3 sample_square() const {
        // Returns the vector to a random point in the [-.5,-.5]-[+.5,+.5] unit square.
        return vec3(random_double() - 0.5, random_double() - 0.5, 0);
    }

    vec3 sample_disk(double radius) const {
        // Returns a random point in the unit (radius 0.5) disk centered at the origin.
        return radius * random_in_unit_disk();
    }

    point3 defocus_disk_sample() const {
        // Returns a random point in the camera defocus disk.
        auto p = random_in_unit_disk();
        return center + (p[0] * defocus_disk_u) + (p[1] * defocus_disk_v);
    }

    color ray_color(const ray& r, int depth, const hittable& world) const {
        // If we've exceeded the ray bounce limit, no more light is gathered.
        if (depth <= 0)
            return color(0, 0, 0);

        hit_record rec;

        if (world.hit(r, interval(0.001, infinity), rec)) {
            // Check if the material is emissive
            auto emissive_mat = dynamic_cast<const emissive*>(rec.mat.get());
            if (emissive_mat) {
                return emissive_mat->emmited();
            }

            ray scattered;
            color attenuation;
            if (rec.mat->scatter(r, rec, attenuation, scattered))
                return attenuation * ray_color(scattered, depth - 1, world);
            return color(0, 0, 0);
        }

        vec3 unit_direction = unit_vector(r.direction());
        auto a = 0.5 * (unit_direction.y() + 1.0);
        return (1.0 - a) * color(1.0, 1.0, 1.0) + a * color(0.5, 0.7, 1.0);
    }
};


#endif