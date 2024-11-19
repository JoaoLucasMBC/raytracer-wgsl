const THREAD_COUNT = 16;
const RAY_TMIN = 0.0001;
const RAY_TMAX = 100.0;
const PI = 3.1415927f;
const FRAC_1_PI = 0.31830987f;
const FRAC_2_PI = 1.5707964f;

@group(0) @binding(0)  
  var<storage, read_write> fb : array<vec4f>;

@group(0) @binding(1)
  var<storage, read_write> rtfb : array<vec4f>;

@group(1) @binding(0)
  var<storage, read_write> uniforms : array<f32>;

@group(2) @binding(0)
  var<storage, read_write> spheresb : array<sphere>;

@group(2) @binding(1)
  var<storage, read_write> quadsb : array<quad>;

@group(2) @binding(2)
  var<storage, read_write> boxesb : array<box>;

@group(2) @binding(3)
  var<storage, read_write> trianglesb : array<triangle>;

@group(2) @binding(4)
  var<storage, read_write> meshb : array<mesh>;

struct ray {
  origin : vec3f,
  direction : vec3f,
};

struct sphere {
  transform : vec4f,
  color : vec4f,
  material : vec4f,
};

struct quad {
  Q : vec4f,
  u : vec4f,
  v : vec4f,
  color : vec4f,
  material : vec4f,
};

// This can also be a cylinder, if the radius.w is 1.0!
struct box {
  center : vec4f,
  radius : vec4f,
  rotation: vec4f,
  color : vec4f,
  material : vec4f,
};

struct triangle {
  v0 : vec4f,
  v1 : vec4f,
  v2 : vec4f,
};

struct mesh {
  transform : vec4f,
  scale : vec4f,
  rotation : vec4f,
  color : vec4f,
  material : vec4f,
  min : vec4f,
  max : vec4f,
  show_bb : f32,
  start : f32,
  end : f32,
};

struct material_behaviour {
  scatter : bool,
  direction : vec3f,
};

struct camera {
  origin : vec3f,
  lower_left_corner : vec3f,
  horizontal : vec3f,
  vertical : vec3f,
  u : vec3f,
  v : vec3f,
  w : vec3f,
  lens_radius : f32,
};

struct hit_record {
  t : f32,
  p : vec3f,
  normal : vec3f,
  object_color : vec4f,
  object_material : vec4f,
  frontface : bool,
  hit_anything : bool,
};

fn ray_at(r: ray, t: f32) -> vec3f
{
  return r.origin + t * r.direction;
}

fn get_ray(cam: camera, uv: vec2f, rng_state: ptr<function, u32>) -> ray
{
  var rd = cam.lens_radius * rng_next_vec3_in_unit_disk(rng_state);
  var offset = cam.u * rd.x + cam.v * rd.y;
  return ray(cam.origin + offset, normalize(cam.lower_left_corner + uv.x * cam.horizontal + uv.y * cam.vertical - cam.origin - offset));
}

fn get_camera(lookfrom: vec3f, lookat: vec3f, vup: vec3f, vfov: f32, aspect_ratio: f32, aperture: f32, focus_dist: f32) -> camera
{
  var camera = camera();
  camera.lens_radius = aperture / 2.0;

  var theta = degrees_to_radians(vfov);
  var h = tan(theta / 2.0);
  var w = aspect_ratio * h;

  camera.origin = lookfrom;
  camera.w = normalize(lookfrom - lookat);
  camera.u = normalize(cross(vup, camera.w));
  camera.v = cross(camera.u, camera.w);

  camera.lower_left_corner = camera.origin - w * focus_dist * camera.u - h * focus_dist * camera.v - focus_dist * camera.w;
  camera.horizontal = 2.0 * w * focus_dist * camera.u;
  camera.vertical = 2.0 * h * focus_dist * camera.v;

  return camera;
}

fn envoriment_color(direction: vec3f, color1: vec3f, color2: vec3f) -> vec3f
{
  var unit_direction = normalize(direction);
  var t = 0.5 * (unit_direction.y + 1.0);
  var col = (1.0 - t) * color1 + t * color2;

  var sun_direction = normalize(vec3(uniforms[13], uniforms[14], uniforms[15]));
  var sun_color = int_to_rgb(i32(uniforms[17]));
  var sun_intensity = uniforms[16];
  var sun_size = uniforms[18];

  var sun = clamp(dot(sun_direction, unit_direction), 0.0, 1.0);
  col += sun_color * max(0, (pow(sun, sun_size) * sun_intensity));

  return col;
}

fn check_ray_collision(r: ray, max: f32) -> hit_record
{
  var spheresCount = i32(uniforms[19]);
  var quadsCount = i32(uniforms[20]);
  var boxesCount = i32(uniforms[21]);
  var trianglesCount = i32(uniforms[22]);
  var meshCount = i32(uniforms[27]);

  var record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
  var closest = record;

  for (var i = 0; i < spheresCount; i++) {
    var s = spheresb[i];
    var new_record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);

    hit_sphere(s.transform.xyz, s.transform.w, r, &new_record, max);

    if (new_record.hit_anything && new_record.t < closest.t) {
      closest = new_record;
      closest.object_color = s.color;
      closest.object_material = s.material;
    }
  }

  for (var i = 0; i < quadsCount; i++) {
    var q = quadsb[i];
    var new_record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
    hit_quad(r, q.Q, q.u, q.v, &new_record, max);

    if (new_record.hit_anything && new_record.t < closest.t) {
      closest = new_record;
      closest.object_color = q.color;
      closest.object_material = q.material;
    }
  }

  for (var i = 0; i < boxesCount; i++) {
    var b = boxesb[i];
    var new_record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);

    var q = quaternion_from_euler(b.rotation.xyz);
    var q_inv = q_inverse(q);
    var rot_ray = rotate_ray_quaternion(r, b.center.xyz, q);

    // If the radius.w is 0.0, it is a box, otherwise it is a cylinder
    if (b.radius.w == 0.0) {
      hit_box(rot_ray, b.center.xyz, b.radius.xyz, &new_record, max);
    }
    else {
      // For the cylinder, it uses the radius.x and radius.y as the radius and height/2 of the cylinder, respectively
      hit_cylinder(r, b.center.xyz, b.radius.x, b.radius.y, &new_record, max);
    }

    new_record.normal = rotate_vector(new_record.normal, q_inv);
    new_record.p = rotate_vector(new_record.p - b.center.xyz, q_inv) + b.center.xyz;

    if (new_record.hit_anything && new_record.t < closest.t) {
      closest = new_record;
      closest.object_color = b.color;
      closest.object_material = b.material;
    }
  }

  for (var i = 0; i < meshCount; i++) {
    var m = meshb[i];

    var q = quaternion_from_euler(m.rotation.xyz);
    var q_inv = q_inverse(q);
    var rot_ray = rotate_ray_quaternion(r, m.transform.xyz, q);

    // Optimization to avoid checking the mesh if it is not in the view
    if (!AABB_intersect(rot_ray, m.min.xyz * m.scale.xyz + m.transform.xyz, m.max.xyz * m.scale.xyz + m.transform.xyz)) {
      continue;
    }

    // FIX: scale and transform for bbox as well
    if (m.show_bb > 0.0) {
      var new_record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
      var center = (m.min.xyz + m.max.xyz) * 0.5 * m.scale.xyz;
      var radius = (m.max.xyz - m.min.xyz) * 0.5 * m.scale.xyz;
      hit_box(rot_ray, center, radius, &new_record, max);

      if (new_record.hit_anything && new_record.t < closest.t) {
        new_record.normal = rotate_vector(new_record.normal, q_inv);
        new_record.p = rotate_vector(new_record.p - m.transform.xyz, q_inv) + m.transform.xyz;
        closest = new_record;
        closest.object_color = m.color;
        closest.object_material = m.material;
      }
    }

    else {
      for (var j = i32(m.start); j < i32(m.end); j++) {
        var t = trianglesb[j];
        var new_record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
        hit_triangle(rot_ray, 
          (t.v0.xyz * m.scale.xyz+ m.transform.xyz), 
          (t.v1.xyz * m.scale.xyz+ m.transform.xyz), 
          (t.v2.xyz * m.scale.xyz+ m.transform.xyz), &new_record, max);

        if (new_record.hit_anything && new_record.t < closest.t) {
          new_record.normal = rotate_vector(new_record.normal, q_inv);
          new_record.p = rotate_vector(new_record.p - m.transform.xyz, q_inv) + m.transform.xyz;
          closest = new_record;
          closest.object_color = m.color;
          closest.object_material = m.material;
        }
      }
    }
  }
  
  return closest;
}

fn lambertian(normal : vec3f, absorption: f32, random_sphere: vec3f, rng_state: ptr<function, u32>) -> material_behaviour
{
  var res = normal + random_sphere;

  return material_behaviour(true, normalize(res));
}

fn metal(normal : vec3f, direction: vec3f, fuzz: f32, random_sphere: vec3f) -> material_behaviour
{
  var reflected = reflect(direction, normal);
  return material_behaviour(true, normalize(reflected + fuzz * random_sphere));
}

fn dielectric(normal : vec3f, r_direction: vec3f, refraction_index: f32, frontface: bool, random_sphere: vec3f, fuzz: f32, rng_state: ptr<function, u32>) -> material_behaviour
{
  // Ratio of refraction indices
  var eta_over_eta_prime = f32(frontface) * (1.0 / refraction_index) + f32(!frontface) * refraction_index;

  var cos_theta = dot(-1.0 * r_direction, normal);
  var sin_theta = sqrt(1.0 - cos_theta * cos_theta);

  // Check for total internal reflection
  var should_reflect = eta_over_eta_prime * sin_theta > 1.0;

  if (should_reflect) {
    return metal(normal, r_direction, fuzz, random_sphere);
  }

  // Schlick's approximation
  var r0 = (1.0 - refraction_index) / (1.0 + refraction_index);
  r0 = r0 * r0;
  var R = r0 + (1.0 - r0) * pow(1.0 - cos_theta, 5.0);

  // Refraction rays
  var r_perpendicular = eta_over_eta_prime * (r_direction + cos_theta * normal);
  var r_parallel = -1.0 * sqrt(1.0 - dot(r_perpendicular, r_perpendicular)) * normal;
  var new_r = r_perpendicular + r_parallel;
  
  // Randomly choose between reflection and refraction according to the fresnel equation
  var should_refract = rng_next_float(rng_state) > R;
  var r_reflect = reflect(r_direction, normal);

  // Return the new direction
  return material_behaviour(true, normalize(mix(r_reflect, new_r, f32(should_refract))));
}

fn emmisive(color: vec3f, light: f32) -> material_behaviour
{
  // using the direction as the color to make it easier in the future
  return material_behaviour(false, light * color);
}

fn trace(r: ray, rng_state: ptr<function, u32>) -> vec3f
{
  var maxbounces = i32(uniforms[2]);
  var light = vec3f(0.0);
  var color = vec3f(1.0);
  var r_ = r;
  
  var backgroundcolor1 = int_to_rgb(i32(uniforms[11]));
  var backgroundcolor2 = int_to_rgb(i32(uniforms[12]));
  var behaviour = material_behaviour(true, vec3f(0.0));

  for (var j = 0; j < maxbounces; j = j + 1)
  {

    if (!behaviour.scatter) {
      break;
    }

    // check ray collision
    var closest_record = check_ray_collision(r_, RAY_TMAX);

    if (closest_record.hit_anything) {
      // calculate new direction
      var normal = closest_record.normal;

      // tenho um valor 0 ou 1 que eh se o random eh maior que a prob, multiplica isso pelo smoothness e manda bala
      var prob = closest_record.object_material.z;
      var rand = rng_next_float(rng_state);

      var smoothness = closest_record.object_material.x * f32(rand < prob);

      var lamb_behaviour = lambertian(normal, closest_record.object_material.y, rng_next_vec3_in_unit_sphere(rng_state), rng_state);
      var metal_behaviour = metal(normal, r_.direction, closest_record.object_material.y, rng_next_vec3_in_unit_sphere(rng_state));
      var emissive_behaviour = emmisive(closest_record.object_color.xyz, closest_record.object_material.w);
      var dielectric_behaviour = dielectric(normal, r_.direction, closest_record.object_material.z, closest_record.frontface, rng_next_vec3_in_unit_sphere(rng_state), closest_record.object_material.y, rng_state);

      var is_emissive = closest_record.object_material.w > 0;
      var is_dielectric = closest_record.object_material.x < 0;

      // include checking if it is emissive, then kill the ray
      behaviour.scatter = lamb_behaviour.scatter && !is_emissive;

      // LERP between the directions
      if (is_dielectric) {
          behaviour.direction = dielectric_behaviour.direction;
      } else {
          behaviour.direction = mix(lamb_behaviour.direction, metal_behaviour.direction, smoothness);
      }

      var new_dir = behaviour.direction;

      // calculate ray color
      var new_color = mix(closest_record.object_color.xyz, vec3(1.0), smoothness) * f32(!is_emissive) + emissive_behaviour.direction * f32(is_emissive);

      if (is_dielectric) {
        new_color = closest_record.object_color.xyz;
      }

      color *= new_color;

      // calculate new origin
      r_ = ray(closest_record.p, new_dir);
      
      // add light if it is emissive
      light += color * f32(is_emissive);
      
    } else {
      // return background color and object color
      light += envoriment_color(r_.direction, backgroundcolor1, backgroundcolor2) * color;
      break;
    }
  }

  return light;
}

@compute @workgroup_size(THREAD_COUNT, THREAD_COUNT, 1)
fn render(@builtin(global_invocation_id) id : vec3u)
{
    var rez = uniforms[1];
    var time = u32(uniforms[0]);

    // init_rng (random number generator) we pass the pixel position, resolution and frame
    var rng_state = init_rng(vec2(id.x, id.y), vec2(u32(rez)), time);

    // Get uv
    var fragCoord = vec2f(f32(id.x), f32(id.y));
    var uv = (fragCoord + sample_square(&rng_state)) / vec2(rez);

    // Camera
    var lookfrom = vec3(uniforms[7], uniforms[8], uniforms[9]);
    var lookat = vec3(uniforms[23], uniforms[24], uniforms[25]);

    // Get camera
    var cam = get_camera(lookfrom, lookat, vec3(0.0, 1.0, 0.0), uniforms[10], 1.0, uniforms[6], uniforms[5]);
    var samples_per_pixel = i32(uniforms[4]);

    var color = vec3(0.0);

    // Steps:
    // 1. Loop for each sample per pixel
    for (var i = 0; i < samples_per_pixel; i++) {
      // 2. Get ray
      var ray = get_ray(cam, uv, &rng_state);
      // 3. Call trace function
      color += trace(ray, &rng_state);
    }
    // 4. Average the color
    color /= f32(samples_per_pixel);

    // Em algum momento a cor veio negativa??
    color = saturate(color);

    var color_out = vec4(linear_to_gamma(color), 1.0);
    var map_fb = mapfb(id.xy, rez);
    
    // 5. Accumulate the colorr
    var should_accumulate = uniforms[3];

    var accumulated_color = rtfb[map_fb] * should_accumulate + color_out;
    var peso = accumulated_color.w;

    // Set the color to the framebuffer
    rtfb[map_fb] = accumulated_color;
    fb[map_fb] = accumulated_color / peso;
}