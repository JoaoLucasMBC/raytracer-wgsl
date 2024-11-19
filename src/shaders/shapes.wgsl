fn hit_sphere(center: vec3f, radius: f32, r: ray, record: ptr<function, hit_record>, max: f32)
{
  var a = dot(r.direction, r.direction);
  var h = dot(r.direction, r.origin - center); // h is b divided by 2.0
  var c = dot(r.origin - center, r.origin - center) - radius * radius;

  var delta = h * h - a * c;

  if (delta < 0.0) {
    record.hit_anything = false;
    return;
  }

  var t1 = (-1.0 * h + sqrt(delta)) / a; // t1 is the first intersection point
  var t2 = (-1.0 * h - sqrt(delta)) / a; // t2 is the second intersection point

  var t = min(t1, t2);

  if (t < RAY_TMIN || t > max) {
    record.hit_anything = false;
    return;
  }

  var p = ray_at(r, t);
  var normal = normalize(p - center);
  record.frontface = dot(r.direction, normal) < 0.0; // if the dot is negative, the ray is hitting the outside of the sphere

  record.t = t;
  record.p = p;
  record.normal = mix(-normal, normal, f32(record.frontface)); // if the ray is hitting the inside of the sphere, the normal is inverted
  record.hit_anything = true;
}

fn hit_quad(r: ray, Q: vec4f, u: vec4f, v: vec4f, record: ptr<function, hit_record>, max: f32)
{
  var n = cross(u.xyz, v.xyz);
  var normal = normalize(n);
  var D = dot(normal, Q.xyz);
  var w = n / dot(n.xyz, n.xyz);

  var denom = dot(normal, r.direction);
  if (abs(denom) < 0.0001)
  {
    record.hit_anything = false;
    return;
  }

  var t = (D - dot(normal, r.origin)) / denom;
  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  var intersection = ray_at(r, t);
  var planar_hitpt_vector = intersection - Q.xyz;
  var alpha = dot(w, cross(planar_hitpt_vector, v.xyz));
  var beta = dot(w, cross(u.xyz, planar_hitpt_vector));

  if (alpha < 0.0 || alpha > 1.0 || beta < 0.0 || beta > 1.0)
  {
    record.hit_anything = false;
    return;
  }

  if (dot(normal, r.direction) > 0.0)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = intersection;
  record.normal = normal;
  record.hit_anything = true;
}

fn hit_triangle(r: ray, v0: vec3f, v1: vec3f, v2: vec3f, record: ptr<function, hit_record>, max: f32)
{
  var v1v0 = v1 - v0;
  var v2v0 = v2 - v0;
  var rov0 = r.origin - v0;

  var n = cross(v1v0, v2v0);
  var q = cross(rov0, r.direction);

  var d = 1.0 / dot(r.direction, n);

  var u = d * dot(-q, v2v0);
  var v = d * dot(q, v1v0);
  var t = d * dot(-n, rov0);

  if (u < 0.0 || u > 1.0 || v < 0.0 || (u + v) > 1.0)
  {
    record.hit_anything = false;
    return;
  }

  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = ray_at(r, t);
  record.normal = normalize(n);
  record.hit_anything = true;
}

fn hit_box(r: ray, center: vec3f, rad: vec3f, record: ptr<function, hit_record>, t_max: f32)
{
  var m = 1.0 / r.direction;
  var n = m * (r.origin - center);
  var k = abs(m) * rad;

  var t1 = -n - k;
  var t2 = -n + k;

  var tN = max(max(t1.x, t1.y), t1.z);
  var tF = min(min(t2.x, t2.y), t2.z);

  if (tN > tF || tF < 0.0)
  {
    record.hit_anything = false;
    return;
  }

  var t = tN;
  if (t < RAY_TMIN || t > t_max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = ray_at(r, t);
  record.normal = -sign(r.direction) * step(t1.yzx, t1.xyz) * step(t1.zxy, t1.xyz);
  record.hit_anything = true;

  return;
}

fn hit_cylinder(r: ray, center: vec3f, rad: f32, height_2: f32, record: ptr<function, hit_record>, max: f32) {
  // Precompute some values for optimization
  var Dx = r.direction.x;
  var Dz = r.direction.z;
  var Ox = r.origin.x;
  var Oz = r.origin.z;
  var Cx = center.x;
  var Cz = center.z;

  // Check if it hits the cylindric surface
  var a = Dx * Dx + Dz * Dz;
  var h = Dx * (Ox - Cx) + Dz * (Oz - Cz); // h is b divided by 2.0
  var c = pow(Ox - Cx, 2.0) + pow(Oz - Cz, 2.0) - rad * rad;

  var delta = h * h - a * c;

  if (delta < 0.0) {
    record.hit_anything = false;
    return;
  }

  var t1 = (-1.0 * h + sqrt(delta)) / a; // t1 is the first intersection point
  var t2 = (-1.0 * h - sqrt(delta)) / a; // t2 is the second intersection point

  // Get the closest intersection point
  var t_sup = min(t1, t2);

  // Check if the intersection point is within the height of the cylinder, otherwise discard it
  var p = ray_at(r, t_sup);
  if (p.y < center.y - height_2 || p.y > center.y + height_2) {
    t_sup = 1e8;
  }

  // Check top cap intersection
  var t_top = (center.y + height_2 - r.origin.y) / (r.direction.y);

  var p_top = ray_at(r, t_top);

  // If it is outside the circle, discard it
  if (pow(p_top.x - Cx, 2.0) + pow(p_top.z - Cz, 2.0) > pow(rad, 2.0)) {
    t_top = 1e8;
  }

  // Check bottom cap intersection
  var t_bottom = (center.y - height_2 - r.origin.y) / (r.direction.y);

  var p_bottom = ray_at(r, t_bottom);

  // If it is outside the circle, discard it
  if (pow(p_bottom.x - Cx, 2.0) + pow(p_bottom.z - Cz, 2.0) > (rad * rad)) {
    t_bottom = 1e8;
  }

  // Get the closest intersection point (this code is awful)
  var t = 1e8;
  if (t_sup >= RAY_TMIN && t_sup <= max) {
      t = t_sup;
  }

  if (t_top >= RAY_TMIN && t_top <= max && t_top < t) {
      t = t_top;
  }

  if (t_bottom >= RAY_TMIN && t_bottom <= max && t_bottom < t) {
      t = t_bottom;
  }

  // Check if it is a valid intersection
  if (t < RAY_TMIN || t > max) {
    record.hit_anything = false;
    return;
  }

  p = ray_at(r, t);


  // Choose the normal based on the intersection
  var normal: vec3f;
  if (t == t_sup) {
      // Cylindrical surface
      normal = normalize(vec3f(p.x - Cx, 0.0, p.z - Cz));
  } else if (t == t_top) {
      // Top cap
      normal = vec3f(0.0, 1.0, 0.0);
  } else if (t == t_bottom) {
      // Bottom cap
      normal = vec3f(0.0, -1.0, 0.0);
  }

  record.frontface = dot(r.direction, normal) < 0.0; // if the dot is negative, the ray is hitting the outside of the sphere

  record.t = t;
  record.p = p;
  record.normal = mix(-normal, normal, f32(record.frontface)); // if the ray is hitting the inside of the cylinder, the normal is inverted
  record.hit_anything = true;
}