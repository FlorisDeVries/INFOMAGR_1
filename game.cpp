#include "precomp.h" // include (only) this in every .cpp file

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void Game::Init()
{
	//Setting up the scene
	cam = new Camera(vec3(0, 0, -8), vec3(0, 0, 1), 1.0f / tanf(PI / 4.0f));

	switch (SCENE) {
	case 1:
#pragma region SimpleScene
		// Simple scene
		primitives.push_back(new Sphere(vec3(2, 0, 9), 1.f, vec3(.0f, 1.0f, 0.0f), 0.9f, 0.0f));
		primitives.push_back(new Sphere(vec3(0, 0, 5), 1.5f, vec3(1.f), 0.0f, .0f));
		primitives.push_back(new Plane(vec3(0, -1, 0), 10, vec3(1.f), 0.0f, 0.0f));
		lights.push_back(new PointLight(vec3(LIGHTINTENSITY), vec3(0, 4, 5)));
		lights.push_back(new PointLight(vec3(LIGHTINTENSITY), vec3(0)));
#pragma endregion
		break;

	case 2:
#pragma region ComplexScene
		//Spheres
		primitives.push_back(new Sphere(vec3(0, 0, 0), 2, vec3(1.0f, 0.3f, 0.3f), 0.0f, 1.5f));
		primitives.push_back(new Sphere(vec3(5, 0, 0), 2, vec3(0.3f, 0.3f, 1.0f), 0.1f, 0.0f));
		primitives.push_back(new Sphere(vec3(-5, 4, 2.5f), 3, vec3(0.3f, 1.0f, 0.3f), 0.9f, 0.0f));
		primitives.push_back(new Sphere(vec3(-1, -6.5f, 4), 2.5f, vec3(0.3f, 0.3f, 0.3f), 0.1f, 0.0f));

		//Box
		primitives.push_back(new Plane(vec3(0, -1, 0), 10, vec3(0.3f, 0.3f, 1.0f), 0.0f, 0.0f));
		primitives.push_back(new Plane(vec3(-1, 0, 0), 10, vec3(0.3f, 1.0f, 0.3f), 0.0f, 0.0f));
		primitives.push_back(new Plane(vec3(0, 0, 1), 10, vec3(1.0f, 0.3f, 0.3f), 0.0f, 0.0f));
		primitives.push_back(new Plane(vec3(0, 0, -1), 10, vec3(0.3f, 0.3f, 1.0f), 0.0f, 0.0f));
		primitives.push_back(new Plane(vec3(0, 1, 0), 10, vec3(0.3f, 1.0f, 0.3f), 0.0f, 0.0f));
		primitives.push_back(new Plane(vec3(1, 0, 0), 10, vec3(0.3f, 0.3f, 1.0f), 0.0f, 0.0f));

		//Lights
		lights.push_back(new PointLight(vec3(LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY), vec3(8, 8, 8)));
		lights.push_back(new PointLight(vec3(LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY), vec3(-8, 8, 8)));
		lights.push_back(new PointLight(vec3(LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY), vec3(8, -8, 8)));
		lights.push_back(new PointLight(vec3(LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY), vec3(-8, -8, 8)));
		lights.push_back(new PointLight(vec3(LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY), vec3(8, 8, -8)));
		lights.push_back(new PointLight(vec3(LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY), vec3(8, -8, -8)));
		lights.push_back(new PointLight(vec3(LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY), vec3(-8, 8, -8)));
		lights.push_back(new PointLight(vec3(LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY), vec3(-8, -8, -8)));
		lights.push_back(new PointLight(vec3(10), vec3(0, 0, -8)));
#pragma endregion
		break;
	}


}

// -----------------------------------------------------------
// Close down application
// -----------------------------------------------------------
void Game::Shutdown()
{
}

vec3 Game::Trace(Ray ray, int recursionDepth) {
	Intersection intersection;

	for (auto p : primitives)
	{
		p->Intersect(ray, intersection);
	}

	if (intersection.t < std::numeric_limits<float>::max()) { // Found some primitive
		// Specularity
		if (intersection.primitive->specularity > 0 && recursionDepth > 0) {
			Ray reflectRay = Reflect(ray, intersection);
			vec3 reflectColor = Trace(reflectRay, recursionDepth - 1);
			vec3 directIlluminationColor = DirectIllumination(ray, intersection);
			return intersection.primitive->specularity * reflectColor * directIlluminationColor + (1 - intersection.primitive->specularity) * directIlluminationColor;
		}
		// Refract
		else if (intersection.primitive->refractionIndex > 0) {
			return Refract(ray, intersection, recursionDepth);
		}
		// DirectIllumination
		else {
			return DirectIllumination(ray, intersection);
		}

	}
	else { // Missed all primitives
		return vec3((6.0f * 16 + 4) / 256, (16.0f * 9 + 5) / 256, (16.0f * 14 + 13) / 256); // Cornflower blue
	}
}

Ray Tmpl8::Game::Reflect(Ray & ray, Intersection & intersection)
{
	vec3 direction = (ray.direction - 2 * dot(ray.direction, intersection.normal) * intersection.normal).normalized();
	vec3 origin = intersection.position + intersection.normal * EPSILON;
	return Ray(origin, direction);
}

vec3 Tmpl8::Game::DirectIllumination(Ray & ray, Intersection & intersection)
{
	vec3 color = vec3(0, 0, 0);

	// Shadow rays
	for (auto l : lights)
	{
		vec3 direction = l->position - intersection.position;
		float distance = direction.length();
		direction = direction * (1.f / distance);
		// ez check for spheres
		//if (dot(direction, intersection.normal) < 0) continue;// Works only for spheres

		vec3 origin = intersection.position + direction * EPSILON;
		Ray shadowRay = Ray(origin, direction);

		// Check for obstruction primitives
		bool obstructed = false;
		Intersection shadowIntersect;
		shadowIntersect.t = (l->position - origin).length() - 2 * EPSILON;
		for (auto p : primitives) {
			obstructed = p->Intersect(shadowRay, shadowIntersect);
			if (obstructed)
				break;
		}

		if (obstructed)
			continue;

		color += intersection.primitive->color * l->color * dot(intersection.normal, direction) * (1 / pow(distance, 2));
	}

	return color + vec3(.01f);
}

vec3 Tmpl8::Game::Refract(Ray & ray, Intersection & intersection, int recursionDepth)
{
	return vec3(1);
}

static int frame = 0;

// -----------------------------------------------------------
// Main application tick function
// -----------------------------------------------------------
void Game::Tick(float deltaTime)
{
	//printf("Frame: %i\n", frame);
	if (ONRAILS) {
		//rotate the camera
		float degrees = 2.0f * PI / 120.0f * frame;
		cam->direction.x = cosf(degrees);// *cam->direction.x - sinf(degrees) * cam->direction.y;
		cam->direction.z = sinf(degrees);// *cam->direction.x - cosf(degrees) * cam->direction.y;
		//cam->direction.z = sinf(2.0f * PI / 360 * (frame + 180));
		cam->direction.normalize();
		cam->ResetBounds();
		//printf("Direction: %f, %f, %f\n", cam->direction.x, cam->direction.y, cam->direction.z);
	}

	// clear the graphics window
	screen->Clear(0);
	frame++;
	int xlim = screen->GetWidth(), ylim = screen->GetHeight();
	Pixel *pointer = screen->GetBuffer();

	for (int y = 0; y < ylim; y++)
	{
		for (int x = 0; x < xlim; x++)
		{
			if (x == 400 && y == 400) {
				//printf("middle\n");
			}
			Ray ray = cam->GetRay(x, y);
			vec3 color = Trace(ray, MAX_DEPTH);

			uint max = 255;
			uint red = sqrt(min(1.0f, color.x)) * 255.0f;
			uint green = sqrt(min(1.0f, color.y)) * 255.0f;
			uint blue = sqrt(min(1.0f, color.z)) * 255.0f;

			//uint temp = (((min(max, (uint)color.x)) << 16) & REDMASK) + (((min(max, (uint)color.y)) << 8) & GREENMASK) + ((min(max, (uint)color.z)) & BLUEMASK);
			*pointer = (red << 16) + (green << 8) + (blue);
			pointer += 1;
		}
	}
}

Camera::Camera(vec3 pos, vec3 dir, float FOV) : position(pos), direction(dir), FOV(FOV), screenWidth(screenWidth), screenHeight(screenHeight)
{
	screenCenter = pos + dir * FOV;
	screenTopLeft = ScreenCorner(0);
	xinc = (ScreenCorner(1) - ScreenCorner(0)) * (1.0f / SCRWIDTH);
	yinc = (ScreenCorner(2) - ScreenCorner(0)) * (1.0f / SCRHEIGHT);

}

Ray Tmpl8::Camera::GetRay(int x, int y)
{
	vec3 rayDirection = ((screenTopLeft + x * xinc + y * yinc) - position).normalized();

	Ray r = Ray(position, rayDirection);

	return r;
}

vec3 Tmpl8::Camera::ScreenCorner(int corner)
{
	vec3 left = cross(direction, vec3(0, 1, 0)).normalized();
	vec3 up = cross(left, direction).normalized();

	switch (corner)
	{
	case 0:
		return screenCenter + left + up;
	case 1:
		return screenCenter - left + up;
	case 2:
		return screenCenter + left - up;
	case 3:
		return screenCenter - left - up;

	}
	return vec3();
}

void Tmpl8::Camera::ResetBounds() {
	screenCenter = position + direction * FOV;
	screenTopLeft = ScreenCorner(0);
	xinc = (ScreenCorner(1) - ScreenCorner(0)) * (1.0f / SCRWIDTH);
	yinc = (ScreenCorner(2) - ScreenCorner(0)) * (1.0f / SCRHEIGHT);
}

bool Tmpl8::Sphere::Intersect(Ray &ray, Intersection &intersection)
{
	vec3 C = position - ray.origin;
	float t = dot(C, ray.direction);
	vec3 Q = C - t * ray.direction;
	float p2 = dot(Q, Q);
	if (p2 > r2) return false; // r2 = r * r
	t -= sqrt((r2 - p2));
	if ((t < intersection.t) && (t > 0)) {
		intersection.primitive = this;
		intersection.position = ray.origin + ray.direction * t;
		intersection.normal = (intersection.position - position).normalized();
		intersection.t = t;
		return true;
	}
	return false;
}

//Taken from https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
bool Tmpl8::Plane::Intersect(Ray & ray, Intersection &intersection)
{
	float denom = dot(ray.direction, normal);
	if (denom > EPSILON) {
		vec3 p0l0 = normal * dist - ray.origin;
		float t = dot(normal, p0l0) / denom;
		if (t >= 0 && t < intersection.t) {
			intersection.primitive = this;
			intersection.position = ray.origin + ray.direction * t;
			intersection.normal = 1.0f * normal;
			intersection.t = t;
			return true;
		}
	}
	return false;
}
