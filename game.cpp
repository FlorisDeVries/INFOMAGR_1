#include "precomp.h" // include (only) this in every .cpp file

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void Game::Init()
{
	//Setting up the scene
	cam = new Camera( vec3( 0, 0, -8 ), vec3( 0, 0, 1 ), 1.0f / tanf( PI / 4.0f ) );

	//Spheres
	primitives.push_back(new Sphere(vec3(0, 0, 0), 2, vec3(1.0f, 0.3f, 0.3f), 0.1f));
	primitives.push_back(new Sphere(vec3(5, 0, 0), 2, vec3(0.3f, 0.3f, 1.0f), 0.1f));
	primitives.push_back(new Sphere(vec3(-5, 4, 2.5f), 3, vec3(0.3f, 1.0f, 0.3f), 0.9f));
	primitives.push_back(new Sphere(vec3(-1, -6.5f, 4), 2.5f, vec3(0.3f, 0.3f, 0.3f), 0.1f));

	//Box
	primitives.push_back(new Plane(vec3(0, -1, 0), 10, vec3(0.3f, 0.3f, 1.0f), 0.01f));
	primitives.push_back(new Plane(vec3(-1, 0, 0), 10, vec3(0.3f, 1.0f, 0.3f), 0.01f));
	primitives.push_back(new Plane(vec3(0, 0, 1), 10, vec3(1.0f, 0.3f, 0.3f), 0.01f));
	primitives.push_back(new Plane(vec3(0, 0, -1), 10, vec3(0.3f, 0.3f, 1.0f), 0.01f));
	primitives.push_back(new Plane(vec3(0, 1, 0), 10, vec3(0.3f, 1.0f, 0.3f), 0.01f));
	primitives.push_back(new Plane(vec3(1, 0, 0), 10, vec3(0.3f, 0.3f, 1.0f), 0.01f));

	//Lights
	lights.push_back(new PointLight(vec3(3000, 3000, 3000), vec3(8, 8, 8)));
	lights.push_back(new PointLight(vec3(3000, 3000, 3000), vec3(-8, 8, 8)));
	lights.push_back(new PointLight(vec3(3000, 3000, 3000), vec3(8, -8, 8)));
	lights.push_back(new PointLight(vec3(3000, 3000, 3000), vec3(-8, -8, 8)));
	lights.push_back(new PointLight(vec3(3000, 3000, 3000), vec3(8, 8, -8)));
	lights.push_back(new PointLight(vec3(3000, 3000, 3000), vec3(8, -8, -8)));
	lights.push_back(new PointLight(vec3(3000, 3000, 3000), vec3(-8, 8, -8)));
	lights.push_back(new PointLight(vec3(3000, 3000, 3000), vec3(-8, -8, -8)));
}

// -----------------------------------------------------------
// Close down application
// -----------------------------------------------------------
void Game::Shutdown()
{
}

vec3 Game::Trace(Ray ray, int recursionDepth) {
	vec3 color = vec3(0,0,0);
	// Cast a ray and plot result
	Intersection intersection;
	for (auto p : primitives)
	{
		p->Intersect(ray, intersection);
	}

	if (intersection.t < std::numeric_limits<float>::max())
	{
		//Non-specular part
		//Shadow rays
		for (auto l : lights)
		{
			vec3 or = intersection.position + intersection.normal * EPSILON;
			vec3 dir = l->position - or ;
			float maxL = dir.length();
			Ray r(or , dir * (1.0f / maxL));
			Intersection shadow;
			for (auto p : primitives)
			{
				p->Intersect(r, shadow);
				if (shadow.t == 0 || shadow.t < maxL)
					break;
			}
			if (shadow.t == 0 || shadow.t < maxL)
				continue;
			color += 
				l->color *													// Light color
				intersection.primitive->color *								// Primitive color
				max(0.0f, dot(intersection.normal, dir * (1.0f / maxL))) *	// L dot N
				(1.0f / (maxL * maxL)) *									// Distance attentuation
				(1.0f - intersection.primitive->specularity);				// Non-specularity
		}
		//Specular part
		float specularity = intersection.primitive->specularity;
		if (recursionDepth > 0 && specularity > 0.0f) {
			vec3 or = intersection.position + intersection.normal * EPSILON;
			vec3 dir = (ray.direction - 2 * dot(ray.direction, intersection.normal) * intersection.normal).normalized();

			Ray r = Ray(or , dir);
			vec3 reflectiveColor = Trace(r, recursionDepth - 1);
			color += specularity * reflectiveColor * intersection.primitive->color;
		}
	}
	else {
		//TODO: HDR, other no-hit methods?
		color = vec3(6 * 16 + 4, 16 * 9 + 5, 16 * 14 + 13); // Cornflower blue
	}

	return color;
}

static int frame = 0;

// -----------------------------------------------------------
// Main application tick function
// -----------------------------------------------------------
void Game::Tick( float deltaTime )
{
	// clear the graphics window
	screen->Clear( 0 );
	frame++;
	int xlim = screen->GetWidth(), ylim = screen->GetHeight();
	Pixel *pointer = screen->GetBuffer();

	for ( int y = 0; y < ylim; y++ )
	{
		for ( int x = 0; x < xlim; x++ )
		{
			if (x == 400 && y == 400) {
				//printf("middle\n");
			}
			Ray ray = cam->GetRay(x, y);
			vec3 color = Trace(ray, MAX_DEPTH);

			uint max = 255;
			uint temp = (((min(max, (uint)color.x)) << 16) & REDMASK) + (((min(max, (uint)color.y)) << 8) & GREENMASK) + ((min(max, (uint)color.z)) & BLUEMASK);
			*pointer = temp;
			pointer += 1;
		}
	}
}

Camera::Camera( vec3 pos, vec3 dir, float FOV ) : position( pos ), direction( dir ), FOV( FOV ), screenWidth( screenWidth ), screenHeight( screenHeight )
{
	screenCenter = pos + dir * FOV;
	screenTopLeft = ScreenCorner(0);
	xinc = ( ScreenCorner( 1 ) - ScreenCorner( 0 ) ) * (1.0f / SCRWIDTH);
	yinc = ( ScreenCorner( 2 ) - ScreenCorner( 0 ) ) * (1.0f / SCRHEIGHT);
	
}

Ray Tmpl8::Camera::GetRay( int x, int y )
{
	vec3 rayDirection = ( ( screenTopLeft + x * xinc + y * yinc ) - position ).normalized();

	Ray r = Ray(position, rayDirection);

	return r;
}

vec3 Tmpl8::Camera::ScreenCorner( int corner )
{
	vec3 left = cross( direction, vec3( 0, 1, 0 ) ).normalized();
	vec3 up = cross( left, direction ).normalized();

	switch ( corner )
	{
		case 0 : 
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

void Tmpl8::Sphere::Intersect( Ray &ray, Intersection &intersection)
{
	vec3 C = position - ray.origin;
	float t = dot( C, ray.direction );
	vec3 Q = C - t * ray.direction;
	float p2 = dot( Q, Q );
	if ( p2 > r2 ) return; // r2 = r * r
	t -= sqrt( ( r2 - p2 ) );
	if ((t < intersection.t) && (t > 0)) { 
		intersection.primitive = this;
		intersection.position = ray.origin + ray.direction * t;
		intersection.normal = (intersection.position - position).normalized();
		intersection.t = t;
	}
}

//Taken from https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
void Tmpl8::Plane::Intersect(Ray & ray, Intersection &intersection)
{
	float denom = dot(ray.direction, normal);
	if (denom > EPSILON) {
		vec3 p0l0 = normal * dist - ray.origin;
		float t = dot(normal, p0l0) / denom;
		if (t >= 0 && t < intersection.t) {
			intersection.primitive = this;
			intersection.position = ray.origin + ray.direction * t;
			intersection.normal = -1.0f * normal;
			intersection.t = t;
		}
	}
}
