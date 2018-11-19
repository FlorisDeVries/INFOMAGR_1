#include "precomp.h" // include (only) this in every .cpp file

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void Game::Init()
{
	cam = new Camera( vec3( 0, 0, -20 ), vec3( 0, 0, 1 ), 1.0f / tanf( PI / 4.0f ) );
	primitives.push_back(new Sphere(vec3(0, 0, 0), 4, vec3(0.1f, 0.1f, 1.0f)));
	primitives.push_back(new Sphere(vec3(10, 0, 0), 4, vec3(1.0f, 0.1f, 0.1f)));
	primitives.push_back(new Sphere(vec3(-10, 8, 5), 8, vec3(0.1f, 1.0f, 0.1f)));
	primitives.push_back(new Sphere(vec3(-2, -13, 8), 6, vec3(1.0f, 1.0f, 1.0f)));
	primitives.push_back(new Plane(vec3(0, -1, 0), 10, vec3(1.0f, 1.0f, 1.0f)));

	lights.push_back(new PointLight(vec3(15000000, 15000000, 15000000), vec3(0, 200, 0)));
}

// -----------------------------------------------------------
// Close down application
// -----------------------------------------------------------
void Game::Shutdown()
{
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

	int spherepixels = 0;

	for ( int y = 0; y < ylim; y++ )
	{
		for ( int x = 0; x < xlim; x++ )
		{
			// Cast a ray and plot result
			Intersection intersection;
			Ray ray = cam->GetRay( x, y );
			for ( auto p : primitives )
			{
				p->Intersect( ray , intersection);
			}

			if ( intersection.t < std::numeric_limits<float>::max())
			{
				//Shadow rays

				if (x == 400 && y == 400) {
					printf("middle\n");
				}

				vec3 color = vec3(0, 0, 0);
				for (auto l : lights)
				{
					vec3 or = intersection.position + intersection.normal * 0.001f;
					vec3 dir = l->position - or;
					float maxL = dir.length();
					Ray r(or, dir * (1.0f / maxL));
					Intersection shadow;
					for (auto p: primitives)
					{
						p->Intersect(r, shadow);
						if (shadow.t == 0 || shadow.t < maxL)
							break;
					}
					if (shadow.t == 0 || shadow.t < maxL)
						break;
					vec3 help = intersection.primitive->color;
					color += l->color * help * max(0.0f, dot(intersection.normal, dir * (1.0f / maxL))) * (1.0f / (maxL * maxL));
				}
				spherepixels++;
				uint max = 255;
				uint temp = (((min(max, (uint)color.x)) << 16) & REDMASK) + (((min(max, (uint)color.y)) << 8) & GREENMASK) + ((min(max, (uint)color.z)) & BLUEMASK);
				*pointer = temp;
			}
			else {
				*pointer = 0x6495ED;
			}
			pointer += 1;
		}
	}

	printf("Spherepixels for frame %i: %i (of %i)\n", frame, spherepixels, xlim * ylim);
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
	if (denom > 0.00001f) {
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
