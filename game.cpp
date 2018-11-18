#include "precomp.h" // include (only) this in every .cpp file

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void Game::Init()
{
	cam = new Camera( vec3( 0, 0, -20 ), vec3( 0, 0, 1 ), 1.0f / tanf( PI / 4.0f ), 1, 1 );
	primitives[0] = new Sphere(vec3(0, 0, 0), 4);
	primitives[1] = new Sphere(vec3(10, 0, 0), 4);
	primitives[2] = new Plane(vec3(0, -1, 0), 1);
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
			Ray ray = cam->GetRay( x, y );
			for ( int i = 0; i < 3; i++ )
			{
				primitives[i]->Intersect( ray );
			}
			if ( ray.t < 1000.0f)
			{
				spherepixels++;
				*pointer = 0xFFFFFF;
			}
			pointer += 1;
		}
	}

	printf("Spherepixels for frame %i: %i (of %i)\n", frame, spherepixels, xlim * ylim);
}

Camera::Camera( vec3 pos, vec3 dir, float FOV, int screenWidth, int screenHeight ) : position( pos ), direction( dir ), FOV( FOV ), screenWidth( screenWidth ), screenHeight( screenHeight )
{
	screenCenter = pos + dir * FOV;
	screenTopLeft = ScreenCorner(0);
	xinc = ( ScreenCorner( 1 ) - ScreenCorner( 0 ) ) * (1.0f / SCRWIDTH);
	yinc = ( ScreenCorner( 2 ) - ScreenCorner( 0 ) ) * (1.0f / SCRHEIGHT);
	
}

Ray Tmpl8::Camera::GetRay( int x, int y )
{
	if (x == 400 && y == 400) {
		printf("help\n");
	}
	vec3 rayDirection = ( ( screenTopLeft + x * xinc + y * yinc ) - position ).normalized();

	Ray r = Ray(position, rayDirection);
	r.t = 1000.0f;

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

Sphere::Sphere( vec3 pos, float r ) : position( pos ), r2( r * r )
{
}

void Tmpl8::Sphere::Intersect( Ray &ray )
{
	vec3 C = position - ray.origin;
	float t = dot( C, ray.direction );
	vec3 Q = C - t * ray.direction;
	float p2 = dot( Q, Q );
	if ( p2 > r2 ) return; // r2 = r * r
	t -= sqrt( ( r2 - p2 ) );
	if ((t < ray.t) && (t > 0)) { 
		ray.t = t;
	}
}

PointLight::PointLight( vec3 col, vec3 pos ) : Light( col, pos )
{
}

Light::Light( vec3 col, vec3 pos ) : color( col ), position( pos )
{
}

Tmpl8::Ray::Ray( vec3 origin, vec3 direction ) : origin( origin ), direction(direction)
{
	t = 1000.0f;
}

void Tmpl8::Ray::SetT(float haha)
{
	t = haha;
}

Tmpl8::Plane::Plane(vec3 normal, float dist) : normal(normal.normalized()), dist(dist)
{
}

//Taken from https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
void Tmpl8::Plane::Intersect(Ray & ray)
{
	float denom = dot(normal, ray.direction);
	if (denom > 0.00001f) {
		vec3 p0l0 = normal * dist - ray.origin;
		float t = dot(normal, p0l0);
		if (t >= 0 && t < ray.t) 
			ray.t = t;
	}
}
