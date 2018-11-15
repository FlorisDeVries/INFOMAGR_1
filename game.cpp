#include "precomp.h" // include (only) this in every .cpp file

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void Game::Init()
{
	cam = new Camera( vec3( 0, 0, -10 ), vec3( 0, 0, 1 ), tanf( 1.55f ) / ( screen->GetWidth() / 2 ), screen->GetWidth(), screen->GetHeight() );
	primitives[0] = new Sphere( vec3( 0, 0, 0 ), 8 );
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
	for ( int y = 0; y < ylim; y++ )
	{
		for ( int x = 0; x < xlim; x++ )
		{
			// Cast a ray and plot result
			Ray ray = cam->GetRay( x, y );
			for ( int i = 0; i < 1; i++ )
			{
				primitives[i]->Intersect( ray );
				//printf( "Final t %f\n", ray.t );
			}
			if ( ray.t < 1000000.0f )
			{
				*pointer = 0xFFFFFF;
			}
			pointer += 1;
		}
	}
}

Camera::Camera( vec3 pos, vec3 dir, float FOV, int screenWidth, int screenHeight ) : position( pos ), direction( dir ), FOV( FOV ), screenWidth( screenWidth ), screenHeight( screenHeight )
{
	screenCenter = pos + dir * FOV;
}

Ray Tmpl8::Camera::GetRay( int x, int y )
{
	vec3 xinc = ( ScreenCorner( 1 ) - ScreenCorner( 0 ) );
	xinc.x /= ( screenWidth / 2.0f );
	xinc.y /= ( screenWidth / 2.0f );
	xinc.z /= ( screenWidth / 2.0f );
	vec3 yinc = ( ScreenCorner( 2 ) - ScreenCorner( 0 ));
	yinc.x /= screenHeight;
	yinc.y /= screenHeight;
	yinc.z /= screenHeight;
	vec3 rayDirection = vec3( ( ScreenCorner( 0 ) + x * xinc + y * yinc ) - position).normalized();

	return Ray(position, rayDirection);
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

Sphere::Sphere( vec3 pos, float r ) : position( pos ), r2( r )
{
}

void Tmpl8::Sphere::Intersect( Ray ray )
{
	vec3 C = position - ray.origin;
	float t = dot( C, ray.direction );
	vec3 Q = C - t * ray.direction;
	float p2 = dot( Q, Q );
	if ( p2 > r2 ) return; // r2 = r * r
	t -= sqrt( ( r2 - p2 ) );
	printf("New t %f\n", t);
	if ( ( t < ray.t ) && ( t > 0 ) ) ray.t = t;
}

PointLight::PointLight( vec3 col, vec3 pos ) : Light( col, pos )
{
}

Light::Light( vec3 col, vec3 pos ) : color( col ), position( pos )
{
}

Tmpl8::Ray::Ray( vec3 origin, vec3 direction ) : origin( origin ), direction(direction)
{
}
