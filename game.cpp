#include "precomp.h" // include (only) this in every .cpp file

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void Game::Init()
{
	//Setting up the scene
	cam = new Camera( vec3( -6, -3, -8 ), vec3( 0, 0, 1 ), 1.0f / tanf( PI / 4.0f ) );

	switch ( SCENE )
	{
	case 1:
#pragma region SimpleScene
		// Simple scene
		primitives.push_back( new Sphere( vec3( 0, 1, 1 ), 2.f, vec3( 1.0f ), 0.f, 1.5f ) );
		primitives.push_back( new Sphere( vec3( 0, -3.5f, 1 ), 2.f, vec3( 1.f ), 0.0f, 1.54f ) );
		primitives.push_back( new Sphere( vec3( 2, 0, 3 ), 2.f, vec3( 1.f ), .9f, 0.0f ) );
		//primitives.push_back(new Sphere(vec3(-5, 0, 5), 1.5f, vec3(1.f), 1.f, .0f));
		primitives.push_back( new Plane( vec3( 0, -1, 0 ), 5, vec3( 1.f, .2f, .2f ), .0f, 0.0f ) );
		lights.push_back( new PointLight( vec3( LIGHTINTENSITY * 10 ), vec3( 0, 20, 0 ) ) );
		lights.push_back( new PointLight( vec3( LIGHTINTENSITY ), vec3( 2, 0, -2 ) ) );
		//lights.push_back(new PointLight(vec3(LIGHTINTENSITY), vec3(0)));
#pragma endregion
		break;

	case 2:
#pragma region ComplexScene
		//Spheres
		primitives.push_back( new Sphere( vec3( 0, 0, 0 ), 2, vec3( 1.0f, 0.3f, 0.3f ), 0.0f, 1.5f ) );
		primitives.push_back( new Sphere( vec3( 5, 0, 0 ), 2, vec3( 0.3f, 0.3f, 1.0f ), 0.1f, 0.0f ) );
		primitives.push_back( new Sphere( vec3( -5, 4, 2.5f ), 3, vec3( 0.3f, 1.0f, 0.3f ), 0.9f, 0.0f ) );
		primitives.push_back( new Sphere( vec3( -1, -6.5f, 4 ), 2.5f, vec3( 0.3f, 0.3f, 0.3f ), 0.1f, 0.0f ) );

		//Box
		primitives.push_back( new Plane( vec3( 0, -1, 0 ), 10, vec3( 0.3f, 0.3f, 1.0f ), 0.0f, 0.0f ) );
		primitives.push_back( new Plane( vec3( -1, 0, 0 ), 10, vec3( 0.3f, 1.0f, 0.3f ), 0.0f, 0.0f ) );
		primitives.push_back( new Plane( vec3( 0, 0, 1 ), 10, vec3( 1.0f, 0.3f, 0.3f ), 0.0f, 0.0f ) );
		primitives.push_back( new Plane( vec3( 0, 0, -1 ), 10, vec3( 0.3f, 0.3f, 1.0f ), 0.0f, 0.0f ) );
		primitives.push_back( new Plane( vec3( 0, 1, 0 ), 10, vec3( 0.3f, 1.0f, 0.3f ), 0.0f, 0.0f ) );
		primitives.push_back( new Plane( vec3( 1, 0, 0 ), 10, vec3( 0.3f, 0.3f, 1.0f ), 0.0f, 0.0f ) );

		//Lights
		lights.push_back( new PointLight( vec3( LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY ), vec3( 8, 8, 8 ) ) );
		lights.push_back( new PointLight( vec3( LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY ), vec3( -8, 8, 8 ) ) );
		lights.push_back( new PointLight( vec3( LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY ), vec3( 8, -8, 8 ) ) );
		lights.push_back( new PointLight( vec3( LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY ), vec3( -8, -8, 8 ) ) );
		lights.push_back( new PointLight( vec3( LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY ), vec3( 8, 8, -8 ) ) );
		lights.push_back( new PointLight( vec3( LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY ), vec3( 8, -8, -8 ) ) );
		lights.push_back( new PointLight( vec3( LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY ), vec3( -8, 8, -8 ) ) );
		lights.push_back( new PointLight( vec3( LIGHTINTENSITY, LIGHTINTENSITY, LIGHTINTENSITY ), vec3( -8, -8, -8 ) ) );
		lights.push_back( new PointLight( vec3( 10 ), vec3( 0, 0, -8 ) ) );
#pragma endregion
		break;
	case 3:
		primitives.push_back( new Sphere( vec3( -9, -1.f, 1 ), 4.f, vec3( 1.f, .2f, .2f ), 0.0f, 1.54f ) );
		primitives.push_back( new Sphere( vec3( -3, -3.f, 1 ), 2.f, vec3( 1.f, .2f, .2f ), 0.0f, 1.54f ) );
		primitives.push_back( new Sphere( vec3( 0, -4.f, 1 ), 1.f, vec3( 1.f, .2f, .2f ), 0.0f, 1.54f ) );
		primitives.push_back( new Sphere( vec3( 2, -4.5f, 1 ), .5f, vec3( 1.f, .2f, .2f ), 0.0f, 1.54f ) );

		primitives.push_back( new Plane( vec3( 0, -1, 0 ), 5, vec3( 0), .0f, 0.0f ) );

		//lights.push_back( new PointLight( vec3( LIGHTINTENSITY * 10 ), vec3( 0, 20, 0 ) ) );
		//lights.push_back( new PointLight( vec3( LIGHTINTENSITY ), vec3( 2, 0, -2 ) ) );
		break;
	}
}

// -----------------------------------------------------------
// Close down application
// -----------------------------------------------------------
void Game::Shutdown()
{
}

vec3 Game::Trace( Ray ray, int recursionDepth, Intersection &intersection2 )
{
	Intersection intersection;

	


	for ( auto p : primitives )
	{
		p->Intersect( ray, intersection );
	}
	intersection2.t = intersection.t;
	if ( intersection.t < std::numeric_limits<float>::max() )
	{ // Found some primitive
		// Specularity
		if ( intersection.primitive->specularity > 0 && recursionDepth > 0 )
		{
			Ray reflectRay = Reflect( ray, intersection );
			vec3 reflectColor = Trace( reflectRay, recursionDepth - 1);
			vec3 directIllumination = DirectIllumination( ray, intersection );
			return intersection.primitive->GetColor( intersection.position ) * ( intersection.primitive->specularity * reflectColor + ( 1 - intersection.primitive->specularity ) * directIllumination ); //intersection.primitive->GetColor(intersection.position)
		}
		// Refract
		else if ( intersection.primitive->refractionIndex > 0 && recursionDepth > 0 )
		{
			return intersection.primitive->GetColor( intersection.position ) * Refract( ray, intersection, recursionDepth);
		}
		// DirectIllumination
		else
		{
			return intersection.primitive->GetColor( intersection.position ) * DirectIllumination( ray, intersection );
		}
	}
	else
	{																								// Missed all primitives
		return vec3( ( 6.0f * 16 + 4 ) / 256, ( 16.0f * 9 + 5 ) / 256, ( 16.0f * 14 + 13 ) / 256 ); // Cornflower blue
	}
}

Ray Tmpl8::Game::Reflect( Ray &ray, Intersection &intersection )
{
	vec3 direction = ( ray.direction - 2 * dot( ray.direction, intersection.normal ) * intersection.normal ).normalized();
	vec3 origin = intersection.position + intersection.normal * EPSILON;
	return Ray( origin, direction );
}

vec3 Tmpl8::Game::DirectIllumination( Ray &ray, Intersection &intersection )
{
	vec3 color = vec3( 0, 0, 0 );

	// Shadow rays
	for ( auto l : lights )
	{
		vec3 direction = l->position - intersection.position;
		float distance = direction.length();
		direction = direction * ( 1.f / distance );
		//if (dot(direction, intersection.normal) < 0) continue;// Works only for spheres

		vec3 origin = intersection.position + direction * EPSILON;
		Ray shadowRay = Ray( origin, direction );

		// Check for obstruction primitives
		bool obstructed = false;
		Intersection shadowIntersect;
		shadowIntersect.t = ( l->position - origin ).length() - 2 * EPSILON;
		for ( auto p : primitives )
		{
			obstructed = p->Intersect( shadowRay, shadowIntersect );
			if ( obstructed )
				break;
		}

		if ( obstructed )
			continue;

		// Get and combine color
		color += l->color * dot( intersection.normal, direction ) * ( 1 / pow( distance, 2 ) );
	}

	return color + .1f;
}

vec3 Tmpl8::Game::Refract( Ray &ray, Intersection &intersection, int recursionDepth )
{
	if ( recursionDepth <= 0 )
		return vec3( 0 );

	// Prepare some values
	float cosI = clamp( -1.f, 1.f, dot( ray.direction, intersection.normal ) );
	vec3 n = intersection.normal;
	float n1 = 1, n2 = intersection.primitive->refractionIndex;

	bool outside = true;

	if ( cosI > 0 )
	{
		std::swap( n1, n2 );
		n = -intersection.normal;
		outside = false;
	}
	else
	{
		cosI *= -1.f;
	}

	// Calculate the reflect ratio
	float sinT = n1 / n2 * sqrtf( std::max( 0.f, 1 - cosI * cosI ) );
	float refractRatio;
	if ( sinT >= 1 )
	{
		refractRatio = 1;
	}
	else
	{
		float cosT = sqrtf( std::max( 0.f, 1 - sinT * sinT ) );
		cosI = fabsf( cosI );

		float Rs = ( ( n2 * cosI ) - ( n1 * cosT ) ) / ( ( n2 * cosI ) + ( n1 * cosT ) );
		;
		float Rp = ( ( n1 * cosI ) - ( n2 * cosT ) ) / ( ( n1 * cosI ) + ( n2 * cosT ) );
		;
		refractRatio = ( Rs * Rs + Rp * Rp ) / 2;
	}

	// Calculate k
	float n1n2 = n1 / n2;
	float k = 1 - n1n2 * n1n2 * ( 1 - cosI * cosI );

	// If there is some refraction, get refraction color
	vec3 refractColor = 0;
	vec3 offset = EPSILON * intersection.normal;
	if ( refractRatio < 1 )
	{
		vec3 direction = k < 0 ? 0 : n1n2 * ray.direction + ( n1n2 * cosI - sqrtf( k ) ) * n;
		vec3 origin = outside ? intersection.position - offset : intersection.position + offset;
		Intersection refractIntersect;
		refractColor = Trace( Ray( origin, direction ), recursionDepth - 1, refractIntersect );
		if ( !outside )
		{
			float beersLaw = exp( -10e-6 * ( refractIntersect.position - intersection.position ).length() ); // absorption rate * depth
			refractColor *= 1 / ( refractIntersect.position - intersection.position ).length();
			//printf("HA BIER!");
		}
	}

	// Get reflect color
	Ray reflectRay = Reflect( ray, intersection );
	// reflectRay.origin = outside ? intersection.position + offset : intersection.position - offset; -> Offsets? Already included in reflect?
	vec3 reflectColor = Trace( reflectRay, recursionDepth - 1 );

	// Combine reflect and refract
	return reflectColor * refractRatio + refractColor * ( 1 - refractRatio );
}

static int frame = 0;

// -----------------------------------------------------------
// Main application tick function
// -----------------------------------------------------------
void Game::Tick( float deltaTime )
{
	//printf("Frame: %i\n", frame);
	if ( ONRAILS )
	{
		//rotate the camera
		float degrees = 2.0f * PI / 120.0f * frame;
		cam->direction.x = cosf( degrees ); // *cam->direction.x - sinf(degrees) * cam->direction.y;
		cam->direction.z = sinf( degrees ); // *cam->direction.x - cosf(degrees) * cam->direction.y;
		//cam->direction.z = sinf(2.0f * PI / 360 * (frame + 180));
		cam->direction.normalize();
		cam->ResetBounds();
		//printf("Direction: %f, %f, %f\n", cam->direction.x, cam->direction.y, cam->direction.z);
	}

	// clear the graphics window
	screen->Clear( 0 );
	frame++;
	int xlim = screen->GetWidth(), ylim = screen->GetHeight();
	Pixel *pointer = screen->GetBuffer();

	for ( int y = 0; y < ylim; y++ )
	{
		for ( int x = 0; x < xlim; x++ )
		{
			if ( x == 400 && y == 400 )
			{
				//printf("middle\n");
			}
			Ray ray = cam->GetRay( x, y );
			vec3 color = Trace( ray, MAX_DEPTH );

			uint max = 255;
			uint red = sqrt( min( 1.0f, color.x ) ) * 255.0f;
			uint green = sqrt( min( 1.0f, color.y ) ) * 255.0f;
			uint blue = sqrt( min( 1.0f, color.z ) ) * 255.0f;

			//uint temp = (((min(max, (uint)color.x)) << 16) & REDMASK) + (((min(max, (uint)color.y)) << 8) & GREENMASK) + ((min(max, (uint)color.z)) & BLUEMASK);
			*pointer = ( red << 16 ) + ( green << 8 ) + ( blue );
			pointer += 1;
		}
	}
}

Camera::Camera( vec3 pos, vec3 dir, float FOV ) : position( pos ), direction( dir ), FOV( FOV ), screenWidth( screenWidth ), screenHeight( screenHeight )
{
	screenCenter = pos + dir * FOV;
	screenTopLeft = ScreenCorner( 0 );
	xinc = ( ScreenCorner( 1 ) - ScreenCorner( 0 ) ) * ( 1.0f / SCRWIDTH );
	yinc = ( ScreenCorner( 2 ) - ScreenCorner( 0 ) ) * ( 1.0f / SCRHEIGHT );
}

Ray Tmpl8::Camera::GetRay( int x, int y )
{
	vec3 rayDirection = ( ( screenTopLeft + x * xinc + y * yinc ) - position ).normalized();

	Ray r = Ray( position, rayDirection );

	return r;
}

vec3 Tmpl8::Camera::ScreenCorner( int corner )
{
	vec3 left = cross( direction, vec3( 0, 1, 0 ) ).normalized();
	vec3 up = cross( left, direction ).normalized();

	switch ( corner )
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

void Tmpl8::Camera::ResetBounds()
{
	screenCenter = position + direction * FOV;
	screenTopLeft = ScreenCorner( 0 );
	xinc = ( ScreenCorner( 1 ) - ScreenCorner( 0 ) ) * ( 1.0f / SCRWIDTH );
	yinc = ( ScreenCorner( 2 ) - ScreenCorner( 0 ) ) * ( 1.0f / SCRHEIGHT );
}

bool Tmpl8::Sphere::Intersect( Ray &ray, Intersection &intersection )
{
	float t;
	//vec3 C = position - ray.origin;
	//float t = dot( C, ray.direction );
	//vec3 Q = C - t * ray.direction;
	//float p2 = dot( Q, Q );
	//if ( p2 > r2 ) return false; // r2 = r * r
	//t -= sqrt( ( r2 - p2 ) );
	vec3 oc = ray.origin - position;
	float a = dot( ray.direction, ray.direction );
	float b = dot( ray.direction * 2.0f, oc );
	float c = dot( oc, oc ) - r2;
	float disc = sqrt( b * b - 4.0 * a * c );
	float t1 = ( -1.0f * b + disc ) / ( 2.0f * a );
	float t2 = ( -1.0f * b - disc ) / ( 2.0f * a );
	if ( t1 < 0 && t2 < 0 ) return false;
	if ( t1 < 0 )
		t = t2;
	else if ( t2 < 0 )
		t = t1;
	else
		t = min( t1, t2 );
	if ( ( t < intersection.t ) && ( t > 0 ) )
	{
		intersection.primitive = this;
		intersection.position = ray.origin + ray.direction * t;
		intersection.normal = ( intersection.position - position ).normalized();
		intersection.t = t;
		return true;
	}
	return false;
}

vec3 Tmpl8::Sphere::GetColor( vec3 pos )
{
	return color;
}

//Taken from https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
bool Tmpl8::Plane::Intersect( Ray &ray, Intersection &intersection )
{
	float denom = dot( ray.direction, normal );
	if ( denom > EPSILON )
	{
		vec3 p0l0 = normal * dist - ray.origin;
		float t = dot( normal, p0l0 ) / denom;
		if ( t >= 0 && t < intersection.t )
		{
			intersection.primitive = this;
			intersection.position = ray.origin + ray.direction * t;
			intersection.normal = -1.0f * normal;
			intersection.t = t;
			return true;
		}
	}
	return false;
}

vec3 Tmpl8::Plane::GetColor( vec3 pos )
{
	pos += 2000;
	return ( abs( (int)pos.x - 100 ) % 2 > 0 ^ abs( (int)pos.z - 100 ) % 2 > 0 ^ abs( (int)pos.y - 100 ) % 2 > 0 ) ? color : 1;
}
