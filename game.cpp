#include "precomp.h" // include (only) this in every .cpp file
//#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#include "tiny_obj_loader.h"

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void Game::Init()
{
	//Setting up the scene
	cam = new Camera( vec3( 0, 0, -8 ), vec3( 0, 0, 1 ), 4.0f, 1.0f, SCRWIDTH, SCRHEIGHT );

	Surface *planeTexture = new Surface( "assets/Textures/PlaneTexture.jpg" );
	Surface *earth = new Surface( "assets/Textures/earth.jpg" );
	Surface *eye = new Surface( "assets/Textures/oog.jpg" );

	const char *testCubePath = "assets/Obj/Test.obj";
	const char *teapotPath = "assets/Obj/teapot.obj";
	const char *dragonPath = "assets/Obj/dragon.obj";

	switch ( SCENE )
	{
	case 1:
#pragma region SimpleScene
		// Simple scene
		//primitives.push_back( new Sphere( vec3( 0, 1, 1 ), 2.f, vec3( 1.0f ), 0.f, 1.5f ) );
		//primitives.push_back( new Sphere( vec3( 0, -3.5f, 1 ), 2.f, vec3( 1.f ), 0.0f, 1.54f ) );
		//primitives.push_back( new Sphere( vec3( 2, 0, 5 ), 2.f, vec3( 1.f, 0.3f, 0.3f ), .6f, 0.0f ) );
		//primitives.push_back( new Sphere( vec3( -5, -1, 5 ), 1.5f, vec3( 1.f ), 1.f, .0f ) );
		//primitives.push_back(new Sphere(vec3(18, 6, 15), 6.f, vec3(1.f), 0.f, .0f, 1, earth));

		// Light sphere
		primitives.push_back( new Sphere( vec3( 0, 2, 4 ), 4.f, vec3( 1.f ), 0.f, .0f, 1, 0, true ) );

		nonBVHprimitives.push_back( new Plane( vec3( 0, -1, 0 ), 5, vec3( 1.f, .2f, .2f ), .0f, 0.0f, 1, planeTexture ) );
		//nonBVHprimitives.push_back( new Plane( vec3( -1, 0, 0 ), 15, vec3( 1.f, .2f, .2f ), .0f, 0.0f ) );

		//nonBVHprimitives.push_back(new Plane(vec3(0, 1, 0), 10, vec3(1.f, 1.f, 1.f), 0.f, 0.f, 1, 0, true));


		//ReadObj( testCubePath, primitives, vec3( 1.f, .2f, .2f ), vec3( 0 ), .0f );
		//primitives.push_back( new Plane( vec3( 0, 1, 0 ), 5, vec3( 1.f, .2f, .2f ), .0f, 0.0f ) );

		//lights.push_back( new PointLight( vec3( LIGHTINTENSITY * 10 ), vec3( 0, 20, 0 ) ) );
		//lights.push_back( new PointLight( vec3( LIGHTINTENSITY ), vec3( 2, 0, -2 ) ) );
		//lights.push_back( new PointLight( vec3( LIGHTINTENSITY ), vec3( 0 ) ) );
#pragma endregion
		break;

	case 2:
#pragma region ComplexScene
		//Spheres
		primitives.push_back( new Sphere( vec3( 0, 0, 0 ), 2, vec3( 1.0f, 0.3f, 0.3f ), 0.0f, 1.5f, vec3( .5f, 2.3f, 2.3f ) ) );
		primitives.push_back( new Sphere( vec3( 5, 0, 0 ), 2, vec3( 0.3f, 0.3f, 1.0f ), 0.0f, 0.0f, 1, earth ) );
		primitives.push_back( new Sphere( vec3( -5, 4, 2.5f ), 3, vec3( 0.3f, 1.0f, 0.3f ), 0.0f, 0.0f, 1, eye ) );
		primitives.push_back( new Sphere( vec3( -1, -6.5f, 4 ), 2.5f, vec3( 0.3f, 0.3f, 0.3f ), 0.0f, 0.0f, 1, eye ) );

		//Box
		nonBVHprimitives.push_back( new Plane( vec3( 0, -1, 0 ), 10, vec3( 0.3f, 0.3f, 1.0f ), 0.0f, 0.0f ) );
		nonBVHprimitives.push_back( new Plane( vec3( -1, 0, 0 ), 10, vec3( 0.3f, 1.0f, 0.3f ), 0.0f, 0.0f, 1, planeTexture ) );
		nonBVHprimitives.push_back( new Plane( vec3( 0, 0, 1 ), 10, vec3( 1.0f, 0.3f, 0.3f ), 0.0f, 0.0f ) );
		nonBVHprimitives.push_back( new Plane( vec3( 0, 0, -1 ), 10, vec3( 0.3f, 0.3f, 1.0f ), 0.0f, 0.0f ) );
		nonBVHprimitives.push_back( new Plane( vec3( 0, 1, 0 ), 10, vec3( 0.3f, 1.0f, 0.3f ), 0.0f, 0.0f ) );
		nonBVHprimitives.push_back( new Plane( vec3( 1, 0, 0 ), 10, vec3( 0.3f, 0.3f, 1.0f ), 0.0f, 0.0f ) );

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
		// Scene for Beer's law
		primitives.push_back( new Sphere( vec3( -4, -1.f, 1 ), 4.f, vec3( 1.f, .2f, .2f ), 0.0f, 1.54f, vec3( .5f, 2.3f, 2.3f ) ) );
		primitives.push_back( new Sphere( vec3( 2, -3.f, 1 ), 2.f, vec3( 1.f, .2f, .2f ), 0.0f, 1.54f, vec3( .5f, 2.3f, 2.3f ) ) );
		primitives.push_back( new Sphere( vec3( 5, -4.f, 1 ), 1.f, vec3( 1.f, .2f, .2f ), 0.0f, 1.54f, vec3( .5f, 2.3f, 2.3f ) ) );
		primitives.push_back( new Sphere( vec3( 7, -4.5f, 1 ), .5f, vec3( 1.f, .2f, .2f ), 0.0f, 1.54f, vec3( .5f, 2.3f, 2.3f ) ) );

		nonBVHprimitives.push_back( new Plane( vec3( 0, -1, 0 ), 10, vec3( 0.3f, 0.3f, 1.0f ), 0.0f, 0.0f ) );

		lights.push_back( new PointLight( vec3( LIGHTINTENSITY * 10 ), vec3( 0, 20, 0 ) ) );
		lights.push_back( new PointLight( vec3( LIGHTINTENSITY ), vec3( 2, 0, -2 ) ) );
		break;
	case 4:
		// A scene to show the obj loader working, not really interactive
		ReadObj(teapotPath, primitives, vec3(.2f, 1.f, .2f), vec3(0, -1.f, 7));
		ReadObj(teapotPath, primitives, vec3(.2f, 1.f, .2f), vec3(0, -1.f, 14));
		ReadObj(teapotPath, primitives, vec3(.2f, 1.f, .2f), vec3(0, -1.f, 0));
		ReadObj(teapotPath, primitives, vec3(.2f, 1.f, .2f), vec3(7, -1.f, 7));
		ReadObj( teapotPath, primitives, vec3( .2f, 1.f, .2f ), vec3( -7, -1.f, 7 ) );

		printf("Number of primitives: %i\n", primitives.size());

		nonBVHprimitives.push_back( new Plane( vec3( 0, -1, 0 ), 5, vec3( 1.f, .2f, .2f ), .0f, 0.0f, 1 ) );
		lights.push_back( new PointLight( vec3( LIGHTINTENSITY * 10 ), vec3( 3, 10, -5 ) ) );
		break;
	case 5:
		// SAH-testing scene
		ReadObj(dragonPath, primitives, vec3(.4f, 9.f, .2f), vec3(0, -2, 8));
		nonBVHprimitives.push_back(new Plane(vec3(0, -1, 0), 5, vec3(1.f, 1.f, .2f), .4f, 0.0f, 1));

		lights.push_back(new PointLight(vec3(LIGHTINTENSITY * 10), vec3(3, 10, -5)));
		break;
	default:
		break;
	}

#ifdef USE_BVH
	timer t = timer();
	t.reset();
	bvh = BVH();
	bvh.ConstructBVH(&primitives);
	printf("BVH constructed in %f ms\n", t.elapsed());
#endif // USE_BVH
	accumulator[SCRWIDTH * SCRHEIGHT] = { 0 };
}

// -----------------------------------------------------------
// Close down application
// -----------------------------------------------------------
void Game::Shutdown()
{
}

vec3 randomHempsphereReflection(vec3 normal) {
	vec3 randomUnitVector = vec3(2);
	while (randomUnitVector.length() > 1) {
		randomUnitVector.x = ((rand() % 1000) / 1000.0f) * 2 - 1;
		randomUnitVector.y = ((rand() % 1000) / 1000.0f) * 2 - 1;
		randomUnitVector.z = ((rand() % 1000) / 1000.0f) * 2 - 1;
	}
	if (randomUnitVector.dot(normal) < 0) randomUnitVector *= -1.0f;
	return randomUnitVector.normalized();
}

vec3 Game::Sample(Ray r) {
	if (r.recursionDepth == 0)
		return vec3(0);
	Intersection intersection = Intersection();
	Trace(r, r.recursionDepth, intersection, true);
	if (intersection.t == MAXFLOAT) return vec3(0); // ray left scene
	if (intersection.primitive->isLight) { return intersection.primitive->color; printf("Hit light\n"); }  // ray hit light
	vec3 reflectDir = randomHempsphereReflection(intersection.normal);
	Ray reflectRay = Ray(intersection.position + reflectDir * EPSILON, reflectDir);
	reflectRay.recursionDepth = r.recursionDepth - 1;
	vec3 Ei = Sample(reflectRay) * intersection.normal.dot(reflectRay.direction);
	vec3 BRDF = intersection.primitive->color * (1.f / PI);
	return PI * 2.0f * BRDF * Ei;
}

vec3 Game::Trace( Ray ray, int recursionDepth, Intersection &intersection, bool shadowRay )
{
#ifdef USE_BVH
	std::stack<BVHNode> bvhNodeStack;
	bvhNodeStack.push(bvh.root);
	uint depth = 0;
	int rayNodeIntersections = 0;
	while (!bvhNodeStack.empty()) {
		depth++;

		// Get next node
		BVHNode node = bvhNodeStack.top();
		bvhNodeStack.pop();
		float tmin = ray.AABBIntersect(node.bounds);

		if (tmin > intersection.t || tmin < 0 || tmin == std::numeric_limits<float>::max()) continue;

		if (node.isLeaf) {
			//Intersect with primitives
			for (int i = node.first; i < node.first + node.count; i++)
				primitives[i]->Intersect(ray, intersection);
		}
		else {
			//Push child nodes to stack
			if (ray.AABBIntersect(node.left->bounds) > ray.AABBIntersect(node.right->bounds)) {
				bvhNodeStack.push(*node.left);
				bvhNodeStack.push(*node.right);
			}
			else {
				bvhNodeStack.push(*node.right);
				bvhNodeStack.push(*node.left);
			}
		}
}
#ifdef DRAW_DEPTH
	return vec3(min((depth / 10), 256u)/256.0f);
#endif
#else

	for (auto p : primitives)
	{
		p->Intersect(ray, intersection);
	}

#endif
	
	for (auto p : nonBVHprimitives)
	{
		p->Intersect(ray, intersection);
	}

	if (shadowRay) {
		return 0;
	}

	if ( intersection.t < std::numeric_limits<float>::max() )
	{ // Found some primitive
		// Specularity
		if ( intersection.primitive->specularity > 0 && recursionDepth > 0 )
		{
			Ray reflectRay = Reflect( ray, intersection );
			vec3 reflectColor = Trace( reflectRay, recursionDepth - 1 );
			vec3 directIllumination = DirectIllumination( ray, intersection );
			float ratio = intersection.primitive->specularity;
			return reflectColor * ratio + ( 1 - ratio ) * intersection.primitive->color * directIllumination;
		}
		// Refract
		else if ( intersection.primitive->refractionIndex > 0 && recursionDepth > 0 )
		{
			return Refract( ray, intersection, recursionDepth );
		}
		// DirectIllumination
		else
		{
			return intersection.primitive->GetColor( intersection.position ) * DirectIllumination( ray, intersection );
		}
	}
	else
	{
																					// Missed all primitives
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
		if ( dot( direction, intersection.normal ) < 0 ) continue; // Works only for spheres

		vec3 origin = intersection.position + direction * EPSILON;
		Ray shadowRay = Ray( origin, direction );

		// Check for obstruction primitives
		bool obstructed = false;
		Intersection shadowIntersect;
		shadowIntersect.t = ( l->position - origin ).length() - 2 * EPSILON;
		float t = shadowIntersect.t;

		Trace(shadowRay, 0, shadowIntersect, true);		

		if ( shadowIntersect.t < t)
			continue;

		// Get and combine color
		color += l->color * dot( intersection.normal, direction ) * ( 1 / pow( distance, 2 ) );
	}

	return color;
}

vec3 Tmpl8::Game::Refract( Ray &ray, Intersection &intersection, int recursionDepth )
{
	if ( recursionDepth <= 0 )
		return vec3( 0 );

	// Prepare some values
	float cosI = clamp( -1.f, 1.f, dot( ray.direction, intersection.normal ) );
	vec3 n = intersection.normal;
	float n1 = 1, n2 = intersection.primitive->refractionIndex;

	if ( intersection.inside )
	{
		std::swap( n1, n2 );
		n = -intersection.normal;
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
		float Rp = ( ( n1 * cosI ) - ( n2 * cosT ) ) / ( ( n1 * cosI ) + ( n2 * cosT ) );
		refractRatio = ( Rs * Rs + Rp * Rp ) / 2;
	}

	// Calculate k
	float n1n2 = n1 / n2;
	float k = 1 - n1n2 * n1n2 * ( 1 - cosI * cosI );

	// If there is some refraction, get refraction color
	vec3 refractColor = 0;
	if ( refractRatio < 1 )
	{
		vec3 direction = ( k < 0 ? 0 : n1n2 * ray.direction + ( n1n2 * cosI - sqrtf( k ) ) * n ).normalized();
		vec3 origin = intersection.position + ( direction * EPSILON );
		Intersection refractIntersect;
		refractColor = Trace( Ray( origin, direction ), recursionDepth - 1, refractIntersect );

		if ( refractIntersect.inside && refractIntersect.t < std::numeric_limits<float>::max() )
		{
			vec3 color = refractIntersect.primitive->GetColor( refractIntersect.position );
			float r = exp( intersection.primitive->absorptionColor.x * -refractIntersect.t ), g = exp( intersection.primitive->absorptionColor.y * -refractIntersect.t ), b = exp( intersection.primitive->absorptionColor.z * -refractIntersect.t ); // Add rate
			refractColor *= vec3( r, g, b );
		}
	}

	// Get reflect color
	Ray reflectRay = Reflect( ray, intersection );
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
	timer t = timer();
	t.reset();
	// Clear the graphics window
	//screen->Clear( 0 );
	frame++;
	// Reset the tile counter used in the camera for multithreading
	cam->ResetCounter();

	if ( ONRAILS )
	{
		// Rotate the camera
		float degrees = 2.0f * PI / 120.0f * frame;
		cam->direction.x = cosf( degrees );
		cam->direction.z = sinf( degrees );
		cam->direction.normalize();
		cam->ResetBounds();
	}

//Threaded rays using OpenMP
#pragma omp parallel for num_threads( THREADS )
	for (int i = 0; i < THREADS; i++) { ThreadedRays(i); };

	//Copy the accumulator to the screen buffer
	//TODO: Postprocessing?
	uint* pointer = screen->GetBuffer();
	for (int i = 0; i < SCRWIDTH * SCRHEIGHT + 1; i++)
	{
		//printf("%u\n", accumulator[i]);
		*pointer = accumulator[i];
		pointer++;
	}

	HandleInput();

	printf("Frametime (s): %f\n", t.elapsed()/1000.0f);
}

uint EncodeColor(vec3 color) {
	uint red = color.x;
	uint green = color.y;
	uint blue = color.z;

	return (red << 16) + (green << 8) + blue;
}

vec3 DecodeColor(uint color) {
	uint red = (color & REDMASK) >> 16;
	uint green = (color & GREENMASK) >> 8;
	uint blue = color & BLUEMASK;
	
	vec3 colorVec = vec3(red, green, blue);

	return colorVec;
}

vec3 Game::BlendColor(vec3 oldColor, vec3 newColor) {
	vec3 scaledColor = oldColor * (frame - 1);
	vec3 totalColor = scaledColor + newColor;
	vec3 averageColor = totalColor * (1.0f / frame);
	return averageColor;
}

void Game::ThreadedRays( int i )
{
	int tileIdx = 0;

	// While the camera has tiles of rays left to render, get these tiles,
	// trace and render them.
	while ( tileIdx >= 0 )
	{
		//The next set of rays (tiled)
		std::tuple<int, std::vector<Ray>> rayVectorTuple = cam->GetNextRays();

		//When no tiles are left, the tile index passed is -1
		tileIdx = std::get<0>( rayVectorTuple );
		if ( tileIdx < 0 )
			continue;

		//Get the vector of tiles
		std::vector<Ray> rayVector = std::get<1>( rayVectorTuple );

		//Calculate some values for vector and screen iteration

		//The number of tiles in one direction
		int sqTiles = sqrt( TILES );

		//The height and width of a tile
		int xHeight = SCRWIDTH / sqTiles, yHeight = SCRHEIGHT / sqTiles;

		//The starting coordinates of the passed tile
		int xstart = tileIdx % sqTiles * xHeight;
		int ystart = tileIdx / sqTiles * yHeight;

		//The first pixel of the tile passed
		Pixel *pointer = screen->GetBuffer();
		pointer += xstart + ystart * SCRWIDTH;

		uint accumulatorPointer = 0;
		accumulatorPointer += xstart + ystart * SCRWIDTH;

		//Iterate over the tile
		for ( int i = 0; i < yHeight; i++ )
			for ( int j = 0; j < xHeight; j++ )
			{
				//printf("Pixel for pointer %i\n", pointer);
				//vec3 color = Trace( rayVector[i * yHeight + j], MAX_DEPTH );
				vec3 color = Sample(rayVector[i * yHeight + j]);

				uint red = sqrt( min( 1.0f, color.x ) ) * 255.0f;
				uint green = sqrt( min( 1.0f, color.y ) ) * 255.0f;
				uint blue = sqrt( min( 1.0f, color.z ) ) * 255.0f;

				vec3 newColorVec = vec3(red, green, blue);

				// Using the screen buffer
				uint oldColor = *pointer;
				vec3 oldColorVec = DecodeColor(oldColor);

				newColorVec += oldColorVec * (frame - 1);
				newColorVec *= (1.0f / frame);

				*pointer = EncodeColor(newColorVec);
				pointer += ( j + 1 ) % xHeight == 0 ? ( SCRWIDTH - xHeight + 1 ) : 1;

				//Using the accumulator
				uint oldColor2 = accumulator[accumulatorPointer];
				vec3 oldColor2Vec = DecodeColor(oldColor2);

				uint newColor2 = EncodeColor((vec3(red, green, blue) + oldColor2Vec * (frame - 1)) * (1.0f / frame));
				accumulator[accumulatorPointer] = newColor2;
				//printf("acc %u\n", accumulatorPointer);
				accumulatorPointer += (j + 1) % xHeight == 0 ? (SCRWIDTH - xHeight + 1) : 1;

				//if (frame > 1 && oldColorVec.length() > 0) {
				//	//printf("old color: %f %f %f\n", oldColorVec.x, oldColorVec.y, oldColorVec.z);
				//}
				//vec3 sumColor = (oldColorVec * (frame - 1) + newColorVec) * (1.0f / frame);

				//red = sqrt(min(1.0f, oldColorVec.x)) * 255.0f;
				//green = sqrt(min(1.0f, oldColorVec.y)) * 255.0f;
				//blue = sqrt(min(1.0f, oldColorVec.z)) * 255.0f;
				
			}
	}
	//	printf("Tile done!\n");
}

void Game::KeyUp( int key )
{
	//Forward, left, backward, right, up, down
	if ( key == SDL_SCANCODE_W )
		isWDown = false;
	else if ( key == SDL_SCANCODE_A )
		isADown = false;
	else if ( key == SDL_SCANCODE_S )
		isSDown = false;
	else if ( key == SDL_SCANCODE_D )
		isDDown = false;
	else if ( key == SDL_SCANCODE_SPACE )
		isRDown = false;
	else if ( key == SDL_SCANCODE_LSHIFT )
		isFDown = false;
	// Increase/decrease FOV
	else if ( key == SDL_SCANCODE_Y )
		isYDown = false;
	else if ( key == SDL_SCANCODE_H )
		isHDown = false;
}

void Game::KeyDown( int key )
{
	//Forward, left, backward, right, up, down
	if ( key == SDL_SCANCODE_W )
		isWDown = true;
	else if ( key == SDL_SCANCODE_A )
		isADown = true;
	else if ( key == SDL_SCANCODE_S )
		isSDown = true;
	else if ( key == SDL_SCANCODE_D )
		isDDown = true;
	else if ( key == SDL_SCANCODE_SPACE )
		isRDown = true;
	else if ( key == SDL_SCANCODE_LSHIFT )
		isFDown = true;
	// Increase/decrease FOV
	else if ( key == SDL_SCANCODE_Y )
		isYDown = true;
	else if ( key == SDL_SCANCODE_H )
		isHDown = true;
}

void Game::HandleInput()
{
	vec3 forward = cam->direction;
	vec3 left = cross( forward, vec3( 0, 1, 0 ) ).normalized();
	vec3 up = cross( forward, left ).normalized();

	vec3 translation = 0;

	//Forward, left, backward, right, up, down
	if ( isWDown ) translation += cam->direction * MOVEMENTRATE;
	if ( isADown ) translation += left * MOVEMENTRATE;
	if ( isSDown ) translation += -cam->direction * MOVEMENTRATE;
	if ( isDDown ) translation += -left * MOVEMENTRATE;
	if ( isRDown ) translation += -up * MOVEMENTRATE;
	if ( isFDown ) translation += up * MOVEMENTRATE;

	cam->position += translation;

	// FOV
	if ( isYDown && cam->FOV > 2.2f ) cam->FOV -= 0.1f;
	if ( isHDown ) cam->FOV += 0.1f;

	cam->ResetFOV();
	cam->ResetBounds();
}

Camera::Camera( vec3 pos, vec3 dir, float FOV, float aspectRatio, int screenWidth, int screenHeight ) : position( pos ), direction( dir ), FOV( FOV ), aspectRatio( aspectRatio ), screenWidth( screenWidth ), screenHeight( screenHeight )
{
	ResetFOV();
	ResetBounds();
}

bool Tmpl8::Game::ReadObj( const char *path, std::vector<Primitive *> &primitives, vec3 color, vec3 position, float specularity )
{
#ifdef TINYOBJLOADER_IMPLEMENTATION
	// tinyObj implementation, seems to be unable to parse files with "f v v v" format
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;

	std::string warn;
	std::string err;
	bool ret = tinyobj::LoadObj( &attrib, &shapes, &materials, &warn, &err, path );

	if ( !err.empty() )
	{ // `err` may contain warning message.
		std::cerr << err << std::endl;
	}

	if ( !ret )
	{
		exit( 1 );
	}

	vec3 normal;

	int count = 0;
	// Loop over shapes
	for ( size_t s = 0; s < shapes.size(); s++ )
	{
		// Loop over faces(polygon)
		size_t index_offset = 0;
		for ( size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++ )
		{
			int fv = shapes[s].mesh.num_face_vertices[f];

			std::vector<vec3 *> temp_vertices;
			std::vector<vec3 *> temp_normals;

			// Loop over vertices in the face.
			for ( size_t v = 0; v < fv; v++ )
			{
				// access to vertex
				tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
				tinyobj::real_t vx = attrib.vertices[3 * idx.vertex_index + 0];
				tinyobj::real_t vy = attrib.vertices[3 * idx.vertex_index + 1];
				tinyobj::real_t vz = attrib.vertices[3 * idx.vertex_index + 2];
				tinyobj::real_t nx = attrib.normals[3 * idx.normal_index + 0];
				tinyobj::real_t ny = attrib.normals[3 * idx.normal_index + 1];
				tinyobj::real_t nz = attrib.normals[3 * idx.normal_index + 2];
				tinyobj::real_t tx = attrib.texcoords[2 * idx.texcoord_index + 0];
				tinyobj::real_t ty = attrib.texcoords[2 * idx.texcoord_index + 1];

				temp_vertices.push_back( new vec3( vx, vy, vz ) );
				temp_normals.push_back( new vec3( nx, ny, nz ) );
				// Optional: vertex colors
				// tinyobj::real_t red = attrib.colors[3*idx.vertex_index+0];
				// tinyobj::real_t green = attrib.colors[3*idx.vertex_index+1];
				// tinyobj::real_t blue = attrib.colors[3*idx.vertex_index+2];
			}
			index_offset += fv;

			if ( fv == 3 )
			{
				vec3 vx = *temp_vertices[0];
				vec3 vy = *temp_vertices[1];
				vec3 vz = *temp_vertices[2];

				vec3 v1v0 = vy - vx;
				vec3 v2v0 = vz - vx;
				normal = cross( v1v0, v2v0 );

				if ( dot( normal, *temp_normals[0] ) < 0 )
					normal = -normal;

				primitives.push_back( new Triangle( vx + position, vy + position, vz + position, normal, color, specularity ) );
			}

			// per-face material
			shapes[s].mesh.material_ids[f];
		}
	}
	return true;
#endif // TINYOBJLOADER_IMPLEMENTATION
	//http://www.opengl-tutorial.org/beginners-tutorials/tutorial-7-model-loading/

	std::vector<unsigned int> vertexIndices;
	std::vector<vec3> temp_vertices;

	FILE *file = fopen( path, "r" );
	if ( file == NULL )
	{
		printf( "Failed to read file %s!\n", path );
		return false;
	}

	while ( 1 )
	{
		char lineHeader[128];

		// Readline
		int res = fscanf( file, "%s", lineHeader );
		if ( res == EOF )
		{
			break; // End of file, break
		}
		if ( strcmp( lineHeader, "v" ) == 0 )
		{
			vec3 vertex;
			fscanf( file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z );
			temp_vertices.push_back( vertex );
		}
		else if ( strcmp( lineHeader, "f" ) == 0 )
		{
			std::string vertex1, vertex2, vertex3;
			unsigned int vertexIndex[3];
			int matches = fscanf( file, "%d %d %d\n", &vertexIndex[0], &vertexIndex[1], &vertexIndex[2] );
			if ( matches != 3 )
			{
				printf( "File can't be read by our simple parser : ( Try exporting with other options\n" );
				return false;
			}
			vertexIndices.push_back( vertexIndex[0] );
			vertexIndices.push_back( vertexIndex[1] );
			vertexIndices.push_back( vertexIndex[2] );

			vec3 vx = temp_vertices[vertexIndex[0] - 1];
			vec3 vy = temp_vertices[vertexIndex[1] - 1];
			vec3 vz = temp_vertices[vertexIndex[2] - 1];

			vec3 v1v0 = vy - vx;
			vec3 v2v0 = vz - vx;
			vec3 normal = cross( v1v0, v2v0 );

			//if (dot(normal, *temp_normals[0]) < 0)
			//	normal = -normal;

			primitives.push_back( new Triangle( vx + position, vy + position, vz + position, normal, color, specularity ) );
		}
	}
	return true;
}

Ray Tmpl8::Camera::GetRay( int x, int y )
{
	vec3 rayDirection = ( ( screenTopLeft + x * xinc + y * yinc ) - position ).normalized();

	Ray r = Ray( position, rayDirection );

	return r;
}

// Returns a tuple containing the tile index
// and a vector of the rays for that tile
std::tuple<int, std::vector<Ray>> Tmpl8::Camera::GetNextRays()
{
	// <rayCounter> is an atomic int
	int tileIdx = rayCounter++;
	std::vector<Ray> rayVector;

	// If the counter is past the number of tiles, let the thread know it can stop
	// requesting tiles
	if ( tileIdx >= TILES ) return std::make_tuple( -1, rayVector );

	// Some helper numbers for tile traversal
	int sqTiles = sqrt( TILES );
	int xHeight = screenWidth / sqTiles, yHeight = screenHeight / sqTiles;

	int xstart = tileIdx % sqTiles * xHeight;
	int ystart = tileIdx / sqTiles * yHeight;

	for ( int j = 0; j < yHeight; j++ )
	{
		for ( int k = 0; k < xHeight; k++ )
		{
			Ray ray = GetRay( xstart + k, ystart + j );
			ray.recursionDepth = MAX_DEPTH;
			rayVector.push_back( ray );
		}
	}

	return std::make_tuple( tileIdx, rayVector );
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
	screenCenter = position + direction * FOV_Distance;
	screenTopLeft = ScreenCorner( 0 );
	xinc = ( ScreenCorner( 1 ) - ScreenCorner( 0 ) ) * ( 1.0f / SCRWIDTH ) * aspectRatio;
	yinc = ( ScreenCorner( 2 ) - ScreenCorner( 0 ) ) * ( 1.0f / SCRHEIGHT ) * ( 1.0f / aspectRatio );
}

void Tmpl8::Camera::ResetFOV()
{
	FOV_Distance = 1.0f / tanf( PI / FOV );
}

void Tmpl8::Camera::ResetCounter()
{
	rayCounter = 0;
}

bool Tmpl8::Sphere::Intersect( Ray &ray, Intersection &intersection )
{
	vec3 oc = position - ray.origin;
	float a = dot( ray.direction, ray.direction );
	float b = dot( ray.direction, oc );
	float c = dot( oc, oc ) - b * b;
	if ( c > r2 )
		return false;
	float disc = sqrt( r2 - c );
	float t = b - disc;
	float t2 = b + disc;
	bool inside = t < 0;
	if ( inside )
	{
		t = t2;
	}

	if ( t < 0 )
		t = t2;

	if ( ( t < intersection.t ) && ( t > 0 ) )
	{
		intersection.primitive = this;
		intersection.position = ray.origin + ray.direction * t;
		intersection.normal = ( intersection.position - position ).normalized();
		intersection.t = t;
		intersection.inside = inside;
		return true;
	}
	return false;
}

vec3 Tmpl8::Sphere::GetColor( vec3 pos )
{
	if ( !texture )
	{
		return color;
	}
	else
	{

		// http://www.cs.unc.edu/~rademach/xroads-RT/RTarticle.html
		vec3 Vn = vec3( 0, 1, 0 );				   // Vector towards the northpole
		vec3 Ve = vec3( 0, 0, -1 ).normalized();   // Vector towards the equator, use this to rotate spheres
		vec3 Vp = ( position - pos ).normalized(); // Vector towards the intersection point

		float phi = acosf( -dot( Vn, Vp ) );
		float v = phi / PI; // Value between 0 and 1

		float theta = ( acosf( dot( Vp, Ve ) / sin( phi ) ) ) / ( 2 * PI );
		float u;
		if ( dot( cross( Vn, Ve ), Vp ) < 0 )
			u = theta; // Also between 0 and 1
		else
			u = 1 - theta;

		int h = texture->GetHeight();
		int w = texture->GetWidth();
		int x = u * w;
		int y = v * h;

		// Retrieve the pixel and transform it into a vec3 color
		Pixel pixel = texture->GetBuffer()[x + y * w];
		BYTE r = ( pixel & REDMASK ) >> 16;
		BYTE g = ( pixel & GREENMASK ) >> 8;
		BYTE b = ( pixel & BLUEMASK );
		vec3 color = vec3( r / 255.f, g / 255.f, b / 255.f );

		return color;
	}
}

aabb Tmpl8::Sphere::GetBounds()
{
	return aabb(position - sqrtf(r2), position + sqrtf(r2));
}

Tmpl8::Plane::Plane( vec3 normal, float dist, vec3 color, float specularity, float refractionIndex, vec3 absorptionColor, Surface *texture, bool isLight) : normal( normal.normalized() ), dist( dist ), Primitive( color, specularity, refractionIndex, absorptionColor, texture, isLight )
{
	UAxis = vec3( normal.y, normal.z, -normal.x );
	VAxis = UAxis.cross( normal );
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
	if ( !texture )
	{
		pos += 2000;
		return ( abs( (int)pos.x - 100 ) % 2 > 0 ^ abs( (int)pos.z - 100 ) % 2 > 0 ^ abs( (int)pos.y - 100 ) % 2 > 0 ) ? color : 1;
	}
	else
	{
		// Calculate u and v
		int h = texture->GetHeight();
		int w = texture->GetWidth();

		int x = abs( (int)( dot( pos, UAxis ) * w ) % w );
		int y = abs( (int)( dot( pos, VAxis ) * h ) % h );

		// Retrieve the pixel and transform it into a vec3 color
		Pixel pixel = texture->GetBuffer()[x + y * w];
		BYTE r = ( pixel & REDMASK ) >> 16;
		BYTE g = ( pixel & GREENMASK ) >> 8;
		BYTE b = ( pixel & BLUEMASK );
		vec3 color = vec3( r / 255.f, g / 255.f, b / 255.f );

		return color;
	}
}

aabb Tmpl8::Plane::GetBounds()
{
	throw new logic_error( 0 );
}

vec3 Tmpl8::Plane::GetCenter()
{
	throw new logic_error(0);
}

bool Tmpl8::Triangle::Intersect( Ray &ray, Intersection &intersection )
{
	vec3 v1v0 = v1 - v0;
	vec3 v2v0 = v2 - v0;
	vec3 rov0 = ray.origin - v0;

	vec3 q = cross( rov0, ray.direction );
	float d = 1.0 / dot( ray.direction, normal );
	float u = d * dot( -q, v2v0 );
	float v = d * dot( q, v1v0 );
	float t = d * dot( -normal, rov0 );

	if ( u < 0.0 || u > 1.0 || v < 0.0 || ( u + v ) > 1.0 || !( t >= 0 && t < intersection.t ) )
		return false;

	intersection.t = t;
	intersection.normal = normal;
	intersection.position = ray.origin + ray.direction * t;
	intersection.primitive = this;
	return true;
}

vec3 Tmpl8::Triangle::GetColor( vec3 pos )
{
	return color;
}

aabb Tmpl8::Triangle::GetBounds()
{
	return aabb(vec3::min(v2, vec3::min(v0, v1)), vec3::max(v2, vec3::max(v0, v1)));
}

vec3 Tmpl8::Triangle::GetCenter()
{
	return (v0 + v1 + v2) * 0.33333f;
}

void Tmpl8::BVH::ConstructBVH( std::vector<Primitive *> * primitives )
{
	// allocate BVH root node
	root = BVHNode();
	// subdivide root node
	root.first = 0;
	root.count = primitives->size();
	root.Subdivide(primitives);
}

void Tmpl8::BVHNode::Subdivide(std::vector<Primitive *> * primitives)
{
	left = new BVHNode();
	right = new BVHNode();
	bool subFurther = Partition(primitives);
	if (subFurther) {
		left->Subdivide(primitives);
		right->Subdivide(primitives);
		isLeaf = false;
	}
}

// Sort primitives for a given axis
bool sortByAxis(Primitive* prim1, Primitive* prim2, int axis) {
	return prim1->GetCenter()[axis] < prim2->GetCenter()[axis];
}

bool Tmpl8::BVHNode::Partition(std::vector<Primitive *> * primitives)
{
	// 1. Pick an axis
	// 2. Sort all the primitives in the current node on axis coordinate
	// 3. Create n bins
	// 4. Calculate AABB's for bins
	// 5. Create two AABB's
	// 6. Iteratively add the bin AABB's to the left/right box
	// 7. Keep track of lowest split cost
	// 8. Repeat for other axes

	int bestSplitAxis = 0;
	float bestSplitPosition = 0;
	float bestSplitScore = MAXFLOAT;

	// Get the current cost of the node bounds
	bounds = aabb(vec3(MAXFLOAT), vec3(MINFLOAT));
	bounds.Reset();
	for (int i = first; i < first + count; i++)
	{
		bounds.bmin3 = vec3::min(bounds.bmin3, primitives[0][i]->GetBounds().bmin3);
		bounds.bmax3 = vec3::max(bounds.bmax3, primitives[0][i]->GetBounds().bmax3);
	}
	float currentSplitScore = count * bounds.Area();

	// Iterate through axes
	for (int i = 0; i < 3; i++)
	{
		// Sort
		std::sort(primitives->begin() + first, primitives->begin() + first + count, std::bind(sortByAxis, std::placeholders::_1, std::placeholders::_2, i));
		
		// Binning
		float interval = bounds.Extend(i) / min(100.0f, (float)count);
		aabb* bins = new aabb[min(count, 100)];

		int* binCounts = new int[min(count, 100) + 1]{0};
		int counter = first;
		for (int j = 0; j < min(100, count); j++)
		{
			bins[j].Reset();
			float max = (j + 1) * interval + bounds.bmin[i];
			for (int k = counter; k < first + count; k++)
			{
				if (primitives[0][k]->GetCenter()[i] > max) {
					counter = k;
					for (int l = j; l < min(100, count); l++)
					{
						binCounts[l] = binCounts[j];
					}
					break;
				}
				bins[j].bmin3 = vec3::min(bins[j].bmin3, primitives[0][k]->GetBounds().bmin3);
				bins[j].bmax3 = vec3::max(bins[j].bmax3, primitives[0][k]->GetBounds().bmax3);
				binCounts[j]++;
			}
		}
		binCounts[min(100, count)] = count;

		// AABB's for bins
		aabb left = aabb(vec3(MAXFLOAT), vec3(MINFLOAT));
		aabb right = aabb(vec3(MAXFLOAT), vec3(MINFLOAT));
		right.Reset();

		// Calculate best bin split
		for (int j = min(100, count) - 1; j >= 0; j--)
		{
			left.Reset();
			for (int k = 0; k < j; k++)
			{
				left.bmin3 = vec3::min(bins[k].bmin3, left.bmin3);
				left.bmax3 = vec3::max(bins[k].bmax3, left.bmax3);

			}

			right.bmin3 = vec3::min(bins[j].bmin3, right.bmin3);
			right.bmax3 = vec3::max(bins[j].bmax3, right.bmax3);

			// Keep track
			float SAHScore = left.Area() * binCounts[j] + right.Area() * (binCounts[min(100, count)] - binCounts[j]);
			
			if (SAHScore < bestSplitScore) {
				bestSplitAxis = i;
				bestSplitPosition = primitives[0][first + binCounts[j] - 1]->GetCenter()[i];
				bestSplitScore = SAHScore;
			}
		}
	}

	// Check if splitting makes sense
	if (bestSplitScore >= currentSplitScore) return false;

	// Partition sort the primitives based on the found split
	uint leftIndex = first, rightIndex = first + count - 1;

	while (leftIndex < rightIndex) {
		Primitive *p = primitives[0][leftIndex];
		if (p->GetCenter()[bestSplitAxis] <= bestSplitPosition) {
			leftIndex++;
		} else {
			swap(primitives[0][leftIndex], primitives[0][rightIndex]);
			rightIndex--;
		}
	}

	Primitive *p = primitives[0][rightIndex];
	if (p->GetCenter()[bestSplitAxis] <= bestSplitPosition) {
		leftIndex++;
	}

	// Ugly check for resolving ties
	if (leftIndex - first == 0) {
		leftIndex++;
	}
	else if (first + count - leftIndex == 0) {
		leftIndex--;
	}

	// Set data for child nodes
	left->first = first;
	left->count = leftIndex - first;
	right->first = leftIndex;
	right->count = first + count - leftIndex;

	// Debug print
	//printf("Partitioning ready! Left first: %i, count: %i. Right first: %i, count: %i\n", left->first, left->count, right->first, right->count);
	return true;
}

float Tmpl8::Ray::AABBIntersect(aabb aabb)
{
	float tmin = (aabb.bmin3.x - origin.x) / direction.x;
	float tmax = (aabb.bmax3.x - origin.x) / direction.x;

	if (tmin > tmax) swap(tmin, tmax);

	float tymin = (aabb.bmin3.y - origin.y) / direction.y;
	float tymax = (aabb.bmax3.y - origin.y) / direction.y;

	if (tymin > tymax) swap(tymin, tymax);

	if ((tmin > tymax) || (tymin > tmax))
		return std::numeric_limits<float>::max();

	if (tymin > tmin)
		tmin = tymin;

	if (tymax < tmax)
		tmax = tymax;

	float tzmin = (aabb.bmin3.z - origin.z) / direction.z;
	float tzmax = (aabb.bmax3.z - origin.z) / direction.z;

	if (tzmin > tzmax) swap(tzmin, tzmax);

	if ((tmin > tzmax) || (tzmin > tmax))
		return std::numeric_limits<float>::max();

	if (tzmin > tmin)
		tmin = tzmin;

	if (tzmax < tmax)
		tmax = tzmax;

	return max(tmin, 0.0f);
}
