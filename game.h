#pragma once
#include <atomic>
#include <limits>
#include <tuple>
#include <vector>

#pragma region Settings
#define EPSILON 0.001f
#define MAX_DEPTH 3
#define SCENE 1
#define MOVEMENTRATE 1.0f
#define SENSITIVITY 0.003f

#define ONRAILS false
#define LIGHTINTENSITY 10.0f
#define THREADS 8  // 8 threads might melt Windows
#define TILES 1024 // Keep as a power of an even number
#pragma endregion

#define MAXFLOAT std::numeric_limits<float>::max()
#define MINFLOAT std::numeric_limits<float>::min()

namespace Tmpl8
{
class Ray
{
  public:
	Ray( vec3 origin, vec3 direction ) : origin( origin ), direction( direction ){};
	vec3 origin, direction;
};

class Primitive;

class Intersection
{
  public:
	vec3 position, normal;
	float t = std::numeric_limits<float>::max();
	Primitive *primitive;
	bool inside = false;
};

class Primitive
{
  public:
	vec3 color, absorptionColor;
	float specularity, refractionIndex = 0.0f;
	virtual bool Intersect( Ray &ray, Intersection &intersection ) = 0;
	virtual vec3 GetColor( vec3 pos ) = 0;
	virtual aabb GetBounds() = 0;
	virtual vec3 GetCenter() = 0;

  protected:
	Surface *texture;
	Primitive( vec3 color, float specularity = 0, float refractionIndex = 0, vec3 absorptionColor = 1, Surface *texture = 0 ) : color( color ), specularity( specularity ), refractionIndex( refractionIndex ), absorptionColor( absorptionColor ), texture( texture ){};
};

class Sphere : public Primitive
{
  public:
	Sphere( vec3 pos, float r, vec3 color, float specularity = 0, float refractionIndex = 0, vec3 absorptionColor = 1, Surface *texture = 0 ) : position( pos ), r2( r * r ), Primitive( color, specularity, refractionIndex, absorptionColor, texture ){};
	bool Intersect( Ray &ray, Intersection &intersection ) override;
	vec3 GetColor( vec3 pos ) override;
	aabb GetBounds() override;
	vec3 GetCenter() override { return position; };

  private:
	vec3 position;
	float r2;
};

class Plane : public Primitive
{
  public:
	Plane( vec3 normal, float dist, vec3 color, float specularity = 0, float refractionIndex = 0, vec3 absorptionColor = 1, Surface *texture = 0 );
	bool Intersect( Ray &ray, Intersection &intersection ) override;
	vec3 GetColor( vec3 pos ) override;
	aabb GetBounds() override;
	vec3 GetCenter() override;

  private:
	vec3 normal, UAxis, VAxis;
	float dist;
};

class Triangle : public Primitive
{
  public:
	Triangle( vec3 v0, vec3 v1, vec3 v2, vec3 normal, vec3 color = 1, float specularity = 0, float refractionIndex = 0, vec3 absorptionColor = 1, Surface *texture = 0 ) : v0( v0 ), v1( v1 ), v2( v2 ), normal( normal ), Primitive( color, specularity, refractionIndex, absorptionColor, texture ){};
	bool Intersect( Ray &ray, Intersection &intersection ) override;
	vec3 GetColor( vec3 pos ) override;
	aabb GetBounds() override;
	vec3 GetCenter() override;

  private:
	vec3 normal, v0, v1, v2;
};

class Camera
{
  public:
	Camera( vec3 pos, vec3 dir, float FOV, float aspectRatio, int screenWidth, int screenHeight );
	Ray GetRay( int x, int y );
	std::tuple<int, std::vector<Ray>> GetNextRays();
	vec3 position, direction, screenTopLeft;
	float FOV, aspectRatio;
	int screenWidth, screenHeight;
	void ResetBounds();
	void ResetFOV();
	void ResetCounter();

  private:
	std::atomic<int> rayCounter = 0;
	float FOV_Distance;
	vec3 screenCenter;
	vec3 ScreenCorner( int corner = 0 );
	vec3 xinc, yinc;
};

class Light
{
  public:
	Light( vec3 col, vec3 pos ) : color( col ), position( pos ){};
	vec3 position;
	vec3 color;
	virtual bool InLoS() = 0;
};

class PointLight : public Light
{
  public:
	PointLight( vec3 col, vec3 pos ) : Light( col, pos ){};
	inline bool InLoS() override { return true; }
};

struct BVHNode
{
	aabb bounds;
	bool isLeaf;
	BVHNode *left, *right;
	int first, count;
	void Subdivide(std::vector<Primitive *> primitives);
	void Partition(std::vector<Primitive *> primitives);
	aabb CalculateBounds(std::vector<Primitive *> primitives);
};

class BVH
{
public:
	BVHNode root;
	BVHNode *pool;
	uint poolIdx;
	uint *indices;
	void ConstructBVH( std::vector<Primitive *> primitives );
};

class Game
{
  public:
	  std::vector<Primitive *> primitives;
	  std::vector<Primitive *> nonBVHprimitives;
	std::vector<Light *> lights;
	Camera *cam;
	bool mouseTurning;

	void SetTarget( Surface *surface ) { screen = surface; }
	void Init();
	void Shutdown();
	vec3 Trace( Ray ray, int recursionDepth, Intersection &intersection = Intersection() );
	Ray Reflect( Ray &ray, Intersection &intersection );
	vec3 DirectIllumination( Ray &ray, Intersection &intersection );
	vec3 Refract( Ray &ray, Intersection &intersection, int recursionDepth );
	void Tick( float deltaTime );
	void ThreadedRays( int i );
	void HandleInput();
	bool ReadObj( const char *path, std::vector<Primitive *> &primitives, vec3 color, vec3 position, float specularity = 0 );
	void MouseUp( int button )
	{ /* implement if you want to detect mouse button presses */
		if ( button == SDL_BUTTON_RIGHT )
			mouseTurning = false;
	}
	void MouseDown( int button )
	{ /* implement if you want to detect mouse button presses */
		if ( button == SDL_BUTTON_RIGHT )
			mouseTurning = true;
	}
	void MouseMove( int x, int y )
	{ /* implement if you want to detect mouse movement */
		if ( mouseTurning )
		{
			//printf("Turning camera\n");
			vec3 left = cross( cam->direction, vec3( 0, 1, 0 ) ).normalized();
			vec3 up = cross( cam->direction, left ).normalized();
			int xdif = -x;
			int ydif = -y;

			vec3 target = cam->position + cam->direction;
			target += left * xdif * SENSITIVITY + up * -ydif * SENSITIVITY;
			cam->direction = ( target - cam->position ).normalized();
			cam->ResetBounds();
		}
	}

	void KeyUp( int key );
	void KeyDown( int key );

  private:
	BVH bvh;
	Surface *screen;
	bool isWDown = false, isADown = false, isSDown = false, isDDown = false, isRDown = false, isFDown = false, isYDown = false, isHDown = false;
};

}; // namespace Tmpl8