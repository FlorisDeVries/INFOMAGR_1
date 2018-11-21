#pragma once
#include <vector>
#include <limits>

#define EPSILON 0.001f
#define MAX_DEPTH 10
#define SCENE 1

#define ONRAILS false
#define LIGHTINTENSITY 10.0f

namespace Tmpl8
{
class Ray
{
  public:
	  Ray(vec3 origin, vec3 direction) : origin(origin), direction(direction) {};
	vec3 origin, direction;
};

class Primitive;

class Intersection
{
public:
	vec3 position, normal;
	float t = std::numeric_limits<float>::max();
	Primitive* primitive;
};

class Primitive
{
  public:
	vec3 color;
	float specularity, refractionIndex = 0.0f;
	virtual bool Intersect( Ray &ray , Intersection &intersection) = 0;
	virtual vec3 GetColor(vec3 pos) = 0;
protected:
	Primitive(vec3 color, float specularity, float refractionIndex) : color(color), specularity(specularity), refractionIndex(refractionIndex) {};
};

class Sphere : public Primitive
{
  public:
	  Sphere(vec3 pos, float r, vec3 color, float specularity, float refractionIndex) : position(pos), r2(r * r), Primitive(color, specularity, refractionIndex) {};
	bool Intersect( Ray &ray, Intersection &intersection) override;
	vec3 GetColor(vec3 pos) override;
  private:
	vec3 position;
	float r2;
};

class Plane : public Primitive
{
  public:
	  Plane(vec3 normal, float dist, vec3 color, float specularity, float refractionIndex) : normal(normal.normalized()), dist(dist), Primitive(color, specularity, refractionIndex) {};
	bool Intersect(Ray &ray, Intersection &intersection) override;
	vec3 GetColor(vec3 pos) override;

  private:
	vec3 normal;
	float dist;
};

class Camera
{
  public:
	Camera( vec3 pos, vec3 dir, float FOV );
	Ray GetRay(int x, int y);
	vec3 position, direction, screenTopLeft;
	float FOV;
	int screenWidth, screenHeight;
	void ResetBounds();

  private:
	vec3 screenCenter;
	vec3 ScreenCorner(int corner = 0);
	vec3 xinc, yinc;
};

class Light
{
  public:
	  Light(vec3 col, vec3 pos) : color(col), position(pos) {};
	vec3 position;
	vec3 color;
	virtual bool InLoS() = 0;
};

class PointLight : public Light
{
  public:
	  PointLight(vec3 col, vec3 pos) : Light(col, pos) {};
	inline bool InLoS() override { return true; }
};

class Game
{
  public:
	void SetTarget( Surface *surface ) { screen = surface; }
	void Init();
	void Shutdown();
	vec3 Trace(Ray ray, int recursionDepth);
	Ray Reflect(Ray &ray, Intersection &intersection);
	vec3 DirectIllumination(Ray &ray, Intersection &intersection);
	vec3 Refract(Ray &ray, Intersection &intersection, int recursionDepth);
	void Tick( float deltaTime );
	void MouseUp( int button )
	{ /* implement if you want to detect mouse button presses */
	}
	void MouseDown( int button )
	{ /* implement if you want to detect mouse button presses */
	}
	void MouseMove( int x, int y )
	{ /* implement if you want to detect mouse movement */
	}
	void KeyUp( int key )
	{ /* implement if you want to handle keys */
	}
	void KeyDown( int key )
	{ /* implement if you want to handle keys */
	}

	std::vector<Primitive *> primitives;
	std::vector<Light *> lights;
	Camera* cam;

  private:
	Surface *screen;
};

}; // namespace Tmpl8