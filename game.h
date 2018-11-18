#pragma once
#include <array>

namespace Tmpl8
{
class Ray
{
  public:
	Ray(vec3 origin, vec3 direction);
	vec3 origin, direction;
	float t;
	void SetT(float haha);
};

class Primitive
{
  public:
	virtual void Intersect( Ray &ray ) = 0;
};

class Sphere : public Primitive
{
  public:
	Sphere( vec3 pos, float r );
	void Intersect( Ray &ray ) override;

  private:
	vec3 position;
	float r2;
};

class Plane : public Primitive
{
  public:
	Plane(vec3 normal, float dist);
	void Intersect(Ray &ray) override;

  private:
	vec3 normal;
	float dist;
};

class Camera
{
  public:
	Camera( vec3 pos, vec3 dir, float FOV, int screenWidth, int screenHeight );
	Ray GetRay(int x, int y);
	vec3 position, direction, screenTopLeft;
	float FOV;
	int screenWidth, screenHeight;

  private:
	vec3 screenCenter;
	vec3 ScreenCorner(int corner = 0);
	vec3 xinc, yinc;
};

class Light
{
  public:
	Light( vec3 col, vec3 pos );
	vec3 color, position;
	virtual bool InLoS() = 0;
};

class PointLight : public Light
{
  public:
	PointLight( vec3 col, vec3 pos );
	inline bool InLoS() override { return true; }
};

class Game
{
  public:
	void SetTarget( Surface *surface ) { screen = surface; }
	void Init();
	void Shutdown();
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

	std::array<Primitive *, 3> primitives;
	std::array<Light *, 1> lights;
	Camera* cam;

  private:
	Surface *screen;
};

}; // namespace Tmpl8