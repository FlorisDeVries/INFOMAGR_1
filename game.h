#pragma once
#include <limits>
#include <vector>

#define EPSILON 0.001f
#define MAX_DEPTH 3
#define SCENE 3
#define MOVEMENTRATE 1.0f
#define SENSITIVITY 0.003f

#define ONRAILS false
#define LIGHTINTENSITY 10.0f
#define THREADS 64				//Keep as a power of 2

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

  protected:
	Surface* texture;
	Primitive( vec3 color, float specularity, float refractionIndex, vec3 absorptionColor = 1, Surface *texture = 0) : color( color ), specularity( specularity ), refractionIndex( refractionIndex ), absorptionColor( absorptionColor ), texture(texture){};
};

class Sphere : public Primitive
{
  public:
	Sphere( vec3 pos, float r, vec3 color, float specularity, float refractionIndex, vec3 absorptionColor = 1, Surface* texture = 0 ) : position( pos ), r2( r * r ), Primitive( color, specularity, refractionIndex, absorptionColor, texture ){};
	bool Intersect( Ray &ray, Intersection &intersection ) override;
	vec3 GetColor( vec3 pos ) override;

  private:
	vec3 position;
	float r2;
};

class Plane : public Primitive
{
  public:
	  Plane(vec3 normal, float dist, vec3 color, float specularity, float refractionIndex, vec3 absorptionColor = 1, Surface* texture = 0);
	bool Intersect( Ray &ray, Intersection &intersection ) override;
	vec3 GetColor( vec3 pos ) override;

  private:
	vec3 normal, UAxis, VAxis;
	float dist;
};

class Camera
{
  public:
	Camera( vec3 pos, vec3 dir, float FOV, float aspectRatio );
	Ray GetRay(int x, int y);
	vec3 position, direction, screenTopLeft;
	float FOV, aspectRatio;
	int screenWidth, screenHeight;
	void ResetBounds();
	void ResetFOV();

  private:
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

class Game
{
  public:

	  std::vector<Primitive *> primitives;
	  std::vector<Light *> lights;
	  Camera* cam;
	  bool mouseTurning;

	void SetTarget( Surface *surface ) { screen = surface; }
	void Init();
	void Shutdown();
	vec3 Trace( Ray ray, int recursionDepth, Intersection &intersection = Intersection() );
	Ray Reflect( Ray &ray, Intersection &intersection );
	vec3 DirectIllumination( Ray &ray, Intersection &intersection );
	vec3 Refract( Ray &ray, Intersection &intersection, int recursionDepth );
	void Tick( float deltaTime );
	void ThreadedRays(int i);
	void MouseUp( int button )
	{ /* implement if you want to detect mouse button presses */
		if (button == SDL_BUTTON_RIGHT)
			mouseTurning = false;
	}
	void MouseDown( int button )
	{ /* implement if you want to detect mouse button presses */
		if (button == SDL_BUTTON_RIGHT)
			mouseTurning = true;
	}
	void MouseMove( int x, int y )
	{ /* implement if you want to detect mouse movement */
		if (mouseTurning) {
			//printf("Turning camera\n");
			vec3 left = cross(cam->direction, vec3(0, 1, 0)).normalized();
			vec3 up = cross(cam->direction, left).normalized();
			int xdif = -x;
			int ydif = -y;

			vec3 target = cam->position + cam->direction;
			target += left * xdif * SENSITIVITY + up * -ydif * SENSITIVITY;
			cam->direction = (target - cam->position).normalized();
			cam->ResetBounds();
		}
	}
	void KeyUp( int key )
	{ /* implement if you want to handle keys */
	}
	void KeyDown( int key )
	{ 
		vec3 left = cross(cam->direction, vec3(0, 1, 0)).normalized();
		vec3 up = cross(cam->direction, left).normalized();
			//Forward, left, backward, right, up, down
		
		if(key == SDL_SCANCODE_W) cam->position += cam->direction * MOVEMENTRATE;
		if(key == SDL_SCANCODE_A) cam->position += left * MOVEMENTRATE;
		if(key == SDL_SCANCODE_S) cam->position += -cam->direction * MOVEMENTRATE;
		if(key == SDL_SCANCODE_D) cam->position += -left * MOVEMENTRATE;
		if(key == SDL_SCANCODE_R) cam->position += up * MOVEMENTRATE;
		if(key == SDL_SCANCODE_F) cam->position += -up * MOVEMENTRATE;
			// Increase/decrease FOV
		if (key == SDL_SCANCODE_Y && cam->FOV > 2.2f) cam->FOV -= 0.1f;
		if (key == SDL_SCANCODE_H) cam->FOV += 0.1f;

		cam->ResetFOV();
		cam->ResetBounds();
	}

  private:
	Surface *screen;
};

}; // namespace Tmpl8