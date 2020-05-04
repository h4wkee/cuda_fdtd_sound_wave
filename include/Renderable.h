#pragma once

#include <vector>
#include <glad/glad.h>

#include <glm/glm.hpp>

class ResourceSet;

class Renderable
{
public:
	Renderable();
	Renderable(std::vector<glm::vec3> & positions, GLenum type = GL_POINTS, float pointSize = 1.f);
	Renderable(Renderable&& renderable) noexcept;
	~Renderable();

	void draw();
	void init(std::vector<glm::vec3> & positions, GLenum type = GL_POINTS, float pointSize = 1.f);

private:
	Renderable(const Renderable& renderable) = delete;
	Renderable& operator=(const Renderable& renderable) = delete;

	GLuint _vao = 0;
	GLuint _vbo = 0;

	GLenum _type;
	unsigned int _positionsCount = 0;
};
