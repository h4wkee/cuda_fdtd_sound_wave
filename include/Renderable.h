#pragma once

#include <vector>

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <Vertex.h>

class Renderable
{
public:
	Renderable();
	Renderable(std::vector<Vertex> & vertices, GLenum type = GL_POINTS, float pointSize = 1.f);
	Renderable(Renderable&& renderable) noexcept;
	~Renderable();

	void draw();
	void init(std::vector<Vertex> & vertices, GLenum type = GL_POINTS, float pointSize = 1.f);

	GLuint * getVBO() { return & _vbo; }

private:
	Renderable(const Renderable& renderable) = delete;
	Renderable& operator=(const Renderable& renderable) = delete;

	GLuint _vao = 0;
	GLuint _vbo = 0;

	GLenum _type;
	unsigned int _vertexCount = 0;
};
