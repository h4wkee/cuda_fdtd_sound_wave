#pragma once

#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
//#include <helper_cuda.h>
//#include <helper_cuda_gl.h>

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

	struct cudaGraphicsResource * getVertexBufferPointer() { return _vertexPointer; }

private:
	Renderable(const Renderable& renderable) = delete;
	Renderable& operator=(const Renderable& renderable) = delete;

	GLuint _vao = 0;
	GLuint _vbo = 0;
	//void * _vertexPointer = nullptr;
	struct cudaGraphicsResource * _vertexPointer;

	GLenum _type;
	unsigned int _vertexCount = 0;
};
