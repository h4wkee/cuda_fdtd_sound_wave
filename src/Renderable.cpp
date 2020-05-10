#include <Renderable.h>

#include <Window.h>

Renderable::Renderable()
{

}

Renderable::Renderable(std::vector<Vertex> & vertices, GLenum type, float pointSize)
{
	init(vertices, type, pointSize);
}

Renderable::Renderable(Renderable && renderable) noexcept
{
	_vao = renderable._vao;
	_vbo = renderable._vbo;
	renderable._vao = 0;
	renderable._vbo = 0;
}

Renderable::~Renderable()
{
	if (_vbo)
	{
		glDeleteBuffers(1, &_vbo);
	}
	if (_vao)
	{
		glDeleteVertexArrays(1, &_vao);
	}
//	if(_vertexPointer)
//	{
//		cudaGLUnmapbufferObject(_vbo);
//	}
}

void Renderable::draw()
{
	glBindVertexArray(_vao);
	glDrawArrays(_type, 0, _vertexCount);
	glBindVertexArray(0);
}

void Renderable::init(std::vector<Vertex> & vertices, GLenum type, float pointSize)
{
	glGenBuffers(1, &_vbo);
	_vertexCount = vertices.size();
	_type = type;

	glBindBuffer(GL_ARRAY_BUFFER, _vbo);
	glBufferData(GL_ARRAY_BUFFER, _vertexCount * sizeof(Vertex), &vertices[0], GL_STATIC_DRAW); // CUDA: GL_DYNAMIC_COPY?
	//cudaGLRegisterBufferObject(_vbo);
	//cudaGLMapBufferObject(&_vertexPointer, _vbo);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	//checkCudaErrors(cudaGraphicsGLRegisterBuffer(&_vertexPointer, _vbo, cudaGraphicsMapFlagsNone));

	glPointSize(pointSize);

	glGenVertexArrays(1, &_vao);

	glBindVertexArray(_vao);

	glBindBuffer(GL_ARRAY_BUFFER, _vbo);

	//position
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
	glEnableVertexAttribArray(0);

	//color
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
	glEnableVertexAttribArray(1);

	glBindVertexArray(0);
}