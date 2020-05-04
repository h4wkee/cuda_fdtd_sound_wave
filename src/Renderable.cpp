#include <Renderable.h>

#include <Window.h>

Renderable::Renderable()
{

}

Renderable::Renderable(std::vector<glm::vec3> & positions, GLenum type, float pointSize)
{
	init(positions, type);
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
}

void Renderable::draw()
{
	glBindVertexArray(_vao);
	glDrawArrays(_type, 0, _positionsCount);
	glBindVertexArray(0);
}

void Renderable::init(std::vector<glm::vec3> & positions, GLenum type, float pointSize)
{
	glGenBuffers(1, &_vbo);
	_positionsCount = positions.size();
	_type = type;

	glBindBuffer(GL_ARRAY_BUFFER, _vbo);
	glBufferData(GL_ARRAY_BUFFER, _positionsCount * sizeof(glm::vec3), &positions[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glPointSize(pointSize);

	if (_vao == 0)
	{
		glGenVertexArrays(1, &_vao);
	}

	glBindVertexArray(_vao);

	glBindBuffer(GL_ARRAY_BUFFER, _vbo);

	//position
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
	glEnableVertexAttribArray(0);

	glBindVertexArray(0);
}