#pragma once

#include <map>
#include <vector>
#include <string>

#include <glad/glad.h>

class Shader
{
public:
	Shader(const std::string & name);
	~Shader();

	bool reload();
	bool load(const std::string & name);
	bool load(const std::string & vertexFileName, const std::string & fragmentFileName);
	void use();

	void setUniform4m(const std::string & name, float * data, int count = 1);

private:
	static Shader * _active;

	std::string _name;
	GLuint _program = 0;

	bool _readFile(const std::string& path, std::string * shaderSource);
	bool _compile(const std::string& source, GLuint shader);
	bool _link(GLuint vertexShader, GLuint fragmentShader, GLuint geometryShader = 0);
	int _getUniform(const std::string & name);
	void _free();
};
