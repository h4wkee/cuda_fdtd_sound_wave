#include <Shader.h>

#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>

static const std::string shadersDirectory = "resources/shaders/";

Shader::Shader(const std::string & name)
{
	load(name);
}

Shader::~Shader()
{
	_free();
}

Shader * Shader::_active = nullptr;

bool Shader::reload()
{
	_free();
	bool loaded = load(_name + ".vert", _name + ".frag");

	if(loaded)
	{
		std::cout << "Shader '" << _name << "' loaded successfully" << std::endl;
		return true;
	}
	else
	{
		std::cout << "Shader '" << _name << "' cannot be loaded properly" << std::endl;
	}
	return false;
}

bool Shader::load(const std::string & name)
{
	_name = name;
	return reload();
}

bool Shader::load(const std::string & vertexFileName, const std::string & fragmentFileName)
{
	//_uniforms.clear(); // clear uniform map

	std::string vertexPath = shadersDirectory + vertexFileName;
	std::string fragmentPath = shadersDirectory + fragmentFileName;

	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
	GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	GLuint geometryShader = 0;

	if(!_compile(vertexPath, vertexShader) || !_compile(fragmentPath, fragmentShader))
	{
		glDeleteShader(vertexShader);
		glDeleteShader(fragmentShader);

		return false;
	}

	_program = glCreateProgram();

	bool linkResult = _link(vertexShader, fragmentShader, geometryShader);

	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);

	if (linkResult)
	{
		Shader * active = _active;
		use();
		if (active)
		{
			active->use();
		}
	}

	return linkResult;
}

bool Shader::_readFile(const std::string& path, std::string * shaderSource)
{
	std::string & source = *shaderSource;
	if (FILE* file = fopen(path.c_str(), "r"))
	{
		fseek(file, 0, SEEK_END);
		unsigned int fileLength = ftell(file);
		source.resize(fileLength + 1);
		rewind(file);
		int charactersRead = fread(&source[0], sizeof(char), fileLength, file);
		source[charactersRead] = '\0';
		fclose(file);
		return true;
	}
	else
	{
		std::cout << "Cannot read '" << path << "' file" << std::endl;
		return false;
	}
}

bool Shader::_compile(const std::string& path, GLuint shader)
{
	std::vector<std::string> includes;
	std::vector<const char*> sourcePointers;
	std::vector<std::string> sources = { "" };
	std::string & headerSource = sources[0];
	std::string shaderSource = "";

	if(!_readFile(path, &shaderSource))
	{
		return false;
	}

	// bake header - GLSL version and run-time specified defines
	std::string versionString = "#version " + std::to_string(GLVersion.major) + std::to_string(GLVersion.minor) + "0" + "\n";
	headerSource += versionString;

	// parse includes
	{
		std::istringstream sourceReader(shaderSource);
		const std::string includePrefix = "#include <";

		std::string line;
		while (std::getline(sourceReader, line))
		{
			// end parsing includes if line does not start properly
			if (line.rfind(includePrefix, 0) == std::string::npos)
			{
				break;
			}

			// get file name to include
			std::string include = line.substr(includePrefix.length(), line.length() - includePrefix.length() - 1);
			includes.push_back(include);
		}
		// erase all read includes from shader source
		shaderSource.erase(0, size_t(sourceReader.tellg()) - line.length() - 1);
	}

	// read included files
	for (auto & include : includes)
	{
		sources.push_back("");
		if (!_readFile(shadersDirectory + "include/" + include, &sources.back()))
		{
			return false;
		}
	}

	// append shaderSource at the end of sources
	sources.push_back(std::move(shaderSource));

	// store pointers of all sources for OpenGL
	for (auto & source : sources)
	{
		sourcePointers.push_back(source.c_str());
	}

	glShaderSource(shader, GLsizei(sourcePointers.size()), &sourcePointers[0], nullptr);
	glCompileShader(shader);

	GLint status;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
	if (status != GL_TRUE)
	{
		int length = 0;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
		std::string errors;
		errors.resize(length);
		glGetShaderInfoLog(shader, length, nullptr, &errors[0]);
		std::cout << "Compilation of '" << path << "' shader file failed:" << "\n" << errors << std::endl;
		return false;
	}
	return true;
}

bool Shader::_link(GLuint vertexShader, GLuint fragmentShader, GLuint geometryShader)
{
	glAttachShader(_program, vertexShader);
	glAttachShader(_program, fragmentShader);

	glLinkProgram(_program);

	glDetachShader(_program, vertexShader);
	glDetachShader(_program, fragmentShader);

	GLint status;
	glGetProgramiv(_program, GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		int length = 0;
		glGetProgramiv(_program, GL_INFO_LOG_LENGTH, &length);
		std::string errors;
		errors.resize(length);
		glGetProgramInfoLog(_program, length, nullptr, &errors[0]);
		std::cout << "Linking of '" << _name << "' shader failed:" << "\n" << errors << std::endl;
		return false;
	}

	return true;
}

int Shader::_getUniform(const std::string & name)
{
	//auto it = _uniforms.find(name);
	//if (it != _uniforms.end())
	//{
	//	return it->second;
	//}
	int location = glGetUniformLocation(_program, name.c_str());
	//_uniforms[name] = location;
	return location;
}

void Shader::_free()
{
	if(_program != 0)
	{
		GLint program;
		glGetIntegerv(GL_CURRENT_PROGRAM, &program);
		if(_program == program)
		{
			glUseProgram(0);
		}
		glDeleteProgram(_program);
		_program = 0;
	}
}

void Shader::use()
{
	glUseProgram(_program);
	_active = this;
}

void Shader::setUniform4m(const std::string & name, float * data, int count)
{
	int location = _getUniform(name);
	if (location != -1)
	{
		glUniformMatrix4fv(location, count, GL_FALSE, data);
	}
}
