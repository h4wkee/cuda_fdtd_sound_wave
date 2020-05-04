#pragma once

#include <map>
#include <vector>
#include <string>

#include <glad/glad.h>

class Shader
{
public:
	typedef float DefineValue;

	enum class Uniform
	{
		texture_diffuse0 = 0,
		texture_specular0,
		texture_normal0,
	};
	static constexpr unsigned uniformCount = static_cast<unsigned>(Uniform::texture_normal0) + 1;

	Shader(const std::string & name, bool geometryShader = false);
	~Shader();

	bool reload();
	bool load(const std::string & name);
	bool load(const std::string & vertexFileName, const std::string & fragmentFileName, const std::string & geometryFileName = "");
	void use();
	unsigned int getId() { return _program; }

	void define(const std::string & declaration, DefineValue value);
	void setDefines(const std::vector<std::pair<std::string, DefineValue>> & defines) { _defines = defines; }
	void clearDefines() { _defines.clear(); }

	static Shader * getActive() { return _active; };
	static DefineValue getGlobalDefine(const std::string & declaration);
	static void setGlobalDefines(const std::vector<std::pair<std::string, DefineValue>> & globalDefines) { _globalDefines = globalDefines; }
	static void clearGlobalDefines() { _globalDefines.clear(); }

	static std::vector<std::pair<std::string, DefineValue>> getGlobalDefines() { return _globalDefines; }
	static void globalDefine(const std::string & declaration, DefineValue value);
	static std::string getDefineType(const std::string & declaration);
	static std::string getDefineName(std::string declaration);

	void setUniform1(const std::string & name, int value);
	void setUniform1(const std::string & name, float value);
	void setUniformArray(const std::string & name, float * values, unsigned int count);
	void setUniform3(const std::string & name, float x, float y, float z);
	void setUniform4m(const std::string & name, float * data, int count = 1);

	void setUniform1(const Uniform uniform, int value);
	void setUniform1(const Uniform uniform, float value);
	void setUniformArray(const Uniform uniform, float * values, unsigned int count);
	void setUniform3(const Uniform uniform, float x, float y, float z);
	void setUniform4m(const Uniform uniform, float * data, int count = 1);

private:
	static Shader * _active;

	static std::vector<std::pair<std::string, DefineValue>> _globalDefines;
	std::vector<std::pair<std::string, DefineValue>> _defines;
	//std::map<std::string, int> _uniforms;
	int _uniforms[uniformCount];
	std::string _name;
	bool _geometryShader;
	GLuint _program = 0;

	bool _readFile(const std::string& path, std::string * shaderSource);
	bool _compile(const std::string& source, GLuint shader);
	bool _link(GLuint vertexShader, GLuint fragmentShader, GLuint geometryShader = 0);
	int _getUniform(const std::string & name);
	void _free();

	static std::string _getDefineCode(const std::string & declaration, DefineValue value);
};
