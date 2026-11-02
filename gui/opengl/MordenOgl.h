#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>

GLFWwindow * ModernOpenGLInit();
void SetupGeometries();
void RenderScene();
void SetupOpenGL();
void SetupSceneData();
void SetupCallbacks( GLFWwindow * window );
void KeyCallback( GLFWwindow * window, int key, int scancode, int action, int mods );
void WindowSizeCallback( GLFWwindow * window, int width, int height );
bool CheckOpenGLErrors();
void SetupShaders();
unsigned int CreateShader( const char * shaderSource, unsigned int shaderDefine );
unsigned int SetupShadeVertfrag( const char * vertexShaderSource, const char * fragmentShaderSource );
GLuint CheckCompilationShader( GLuint shader );
GLuint CheckLinkStatus( GLuint program );

class ModernOpenGL
{
public:
    ModernOpenGL();
    ~ModernOpenGL();
public:
    static unsigned int currProgram;
};