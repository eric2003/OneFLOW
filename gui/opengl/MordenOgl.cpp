#include "MordenOgl.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
using namespace std;

int CurrentMode = 1;	// Controls what is drawn.
const int NumObjects = 3;
const int iPoints = 0;
const int iLines = 1;
const int iTriangles = 2;

unsigned int myVBO[NumObjects];  // Vertex Buffer Object - holds an array of data
unsigned int myVAO[NumObjects];  // Vertex Array Object - holds info about how the vertex data is formatted

// We create one shader program: it consists of a vertex shader and a fragment shader
unsigned int shaderProgram1;
const unsigned int vertPos_loc = 0;   // Corresponds to "location = 0" in the verter shader definition
const unsigned int vertColor_loc = 1; // Corresponds to "location = 1" in the verter shader definition

unsigned int ModernOpenGL::currProgram = 0;

ModernOpenGL::ModernOpenGL()
{
    ;
}

ModernOpenGL::~ModernOpenGL()
{
    ;
}

void SetupGeometries()
{
    // Allocate Vertex Array Objects (VAOs) and Vertex Buffer Objects (VBOs).
    glGenVertexArrays( NumObjects, &myVAO[0] );
    glGenBuffers( NumObjects, &myVBO[0] );

    float threeVerts[] = {
        -0.6f, -0.3f,  // First point
        0.6f, -0.3f,   // Second point
        0.6f, 0.3f	   // Third point
    };

    glBindVertexArray( myVAO[iPoints] );
    glBindBuffer( GL_ARRAY_BUFFER, myVBO[iPoints] );
    glBufferData( GL_ARRAY_BUFFER, sizeof(threeVerts), threeVerts, GL_STATIC_DRAW );

    glVertexAttribPointer( vertPos_loc, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);	// Info about where positions are in the VBO
    glEnableVertexAttribArray( vertPos_loc );									// Enable the stored vertex positions

    float sixVertsForLines[] = {
        -0.2f, -0.8f,
        0.8f, 0.2f,
        0.6f, 0.8f,
        -0.5f, 0.6f,
        -0.4f, -0.0f,
        -0.7f, -0.3f,
    };
    // Do the same as above, now for the vertices that specify lines.
    glBindVertexArray( myVAO[iLines] );
    glBindBuffer( GL_ARRAY_BUFFER, myVBO[iLines] );
    glBufferData( GL_ARRAY_BUFFER, sizeof( sixVertsForLines ), sixVertsForLines, GL_STATIC_DRAW);
    glVertexAttribPointer( vertPos_loc, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glEnableVertexAttribArray( vertPos_loc );

    float trianglesVerts[] = {
        // x,y,z coordinates	// R,G,B colors
        0.7f, -0.42f, 0.0f,		1.0f, 0.8f, 0.8f, // First triangle
        0.7f, -0.18f, 0.0f,		1.0f, 0.8f, 0.8f,
        -0.7f, -0.3f, 0.5f,		1.0f, 0.0f, 0.0f,

        -0.25f, 0.7f, 0.0f,		0.8f, 1.0f, 0.8f, // Second triangle
        -0.40f, 0.55f, 0.0f,	0.8f, 1.0f, 0.8f,
        0.5f, -0.6f, 0.5f,		0.0f, 1.0f, 0.0f,

        -0.57f, -0.53f, 0.0f,	0.8f,  0.8f, 1.0f,	// Third triangle
        -0.43f, -0.67f, 0.0f,	0.8f,  0.8f, 1.0f,
        0.32f, 0.62f, 0.5f,		0.0f,  0.0f, 1.0f,
    };

    glBindVertexArray( myVAO[iTriangles] );
    glBindBuffer( GL_ARRAY_BUFFER, myVBO[ iTriangles ] );
    glBufferData( GL_ARRAY_BUFFER, sizeof( trianglesVerts ), trianglesVerts, GL_STATIC_DRAW );
    glVertexAttribPointer( vertPos_loc, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray( vertPos_loc );
    glVertexAttribPointer( vertColor_loc, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray( vertColor_loc );

    glBindBuffer( GL_ARRAY_BUFFER, 0 );
    glBindVertexArray( 0 );

    CheckOpenGLErrors();
}

void RenderScene()
{
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glUseProgram( ModernOpenGL::currProgram );

    switch ( CurrentMode )
    {
    case 0:
        glBindVertexArray( myVAO[iTriangles] );
        glDrawArrays( GL_TRIANGLES, 0, 9 );
        break;
    case 1:
        glBindVertexArray( myVAO[iLines] );
        glVertexAttrib3f( vertColor_loc, 0.5f, 1.0f, 0.2f );		// A greenish color (R, G, B values).
        glDrawArrays( GL_LINES, 0, 6 );
        break;
    case 2:
        glBindVertexArray( myVAO[iLines] );
        glVertexAttrib3f( vertColor_loc, 1.0f, 0.2f, 1.0f );		// Magenta color (R, G, B values).
        glDrawArrays( GL_LINE_STRIP, 0, 6 );
        break;
    case 3:
        glBindVertexArray( myVAO[iLines] );
        glVertexAttrib3f( vertColor_loc, 1.0f, 1.0f, 0.2f );		// A yellow-ish color (R, G, B values).
        glDrawArrays( GL_LINE_LOOP, 0, 6 );
        break;
    case 4:
        glBindVertexArray( myVAO[iLines] );
        glVertexAttrib3f( vertColor_loc, 1.0f, 1.0f, 1.0f );		// A white color (R, G, B values).
        glDrawArrays( GL_POINTS, 0, 6 );
        break;
    case 5:
        // Draw three points
        glBindVertexArray( myVAO[iPoints] );
        glVertexAttrib3f( vertColor_loc, 1.0f, 0.5f, 0.2f );		// An orange-red color (R, G, B values).
        glDrawArrays( GL_POINTS, 0, 3 );
        break;
    }

    glBindVertexArray( 0 );
    CheckOpenGLErrors();
}

void KeyCallback( GLFWwindow * window, int key, int scancode, int action, int mods )
{
    if ( action == GLFW_RELEASE )
    {
        return;
    }
    if ( key == GLFW_KEY_ESCAPE || key == GLFW_KEY_X )
    {
        glfwSetWindowShouldClose(window, true);
    }
    else if ( key == GLFW_KEY_SPACE )
    {
        CurrentMode = ( CurrentMode + 1 ) % 6;	// Takes on values from 0 to 5
    }
}

char errNames[9][36] =
{
    "Unknown OpenGL error",
    "GL_INVALID_ENUM", "GL_INVALID_VALUE", "GL_INVALID_OPERATION",
    "GL_INVALID_FRAMEBUFFER_OPERATION", "GL_OUT_OF_MEMORY",
    "GL_STACK_UNDERFLOW", "GL_STACK_OVERFLOW", "GL_CONTEXT_LOST"
};

bool CheckOpenGLErrors()
{
    int numErrors = 0;
    GLenum err;
    while ((err = glGetError()) != GL_NO_ERROR) {
        numErrors++;
        int errNum = 0;
        switch (err) {
        case GL_INVALID_ENUM:
            errNum = 1;
            break;
        case GL_INVALID_VALUE:
            errNum = 2;
            break;
        case GL_INVALID_OPERATION:
            errNum = 3;
            break;
        case GL_INVALID_FRAMEBUFFER_OPERATION:
            errNum = 4;
            break;
        case GL_OUT_OF_MEMORY:
            errNum = 5;
            break;
        case GL_STACK_UNDERFLOW:
            errNum = 6;
            break;
        case GL_STACK_OVERFLOW:
            errNum = 7;
            break;
        case GL_CONTEXT_LOST:
            errNum = 8;
            break;
        }
        cout << "OpenGL ERROR: " << errNames[errNum] << "\n";
    }
    return ( numErrors != 0 );
}

void SetupOpenGL()
{
    glEnable( GL_DEPTH_TEST );	// Enable depth buffering
    glDepthFunc( GL_LEQUAL );	// Useful for multipass shaders
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    glPointSize( 8 );
    glLineWidth( 5 );

}

void WindowSizeCallback( GLFWwindow * window, int width, int height )
{
    glViewport(0, 0, width, height);// Draw into entire window
}

void SetupCallbacks( GLFWwindow * window )
{
    // Set callback function for resizing the window
    glfwSetFramebufferSizeCallback( window, WindowSizeCallback );

    // Set callback for key up/down/repeat events
    glfwSetKeyCallback( window, KeyCallback );
}

void SetupSceneData()
{
    SetupGeometries();

    SetupShaders();

    CheckOpenGLErrors();
}

// Sets the position and color of a vertex.
//   This implementations just copies the position with no transformations.
//   It copies the color to "theColor" so that the fragment shader can access it.
const char *vertexShader_PosColorOnly =
"#version 330 core\n"
"layout (location = 0) in vec3 vertPos;	   // Position in attribute location 0\n"
"layout (location = 1) in vec3 vertColor;  // Color in attribute location 1\n"
"out vec3 theColor;					       // Output a color to the fragment shader\n"
"void main()\n"
"{\n"
"   gl_Position = vec4(vertPos.x, vertPos.y, vertPos.z, 1.0);\n"
"   theColor = vertColor;\n"
"}\0";

// Set a general color using a fragment shader. (A "fragment" is a "pixel".)
//    The color value is passed in, obtained from the colors on the vertices.
//    Color values range from 0.0 to 1.0.
//    First three values are Red/Green/Blue (RGB).
//    Fourth color value (alpha) is 1.0, meaning there is no transparency.
const char *fragmentShader_ColorOnly =
"#version 330 core\n"
"in vec3 theColor;		// Color value came from the vertex shader (smoothed) \n"
"out vec4 FragColor;	// Color that will be used for the fragment\n"
"void main()\n"
"{\n"
"   FragColor = vec4(theColor, 1.0f);   // Add alpha value of 1.0.\n"
"}\n\0";

/*
 * Build and compile our shader programs
 */ 
void SetupShaders()
{
    ModernOpenGL::currProgram = SetupShadeVertfrag( vertexShader_PosColorOnly, fragmentShader_ColorOnly );
}

unsigned int CreateShader( const char * shaderSource, unsigned int shaderDefine )
{
    unsigned int outShader = glCreateShader( shaderDefine );
    glShaderSource( outShader, 1, &shaderSource, NULL );
    glCompileShader( outShader );
    CheckCompilationShader( outShader );
    return outShader;
}

unsigned int SetupShadeVertfrag( const char* vertexShaderSource, const char* fragmentShaderSource )
{
    // vertex shader
    unsigned int vertexShader = CreateShader( vertexShaderSource, GL_VERTEX_SHADER );
    // fragment shader
    unsigned int fragmentShader = CreateShader( fragmentShaderSource, GL_FRAGMENT_SHADER );

    // link shaders
    unsigned int shaderProgram = glCreateProgram();
    glAttachShader( shaderProgram, vertexShader );
    glAttachShader( shaderProgram, fragmentShader );
    glLinkProgram( shaderProgram );
    CheckLinkStatus( shaderProgram );

    // Deallocate shaders since we do not need to use these for other shader programs.
    glDeleteShader( vertexShader );
    glDeleteShader( fragmentShader );

    return shaderProgram;		// Return the compiled shaders as a single shader program.
}

GLuint CheckCompilationShader( GLuint shader )
{
    if ( !glIsShader( shader ) )
    {
        cout << "ERROR: Not a shader! Possibly an allocation error.\n";
        return 0;
    }

    int success;
    glGetShaderiv( shader, GL_COMPILE_STATUS, &success );
    if ( success )
    {
        return shader;	// Compilation was successful
    }

    // Compilation failed
    int infoLogLength;
    glGetShaderiv( shader, GL_INFO_LOG_LENGTH, &infoLogLength );
    char* infoLog = new char[ infoLogLength ];
    glGetShaderInfoLog( shader, infoLogLength, NULL, infoLog );
    cout << "ERROR::Shader compilation failed!\n" << infoLog << "\n";
    delete[] infoLog;
    return 0;
}

GLuint CheckLinkStatus( GLuint program )
{
    if ( !glIsProgram( program ) ) 
    {
        cout << "ERROR: Not a shader program! Possibly an allocation error.\n";
        return 0;
    }

    int success;
    glGetProgramiv( program, GL_LINK_STATUS, &success );
    if ( success )
    {
        return program;		// Linkage was successful
    }

    int infoLogLength;
    glGetProgramiv( program, GL_INFO_LOG_LENGTH, &infoLogLength );
    char* infoLog = new char[ infoLogLength ];
    glGetProgramInfoLog( program, infoLogLength, NULL, infoLog );
    cout << "ERROR::Shader program link failed!\n" << infoLog << "\n";
    return 0;
}


GLFWwindow * ModernOpenGLInit()
{
    glfwInit();

    const int initWidth = 800;
    const int initHeight = 600;
    GLFWwindow * window = glfwCreateWindow( initWidth, initHeight, "OneFLOW CFD Modern OpenGL Test", NULL, NULL );

    if ( window == NULL )
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return 0;
    }

    glfwMakeContextCurrent( window );
    glfwSetFramebufferSizeCallback( window, WindowSizeCallback );

    if ( ! gladLoadGLLoader((GLADloadproc)glfwGetProcAddress ) )
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return 0;
    }

    SetupCallbacks( window );

    SetupOpenGL();
    SetupSceneData();

    return window;
}