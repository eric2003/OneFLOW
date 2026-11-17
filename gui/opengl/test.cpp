#include "MordenOgl.h"
#include <iostream>
using namespace std;

int main( int argc, char ** argv )
{
    GLFWwindow * window = ModernOpenGLInit();

    // Loop while program is not terminated.
    while ( ! glfwWindowShouldClose( window ) )
    {
        RenderScene();
        glfwSwapBuffers( window );
        glfwWaitEvents();
    }

    glfwTerminate();
    return 0;
}