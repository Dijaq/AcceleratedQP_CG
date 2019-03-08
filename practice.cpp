#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "code/mathelpers.h"
#include "code/utils.h"
#include <glm/glm.hpp>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include "common/loadShader.cpp"
#include <glm/gtc/matrix_transform.hpp>
//#include <igl/viewer/Viewer.h>
//#include "libs/viewer/Viewer.h"
//#include "class/smesh.h"

using namespace std;
using namespace Eigen;
using namespace glm;

GLuint VBO;

void load_camera()
{
    // Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
    glm::mat4 Projection = glm::perspective(glm::radians(45.0f), (float) 1024 / (float)768, 0.1f, 100.0f);
      
    // Or, for an ortho camera :
    //glm::mat4 Projection = glm::ortho(-10.0f,10.0f,-10.0f,10.0f,0.0f,100.0f); // In world coordinates
      
    // Camera matrix
    glm::mat4 View = glm::lookAt(
        glm::vec3(4,3,3), // Camera is at (4,3,3), in World Space
        glm::vec3(0,0,0), // and looks at the origin
        glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
        );
      
    // Model matrix : an identity matrix (model will be at the origin)
    glm::mat4 Model = glm::mat4(1.0f);
    // Our ModelViewProjection : multiplication of our 3 matrices
    glm::mat4 mvp = Projection * View * Model; 
}


int main()
{
    cout << "load shaders" << endl;
    cout << "load shaders" << endl;
    if(!glfwInit())
    { 
        cout << "Error with GLFW" << endl;
        return -1;
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window;

    window = glfwCreateWindow(1024, 768, "Tutorial 01", NULL, NULL);

    if(window == NULL)
    {
        cout << "NO se puede abrir la ventana GLDW" << endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glewExperimental = true;

    if(glewInit() != GLEW_OK)
    {
        cout << "NO inicio GLEW" << endl;
        return -1;
    }

    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

    //Tut 2
    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    cout << "create vertex" << endl;
    static const GLfloat g_vertex_buffer_data[] = {
       -1.0f, -1.0f, 0.0f,
       1.0f, -1.0f, 0.0f,
       0.0f,  1.0f, 0.0f,
    };

    GLuint vertexbuffer;
    // Generar un buffer, poner el resultado en el vertexbuffer que acabamos de crear
    glGenBuffers(1, &vertexbuffer);
    // Los siguientes comandos le darán características especiales al 'vertexbuffer' 
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    // Darle nuestros vértices a  OpenGL.
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);
//Tut 2

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

     // Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
    glm::mat4 Projection = glm::perspective(glm::radians(45.0f), (float) 1024 / (float)768, 0.1f, 100.0f);
      
    // Or, for an ortho camera :
    //glm::mat4 Projection = glm::ortho(-10.0f,10.0f,-10.0f,10.0f,0.0f,100.0f); // In world coordinates
      
    // Camera matrix
    glm::mat4 View = glm::lookAt(
        glm::vec3(4,3,3), // Camera is at (4,3,3), in World Space
        glm::vec3(0,0,0), // and looks at the origin
        glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
        );
      
    // Model matrix : an identity matrix (model will be at the origin)
    glm::mat4 Model = glm::mat4(1.0f);
    // Our ModelViewProjection : multiplication of our 3 matrices
    glm::mat4 mvp = Projection * View * Model; 

    GLuint programID = LoadShaders("common/SimpleVertexShader.vertexshader", "common/SimpleFragmentShader.fragmentshader" );
    do{
        glUseProgram(programID);
        GLuint MatrixID = glGetUniformLocation(programID, "MVP");
  
// Send our transformation to the currently bound shader, in the "MVP" uniform
// This is done in the main loop since each model will have a different MVP matrix (At least for the M part)
        glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &mvp[0][0]);
        
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
        glVertexAttribPointer(
           0,                  // atributo 0. No hay razón particular para el 0, pero debe corresponder en el shader.
           3,                  // tamaño
           GL_FLOAT,           // tipo
           GL_FALSE,           // normalizado?
           0,                    // Paso
           (void*)0            // desfase del buffer
        );
        // Dibujar el triángulo !
        glDrawArrays(GL_TRIANGLES, 0, 3); // Empezar desde el vértice 0S; 3 vértices en total -> 1 triángulo
        glDisableVertexAttribArray(0);

        //load_camera();
            
        glfwSwapBuffers(window);
        glfwPollEvents();

    }
    while(glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) == 0);


    return 0;
}