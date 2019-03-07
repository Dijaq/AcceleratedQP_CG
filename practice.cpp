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
//#include <igl/viewer/Viewer.h>
//#include "libs/viewer/Viewer.h"
//#include "class/smesh.h"

using namespace std;
using namespace Eigen;
using namespace glm;

GLuint VBO;

int main()
{
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

    do{
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
    		
		glfwSwapBuffers(window);
		glfwPollEvents();


    }
    while(glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) == 0);


	


    return 0;
}

/*static void RenderSceneCB()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();
}

static void InitializeGlutCallbacks()
{
    glutDisplayFunc(RenderSceneCB);
    //glutIdleFunc(RenderSceneCB);
}

int main(int argc, char** argv)
{

	Param_State mesh;

	read_mesh_2D("data_gecko/V.csv", "data_gecko/F.csv","data_gecko/eq_lhs.csv", "data_gecko/eq_rhs.csv", mesh);
    update_F(mesh.F);
    MatrixXd V3(mesh.V.rows(), mesh.V.cols()+1);
    create_column_zeros(mesh.V, V3);

    print_dimensions("V3: ", V3);

	vector<unsigned int> vertexIndices, uvIndices, normalIndice;
	vector<glm::vec3> temp_vertices;
	vector<glm::vec2> temp_uvs;
	vector<glm::vec3> temp_normasl;

	//Load vertex
    for(int i=0; i<V3.rows(); i++)
    {
		glm::vec3 vertex;
		vertex.x = V3(i,0);
		vertex.y = V3(i,1);
		vertex.z = V3(i,2);
    	

		temp_vertices.push_back(vertex);
    }

    //Load faces
    for(int i=0; i<mesh.F.rows(); i++)
    {
    	for(int j=0; j<mesh.F.cols(); j++)
    	{
    		vertexIndices.push_back(mesh.F(i,j));
    	}
    }

    vector<glm::vec3> vertices;
    for(int i=0; i<vertexIndices.size(); i++)
    {
    	int vertexIndex = vertexIndices[i];
    	glm::vec3 vertex = temp_vertices[vertexIndex];
    	vertices.push_back(vertex);
    }

    cout << "gl buffer data"	 << endl;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA);
    glutInitWindowSize(1024, 768);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Gecko");

    InitializeGlutCallbacks();
    // Must be done after glut is initialized!
    GLenum res = glewInit();
    if (res != GLEW_OK) {
      fprintf(stderr, "Error: '%s'\n", glewGetErrorString(res));
      return 1;
    }
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);


    glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), &vertices[0], GL_STATIC_DRAW);

    glutMainLoop();

    return 0;
}*/
