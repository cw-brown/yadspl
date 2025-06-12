#ifndef HELPERS_H
#define HELPERS_H

#include <GL/gl.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

GLFWwindow* glfw_makeNewWindow(int, int, const char*, bool = false, bool = true, bool = true);

void glfw_cleanup(GLFWwindow*);

void glfw_render(GLFWwindow*);

void glfw_frame();

#endif