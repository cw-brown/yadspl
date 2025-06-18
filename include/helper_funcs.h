/**
 * @file helper_funcs.h
 * @author cw-brown (https://github.com/cw-brown)
 * @brief Helper function file for certain commonly used ideas. Mostly meant for ImGui stuff.
 * @version 0.1
 * @date 2025-06-18 
 */
#ifndef HELPERS_H
#define HELPERS_H

#include <random>

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