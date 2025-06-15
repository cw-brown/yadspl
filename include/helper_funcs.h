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

/**
 * @brief Generate additive white guassian noise.
 * @param N0 The power of the noise relative to the input signal.
 * @return double 
 */
double AWGN(const double& N0, const double& bandwidth);

#endif