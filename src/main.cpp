#include <iostream>
#include <math.h>

#include "implot.h"

#include "helper_funcs.h"
#include "fir_filter.hpp"

void key_callback(GLFWwindow* window, int key, int, int action, int){
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS){
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
}


int main(){
    GLFWwindow* window = glfw_makeNewWindow(1920, 720, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_callback);

    FIR_FILTER filt;
    float betaSlider = 0.5;

    while (!glfwWindowShouldClose(window)){
        glfwPollEvents();
        glfw_frame();

        filt.designTaps(4, betaSlider, 16);
        double* taps = filt.getTaps();
        size_t ntaps = filt.getNumTaps();


        ImGui::Begin("Data Grapher");
        ImGui::SliderFloat("Beta", &betaSlider, 0.1, 1);
        if(ImPlot::BeginPlot("Plots")){
            ImPlot::PlotBars("Taps", taps, ntaps);
            ImPlot::EndPlot();
        }
        ImGui::End();

        glfw_render(window);
    }

    glfw_cleanup(window);
    ImPlot::DestroyContext();

    return 0;
}